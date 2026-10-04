from pathlib import Path
from dataclasses import dataclass
import hashlib,json,math,time,os,importlib.metadata
import numpy as np
import pandas as pd
import pymc as pm
import pytensor.tensor as pt
import arviz as az
import zipfile
from scipy.linalg import helmert
LAYERS=['Protein','PTM']
MODE='worker'
ASSAYS=['phospho','glyco']
MIN_PAIRED_N=102
MIN_PLEXES=16
FEATURE_SCOPE='original_signatures_plus_targets'
SENSITIVITY_TARGETS={'phospho':['NP_644805.1_Y705','NP_644805.1_S727'],'glyco_genes':['ASAH1']}
EXPECTED_SIGNATURE_ZIP_SHA256='8098755eeca1d564a4e5e004aa9bdfd0e16f187b703f21110fe517c45200eeab'
COMPATIBLE_LOG2_RATIO=False


def sha256_file(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for part in iter(lambda: handle.read(2 ** 20), b''):
            h.update(part)
    return h.hexdigest()

def canonical_hash(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, allow_nan=False, separators=(',', ':')).encode()).hexdigest()

def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.partial')
    temp.write_text(json.dumps(obj, indent=2, allow_nan=False, default=str))
    os.replace(temp, path)

def write_csv(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.partial')
    pd.DataFrame(data).to_csv(temp, index=False)
    os.replace(temp, path)

def truth(series):
    return series.astype(str).str.strip().str.lower().isin(['true', '1', 'yes'])

def unpack_inputs(path, work):
    if not path.is_file():
        raise FileNotFoundError(path)
    digest = sha256_file(path)
    dest = work / ('inputs_' + digest[:16])
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        manifests = [s for s in names if s.endswith('/package_manifest.json')]
        if len(manifests) != 1:
            raise ValueError('Expected one package manifest')
        prefix = manifests[0].removesuffix('package_manifest.json')
        manifest = json.loads(z.read(manifests[0]))
        for item in manifest:
            rel = item['path']
            target = (dest / rel).resolve()
            if not target.is_relative_to(dest.resolve()):
                raise ValueError('Unsafe archive path')
            content = z.read(prefix + rel)
            if len(content) != item['bytes'] or hashlib.sha256(content).hexdigest() != item['sha256']:
                raise ValueError('Input checksum failure: ' + rel)
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists() or sha256_file(target) != item['sha256']:
                target.write_bytes(content)
    return (dest, digest)

def load_signature_links(path, input_hash):
    if not path.is_file():
        if MODE == 'preflight':
            return (pd.DataFrame(), None, {'status': 'missing original-signature ZIP; full fits blocked'})
        raise FileNotFoundError(f'Required frozen signature archive is missing: {path}')
    if sha256_file(path) != EXPECTED_SIGNATURE_ZIP_SHA256:
        raise ValueError('Original signature ZIP changed')
    with zipfile.ZipFile(path) as z:
        payload = json.loads(z.read('payload.json'))
        if payload['input_zip_sha256'] != input_hash:
            raise ValueError('Cohort ZIP differs from archived original run')
        content = z.read(payload['signature_file'])
        if hashlib.sha256(content).hexdigest() != payload['signature_sha256']:
            raise ValueError('Signature checksum mismatch')
        signatures = json.loads(content)
    inv = pd.DataFrame(payload['original_features'])
    if len(inv) != 914 or not inv.global_feature.is_unique:
        raise ValueError('Expected 914 original labels')
    rows = []
    for f in inv.to_dict('records'):
        raw = signatures[f['CellType']][f['Signature']]
        if not isinstance(raw, list) or not all((isinstance(g, str) for g in raw)):
            raise ValueError('Unsigned original schema changed')
        genes = sorted(set((g.strip().upper() for g in raw if g.strip())))
        h = hashlib.sha256(('POS:' + ','.join(genes) + '|NEG:').encode()).hexdigest()
        if h != f['source_gene_hash']:
            raise ValueError('Source gene list hash mismatch')
        for gene in genes:
            rows.append({k: f[k] for k in ['global_feature', 'family', 'CellType', 'Signature', 'source_gene_hash']} | {'gene': gene})
    return (pd.DataFrame(rows), payload, {'status': 'verified', 'original_labels': 914})

def load_data(root, payload):
    p = root / 'prepared'
    cohort = pd.read_csv(p / 'model_covariates.csv')
    ids = cohort.base_sample_id.astype(str).tolist()
    ref = pd.read_csv(p / 'main_cohort_127.csv')
    assert len(ids) == len(set(ids)) == 127 and ids == ref.base_sample_id.tolist()
    assert np.isfinite(cohort[['bmi', 'age', 'sex_centered']]).all().all()
    assert cohort.bmi_group.value_counts().to_dict() == {'overweight': 57, 'normal_weight': 52, 'obese': 18}
    contract = json.loads((p / 'preparation_contract.json').read_text())
    assert sha256_file(root / 'provenance/analysis_cohort.csv') == contract['main_cohort_sha256']
    for col in ['adenosquamous_recorded', 'weight_loss_recorded']:
        cohort[col] = truth(cohort[col])
    if payload is not None:
        purity = pd.DataFrame(payload['dna_purity_records'])
        cohort = cohort.merge(purity, on='base_sample_id', how='left', validate='one_to_one', sort=False)
        for col in ['wes_purity', 'wgs_purity']:
            assert cohort[col].dropna().between(0, 1).all()
    mutation = pd.read_csv(p / 'matrices/mutation_tumor_main_cohort.tsv.gz', sep='\t', index_col=0, keep_default_na=False)
    nonsyn = {'missense_variant', 'frameshift_variant', 'stop_gained', 'stop_lost', 'start_lost', 'inframe_deletion', 'inframe_insertion', 'splice_donor_variant', 'splice_acceptor_variant'}
    import re
    for gene in ['TP53', 'CDKN2A', 'SMAD4']:
        calls = mutation.loc[gene, ids]
        if calls.eq('').any():
            raise ValueError('Missing mutation status')
        cohort['mutation_' + gene] = [int(bool(set(re.split('[,;]', s)) & nonsyn)) for s in calls]
    plex = pd.read_csv(p / 'tumor_plex_map.csv')
    assert not plex.duplicated(['layer', 'base_sample_id']).any()
    by_layer = {}
    for layer in ['protein'] + ASSAYS:
        s = plex[plex.layer.eq(layer)].set_index('base_sample_id').loc[ids]
        assert len(s) == 127 and s.plex.notna().all()
        by_layer[layer] = s.plex.to_numpy()
    for layer in ASSAYS:
        a = by_layer[layer]
        b = by_layer['protein']
        if not np.array_equal(a[:, None] == a[None, :], b[:, None] == b[None, :]):
            raise ValueError('Assay batch partitions differ: needs an explicitly separate-batch model')
    matrices = {}
    annotations = {}
    for layer in ['protein'] + ASSAYS:
        with np.load(p / f'matrices/{layer}_tumor_main_cohort.npz', allow_pickle=False) as z:
            assert z['patient_ids'].tolist() == ids
            values = z['values'].astype(float)
            features = z['feature_ids'].astype(str)
        assert values.shape == (len(features), 127) and len(set(features)) == len(features)
        assert not np.isinf(values).any()
        ann = pd.read_csv(p / f'feature_metadata/{layer}_annotations.csv', keep_default_na=False)
        assert ann.feature_id.tolist() == features.tolist()
        matrices[layer] = {'values': values, 'feature_ids': features}
        annotations[layer] = ann
    return dict(cohort=cohort, matrices=matrices, annotations=annotations, batches=by_layer, contract=contract)

def make_design(cohort, mask, adjustment='mundlak', extra=None, categorical=False):
    d = cohort.loc[mask].copy().reset_index(drop=True)
    levels = sorted(d.plex.unique())
    group = pd.Categorical(d.plex, categories=levels).codes
    raw = {'overweight': d.bmi_group.eq('overweight').to_numpy(float), 'obese': d.bmi_group.eq('obese').to_numpy(float)} if categorical else {'bmi_5': d.bmi.to_numpy() / 5}
    focal = list(raw)
    raw.update(age_10y=d.age.to_numpy() / 10, sex=d.sex_centered.to_numpy())
    if extra == 'stage':
        for stage in ['II', 'III', 'IV']:
            raw['stage_' + stage] = d.stage_major.eq(stage).to_numpy(float)
    elif extra == 'mutations':
        for gene in ['TP53', 'CDKN2A', 'SMAD4']:
            raw['mutation_' + gene] = d['mutation_' + gene].to_numpy(float)
    elif extra == 'wes_purity':
        raw['wes_purity_10pct'] = d.wes_purity.to_numpy() / 0.1
    names = []
    columns = []
    between = []
    dropped = []
    for key, x in raw.items():
        if not np.isfinite(x).all():
            raise ValueError('Non-finite covariate ' + key)
        mean = pd.Series(x).groupby(group).transform('mean').to_numpy()
        v = x - mean
        if v.std() < 1e-10:
            if key in focal:
                raise ValueError('No within-plex variation for ' + key)
            dropped.append(key + '_within')
        else:
            names.append(key + '_within')
            columns.append(v)
        if adjustment == 'mundlak':
            between.append((key + '_between', mean - x.mean()))
    for name, value in between:
        if value.std() < 1e-10:
            dropped.append(name)
        else:
            names.append(name)
            columns.append(value)
    X = np.column_stack(columns)
    Q = helmert(len(levels), full=False).T if adjustment == 'fixed' else None
    full = np.column_stack([np.ones(len(d)), X] + ([Q[group]] if Q is not None else []))
    if np.linalg.matrix_rank(full) != full.shape[1]:
        raise ValueError('Rank-deficient design; no confounder deletion')
    return dict(data=d, X=X, names=names, group=group, levels=levels, Q=Q, adjustment=adjustment, extra=extra, dropped_constant_terms=dropped, focal=[s + '_within' for s in focal], categorical=categorical)

def prepare_features(data, links, qc_path=None):
    protein = data['matrices']['protein']
    lookup = {g: i for i, g in enumerate(protein['feature_ids'])}
    genes = set(links.gene) if len(links) else set()
    audits = []
    records = []
    arrays = {}
    qc = None
    if qc_path:
        qc = pd.read_csv(qc_path, keep_default_na=False)
        if not {'assay', 'feature_id', 'pass_qc', 'provenance'} <= set(qc):
            raise ValueError('Invalid QC schema')
        if qc.duplicated(['assay', 'feature_id']).any():
            raise ValueError('Duplicate QC feature')
        if qc.provenance.eq('').any():
            raise ValueError('QC provenance missing')
        qc = {(r.assay, r.feature_id): str(r.pass_qc).lower() in ('true', '1', 'yes') for r in qc.itertuples()}
    for assay in ASSAYS:
        matrix = data['matrices'][assay]
        ann = data['annotations'][assay]
        for row in ann.itertuples(index=True):
            sid = row.feature_id
            gene = row.Gene.strip().upper()
            reason = ''
            n = 0
            nplex = 0
            target = sid in SENSITIVITY_TARGETS['phospho'] if assay == 'phospho' else gene in SENSITIVITY_TARGETS['glyco_genes']
            if FEATURE_SCOPE == 'original_signatures_plus_targets' and gene not in genes and (not target):
                reason = 'outside original genes and named targets' if len(links) else 'signature archive unavailable; scope unresolved'
            elif qc is not None and (not qc.get((assay, sid), False)):
                reason = 'external QC not passed or absent'
            elif gene not in lookup:
                reason = 'parent protein unavailable by exact gene symbol'
            else:
                pair = np.column_stack([protein['values'][lookup[gene]], matrix['values'][row.Index]])
                mask = np.isfinite(pair).all(1)
                n = int(mask.sum())
                nplex = len(set(data['batches'][assay][mask]))
                if n < MIN_PAIRED_N or nplex < MIN_PLEXES:
                    reason = 'paired coverage below threshold'
                elif np.any(pair[mask].std(0, ddof=1) < 1e-08):
                    reason = 'constant paired response'
            rec = dict(assay=assay, feature_id=sid, gene=gene, n_paired=n, n_plexes=nplex, target=target, included=not bool(reason), exclusion_reason=reason, localization_policy='external_QC' if qc is not None else 'source_processed_unverified_localization')
            audits.append(rec)
            if reason:
                continue
            mean = pair[mask].mean(0)
            sd = pair[mask].std(0, ddof=1)
            Y = (pair - mean) / sd
            key = assay + '_' + hashlib.sha256(sid.encode()).hexdigest()[:20]
            arrays[key] = {'Y': Y, 'mask': mask, 'means': mean, 'sds': sd}
            records.append(rec | dict(key=key, protein_row=gene, protein_mean=float(mean[0]), protein_sd=float(sd[0]), ptm_mean=float(mean[1]), ptm_sd=float(sd[1]), parent_resolution='gene_level'))
    return (pd.DataFrame(records), pd.DataFrame(audits), arrays)

def feature_cohort(data, assay):
    d = data['cohort'].copy()
    d['plex'] = data['batches'][assay]
    return d

def variant_spec(d, base, variant):
    mask = base.copy()
    extra = None
    adjustment = 'mundlak'
    categorical = False
    if variant == 'fixed_plex':
        adjustment = 'fixed'
    elif variant in ('stage_subset_base', 'stage_adjusted'):
        mask &= d.stage_major.isin(['I', 'II', 'III', 'IV']).to_numpy()
        if variant == 'stage_adjusted':
            extra = 'stage'
    elif variant in ('wes_subset_base', 'wes_adjusted'):
        if 'wes_purity' not in d:
            raise ValueError('WES purity absent from original payload')
        mask &= np.isfinite(d.wes_purity.to_numpy())
        if variant == 'wes_adjusted':
            extra = 'wes_purity'
    elif variant == 'mutation_adjusted':
        extra = 'mutations'
    elif variant == 'exclude_recorded_weight_loss':
        mask &= ~d.weight_loss_recorded.to_numpy(bool)
    elif variant == 'exclude_adenosquamous':
        mask &= ~d.adenosquamous_recorded.to_numpy(bool)
    elif variant == 'categorical':
        categorical = True
    elif variant != 'primary':
        raise ValueError(variant)
    if mask.sum() < 30:
        raise ValueError('Fewer than 30 patients in sensitivity subset')
    return (mask, adjustment, extra, categorical)

def build_model(design, observed, config):
    import pymc as pm
    import pytensor.tensor as pt
    coords = {'patient': design['data'].base_sample_id.tolist(), 'layer': LAYERS, 'term': design['names'], 'plex': design['levels']}
    if design['adjustment'] == 'fixed':
        coords['plex_contrast'] = list(range(len(design['levels']) - 1))
    prior = np.full(len(design['names']), config['COVARIATE_PRIOR_SD'])
    for i, name in enumerate(design['names']):
        if name.startswith('bmi_5_'):
            prior[i] = config['BMI_PRIOR_SD']
    with pm.Model(coords=coords) as model:
        alpha = pm.Normal('alpha', 0, 0.5, dims='layer')
        coefficient = pm.Normal('coefficient', 0, prior[:, None], dims=('term', 'layer'))
        if design['adjustment'] == 'mundlak':
            group_scale = pm.HalfNormal('group_scale', config['GROUP_SD_PRIOR'], dims='layer')
            group_z = pm.Normal('group_z', 0, 1, dims=('plex', 'layer'))
            group_effect = pm.Deterministic('group_effect', group_z * group_scale, dims=('plex', 'layer'))
        else:
            group_coefficient = pm.Normal('group_coefficient', 0, config['FIXED_GROUP_PRIOR_SD'], dims=('plex_contrast', 'layer'))
            group_effect = pm.Deterministic('group_effect', pt.dot(design['Q'], group_coefficient), dims=('plex', 'layer'))
        mu = alpha + pt.dot(design['X'], coefficient) + group_effect[design['group']]
        chol, corr, scales = pm.LKJCholeskyCov('resid_chol', n=2, eta=2, sd_dist=pm.HalfNormal.dist(1.0, shape=2), compute_corr=True)
        nu_minus_two = pm.Exponential('nu_minus_two', 0.1)
        nu = pm.Deterministic('nu', 2 + nu_minus_two)
        positions = [design['names'].index(s) for s in design['focal']]
        slopes = coefficient[positions]
        coupling = pm.Deterministic('protein_to_ptm_coupling', corr[0, 1] * scales[1] / scales[0])
        pm.Deterministic('focal_effect', slopes)
        pm.Deterministic('ptm_minus_protein', slopes[:, 1] - slopes[:, 0])
        pm.Deterministic('ptm_conditional_effect', slopes[:, 1] - coupling * slopes[:, 0])
        pm.Deterministic('residual_correlation', corr[0, 1])
        pm.MvStudentT('paired_score', nu=nu, mu=mu, chol=chol, observed=observed, dims=('patient', 'layer'))
    return model

def posterior_array(idata, name):
    da = idata.posterior[name].stack(sample=('chain', 'draw'))
    return da.transpose('sample', *[d for d in da.dims if d != 'sample']).values

def hdi(values, probability=0.95):
    a = np.sort(np.asarray(values, float).ravel())
    if not np.isfinite(a).all() or len(a) < 20:
        raise ValueError('Cannot summarize incomplete posterior draws')
    width = int(np.floor(probability * len(a)))
    left = int(np.argmin(a[width:] - a[:-width]))
    return (float(a[left]), float(a[left + width]))

def numerical_diagnostics(idata, model, config, folder):
    import arviz as az
    names = [v.name for v in model.free_RVs] + ['focal_effect', 'ptm_conditional_effect', 'protein_to_ptm_coupling']
    summary = az.summary(idata, var_names=names, kind='diagnostics', round_to='none')
    summary.rename_axis('parameter').reset_index().to_csv(folder / 'parameter_diagnostics.csv', index=False)
    columns = summary[['r_hat', 'ess_bulk', 'ess_tail']].to_numpy(float)
    finite = bool(np.isfinite(columns).all())
    rhat = float(np.nanmax(summary.r_hat))
    bulk, tail = (float(np.nanmin(summary.ess_bulk)), float(np.nanmin(summary.ess_tail)))
    stats = idata.sample_stats
    divergences = int(stats['diverging'].sum())
    bfmi = float(np.nanmin(az.bfmi(idata)))
    if 'reached_max_treedepth' in stats:
        depth_hits = int(stats['reached_max_treedepth'].sum())
    elif 'tree_depth' in stats:
        depth_hits = int((stats.tree_depth >= config['MAX_TREEDEPTH']).sum())
    elif 'depth' in stats:
        depth_hits = int((stats.depth >= config['MAX_TREEDEPTH']).sum())
    else:
        raise RuntimeError('Sampler returned no tree-depth statistic; cannot verify the gate')
    passed = finite and rhat <= 1.01 and (min(bulk, tail) >= 400) and (bfmi >= 0.3) and (divergences == 0) and (depth_hits == 0)
    return {'max_rhat': rhat if np.isfinite(rhat) else None, 'min_ess_bulk': bulk if np.isfinite(bulk) else None, 'min_ess_tail': tail if np.isfinite(tail) else None, 'min_bfmi': bfmi if np.isfinite(bfmi) else None, 'divergences': divergences, 'maximum_depth_hits': depth_hits, 'all_diagnostics_finite': finite, 'numerical_checks_pass': bool(passed), 'publication_profile': config['PROFILE'] == 'publication'}

def effect_summary(idata, feature, variant, design, passed):
    beta = posterior_array(idata, 'focal_effect')
    conditional = posterior_array(idata, 'ptm_conditional_effect')
    rows = []

    def add(quantity, draws, term, units):
        lo, hi = hdi(draws)
        rows.append(dict(key=feature['key'], assay=feature['assay'], feature_id=feature['feature_id'], gene=feature['gene'], variant=variant, term=term, quantity=quantity, posterior_mean=float(np.mean(draws)), hdi_95_lower=lo, hdi_95_upper=hi, hdi_excludes_zero=bool(lo > 0 or hi < 0), probability_positive=float(np.mean(draws > 0)), n_patients=len(design['data']), n_plexes=len(design['levels']), numerical_checks_pass=bool(passed), effect_units=units, interval_type='pointwise 95% HDI; no multiplicity control'))
    for j, term in enumerate(design['focal']):
        units = 'primary-paired-sample response SD' + (' vs normal weight, within plex' if design['categorical'] else ' per 5 kg/m2, within plex')
        for label, values in [('Protein', beta[:, j, 0]), ('PTM', beta[:, j, 1]), ('PTM_minus_Protein_standardized', beta[:, j, 1] - beta[:, j, 0]), ('PTM_given_Protein', conditional[:, j])]:
            add(label, values, term, units)
        if COMPATIBLE_LOG2_RATIO:
            add('log2_PTM_minus_Protein', feature['ptm_sd'] * beta[:, j, 1] - feature['protein_sd'] * beta[:, j, 0], term, 'log2 ratio' + (' vs normal weight' if design['categorical'] else ' per 5 kg/m2'))
    add('residual_Protein_PTM_correlation', posterior_array(idata, 'residual_correlation'), 'residual', 'correlation')
    add('protein_to_ptm_coupling', posterior_array(idata, 'protein_to_ptm_coupling'), 'residual', 'PTM SD per protein SD')
    return pd.DataFrame(rows)

def compute_ppc(idata, design, observed, config, folder, seed):
    rng = np.random.default_rng(seed + 231)
    a, c, g = [posterior_array(idata, n) for n in ['alpha', 'coefficient', 'group_effect']]
    total = len(a)
    chosen = np.sort(rng.choice(total, min(config['PPC_DRAWS'], total), replace=False))
    mu = a[chosen, None, :] + np.einsum('nt,stl->snl', design['X'], c[chosen]) + g[chosen][:, design['group'], :]
    scales = posterior_array(idata, 'resid_chol_stds')[chosen]
    corr = posterior_array(idata, 'resid_chol_corr')[chosen]
    nu = posterior_array(idata, 'nu')[chosen]
    S = corr * scales[:, :, None] * scales[:, None, :]
    L = np.linalg.cholesky(S)
    z = rng.normal(size=mu.shape)
    gaussian = np.einsum('sij,snj->sni', L, z)
    w = np.sqrt(rng.chisquare(nu[:, None], size=mu.shape[:2]) / nu[:, None])
    replicate = mu + gaussian / w[:, :, None]
    observed_residual = observed[None, :, :] - mu
    replicate_residual = replicate - mu
    checks = []

    def add(name, obs, rep):
        upper = float(np.mean(rep >= obs))
        checks.append({'check': name, 'observed_draw_median': float(np.median(obs)), 'replicated_median': float(np.median(rep)), 'replicated_95_lower': float(np.quantile(rep, 0.025)), 'replicated_95_upper': float(np.quantile(rep, 0.975)), 'posterior_predictive_upper_tail_probability': upper, 'flag': bool(upper < 0.025 or upper > 0.975), 'interpretation': 'Descriptive predictive check; not a calibrated p-value'})
    for j, layer in enumerate(LAYERS):
        add(layer + '_residual_SD', observed_residual[:, :, j].std(1, ddof=1), replicate_residual[:, :, j].std(1, ddof=1))
        add(layer + '_maximum_absolute_residual', np.abs(observed_residual[:, :, j]).max(1), np.abs(replicate_residual[:, :, j]).max(1))

    def corrs(x):
        v = x - x.mean(1, keepdims=True)
        return np.sum(v[:, :, 0] * v[:, :, 1], axis=1) / np.sqrt(np.sum(v[:, :, 0] ** 2, axis=1) * np.sum(v[:, :, 1] ** 2, axis=1))
    add('Protein_PTM_residual_correlation', corrs(observed_residual), corrs(replicate_residual))
    low, high = np.quantile(replicate, [0.05, 0.95], axis=0)
    for j, layer in enumerate(LAYERS):
        checks.append({'check': layer + '_90pct_in_sample_predictive_coverage', 'observed_draw_median': float(np.mean((observed[:, j] >= low[:, j]) & (observed[:, j] <= high[:, j]))), 'interpretation': 'In-sample descriptive coverage; not out-of-sample validation', 'flag': False})
    write_csv(folder / 'predictive_checks.csv', checks)
    write_csv(folder / 'patient_residuals.csv', pd.DataFrame({'base_sample_id': design['data'].base_sample_id, 'Protein_residual': observed[:, 0] - mu[:, :, 0].mean(0), 'PTM_residual': observed[:, 1] - mu[:, :, 1].mean(0)}))
    return int(sum((bool(x.get('flag', False)) for x in checks)))
