from datetime import datetime, timezone
BASE = Path(cfg['base'])
RI_ROOT = BASE / 'normalized_state_bhm_ri_hc_v1'
N_COMPOSITE_PILOT = 500

def ols_hc3(X, y, n_total, ix):
    inv = np.linalg.inv(X.T @ X); A = inv @ X.T
    h = np.einsum('ij,ji->i', X, A)
    B = A @ y; r = y - X @ B
    u = np.zeros(n_total); u[ix] = A[1] * r / (1 - h)
    return B, r, h, u

comp_rows, comp_U, comp_designs = [], [], []
for family in FAMILIES:
    a = arrays[family]; n = len(patients)
    Y = np.full((n, len(a['features'])), np.nan); Y[a['pi'], a['fi']] = a['y']
    scale = 5 / metadata[family]['bmi_sd']
    for c, ct in enumerate(a['celltypes']):
        js = np.flatnonzero(a['f2c'] == c)
        Yc = Y[:, js]
        obs = np.isfinite(Yc).all(axis=1)
        assert (np.isfinite(Yc).any(axis=1) == obs).all(), f'{ct}: partial profiles'
        ix = np.flatnonzero(obs); X = a['X'][ix]
        y = Yc[ix].mean(axis=1)
        B, r, h, u = ols_hc3(X, y, n, ix)
        se = np.sqrt((u ** 2).sum()); df = len(ix) - X.shape[1]
        crit = stats.t.ppf(.975, df)
        C = np.corrcoef(Yc[ix].T); m = len(js)
        rbar = float(C[np.triu_indices(m, 1)].mean())
        comp_rows.append(dict(family=family, celltype=ct, n_patients=len(ix), n_signatures=m,
                              mean_pairwise_r=rbar, design_effect_se_ratio=np.sqrt(1 + (m - 1) * rbar),
                              estimate=B[1] * scale, se_hc3=se * scale,
                              pointwise_low=(B[1] - crit * se) * scale, pointwise_high=(B[1] + crit * se) * scale,
                              p_hc3=2 * stats.t.sf(abs(B[1] / se), df), df=df))
        comp_U.append(u * scale)
        comp_designs.append(dict(ix=ix, X=X, B=B, r=r, h=h, df=df))
composites = pd.DataFrame(comp_rows)
Uc = np.column_stack(comp_U)
maxc = multiplier_max(Uc, composites.se_hc3.to_numpy(), N_MULTIPLIERS, SEED + 5)
crit_c = float(np.quantile(maxc, .95, method='higher'))
composites['simultaneous_low'] = composites.estimate - crit_c * composites.se_hc3
composites['simultaneous_high'] = composites.estimate + crit_c * composites.se_hc3
composites['p_maxT'] = [(1 + np.sum(maxc >= abs(b / s))) / (len(maxc) + 1)
                        for b, s in zip(composites.estimate, composites.se_hc3)]
composites['pointwise_excludes_zero'] = (composites.pointwise_low > 0) | (composites.pointwise_high < 0)
composites['simultaneous_excludes_zero'] = (composites.simultaneous_low > 0) | (composites.simultaneous_high < 0)
print('Composite simultaneous critical value:', round(crit_c, 3))

rng = np.random.default_rng(SEED + 6); pil = []
truth = np.array([d['B'][1] for d in comp_designs])
for rep in tqdm(range(N_COMPOSITE_PILOT), desc='Composite pilot'):
    signs = rng.choice([-1., 1.], size=len(patients))
    bb, ss, uu = [], [], []
    for d in comp_designs:
        y = d['X'] @ d['B'] + signs[d['ix']] * d['r'] / np.sqrt(1 - d['h'])
        B, _, _, u = ols_hc3(d['X'], y, len(patients), d['ix'])
        bb.append(B[1]); ss.append(np.sqrt((u ** 2).sum())); uu.append(u)
    bb, ss, uu = np.array(bb), np.array(ss), np.column_stack(uu)
    mm = multiplier_max(uu, ss, PILOT_MULTIPLIERS, SEED + 2000 + rep)
    t = np.abs((bb - truth) / ss)
    dfs = composites.df.to_numpy()
    pil.append(dict(replicate=rep, pointwise_coverage=float(np.mean(t <= stats.t.ppf(.975, dfs))),
                    simultaneous_covered=bool(np.all(t <= np.quantile(mm, .95, method='higher')))))
pil = pd.DataFrame(pil)
print('Composite pilot pointwise coverage:', round(pil.pointwise_coverage.mean(), 3),
      '| simultaneous:', round(pil.simultaneous_covered.mean(), 3))

bayes, feat_cmp = [], []
for family in FAMILIES:
    prim = f"{e.SPECS[family]['architecture']}__{family}__lowread1000"
    for label, suffix in (('ri_hc_gaussian', '__ri_hc'), ('ri_hc_t', '__ri_hc_t')):
        p = RI_ROOT / f'{prim}{suffix}' / 'comparison'
        if not (p / 'celltype_comparison.csv').is_file():
            print('missing', label, family); continue
        ct = pd.read_csv(p / 'celltype_comparison.csv')
        if label == 'ri_hc_gaussian':
            bayes.append(ct[['celltype', 'prim_mean', 'prim_hdi_low', 'prim_hdi_high']]
                         .set_axis(['celltype', 'mean', 'low', 'high'], axis=1).assign(model='published'))
        bayes.append(ct[['celltype', 'new_mean', 'new_hdi_low', 'new_hdi_high']]
                     .set_axis(['celltype', 'mean', 'low', 'high'], axis=1).assign(model=label))
        if label == 'ri_hc_t':
            fc = pd.read_csv(p / 'feature_comparison.csv')
            m_ = results[results.family.eq(family)].merge(fc, on='global_feature', validate='one_to_one')
            feat_cmp.append(pd.DataFrame(dict(family=family, global_feature=m_.global_feature,
                hc3_estimate=m_.estimate, hc3_width=m_.pointwise_high - m_.pointwise_low,
                hc3_excl=m_.pointwise_excludes_zero, published_width=m_.prim_hdi_high - m_.prim_hdi_low,
                published_excl=m_.prim_excl, t_mean=m_.new_mean, t_width=m_.new_hdi_high - m_.new_hdi_low,
                t_excl=m_.new_excl)))
bayes = pd.concat(bayes, ignore_index=True)
bayes['width'] = bayes.high - bayes.low
wide = bayes.pivot_table(index='celltype', columns='model', values=['mean', 'width'])
wide.columns = [f'{m}_{s}' for s, m in wide.columns]
bench = composites.merge(wide.reset_index(), on='celltype', how='left')
bench['composite_width'] = bench.pointwise_high - bench.pointwise_low
feat_cmp = pd.concat(feat_cmp, ignore_index=True) if feat_cmp else pd.DataFrame()

pd.set_option('display.width', 250)
cols = ['family', 'celltype', 'n_signatures', 'mean_pairwise_r', 'estimate', 'pointwise_low', 'pointwise_high',
        'p_maxT', 'composite_width', 'published_width', 'ri_hc_gaussian_width', 'ri_hc_t_width',
        'published_mean', 'ri_hc_gaussian_mean', 'ri_hc_t_mean']
print(bench[[c for c in cols if c in bench]].round(3).to_string(index=False))
print('\nComposites excluding zero: pointwise', int(bench.pointwise_excludes_zero.sum()),
      '| simultaneous', int(bench.simultaneous_excludes_zero.sum()), 'of', len(bench))
if len(feat_cmp):
    print('\nFeature intervals vs unpooled HC3 (median width ratio; pointwise exclusions)')
    print(feat_cmp.groupby('family').apply(lambda g: pd.Series({
        'published/HC3': (g.published_width / g.hc3_width).median(), 'ri_hc_t/HC3': (g.t_width / g.hc3_width).median(),
        'corr(HC3, ri_hc_t means)': np.corrcoef(g.hc3_estimate, g.t_mean)[0, 1],
        'HC3 excl': int(g.hc3_excl.sum()), 'published excl': int(g.published_excl.sum()), 'ri_hc_t excl': int(g.t_excl.sum())}),
        include_groups=False).round(3).to_string())

stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
OUT7 = SOURCE.parent / 'rna_celltype_composite_benchmark_v1' / stamp
OUT7.mkdir(parents=True, exist_ok=False)
bench.to_csv(OUT7 / 'celltype_composites_vs_bayes.csv', index=False)
pil.to_csv(OUT7 / 'composite_pilot.csv', index=False)
feat_cmp.to_csv(OUT7 / 'feature_intervals_vs_hc3.csv', index=False)
(OUT7 / 'README.txt').write_text(
    'Equal-weight composite of standardised signature scores per patient and cell type (complete profiles only), '
    'regressed on BMI + age + sex. HC3 pointwise intervals; shared-patient Gaussian multiplier simultaneous intervals '
    'across all composites. Per 5 kg/m2. The composite slope is a patient-level benchmark, not the hierarchical m_c. '
    f'Seed {SEED + 5}; pilot {N_COMPOSITE_PILOT} sign-flip datasets.\n')
archive7 = shutil.make_archive(str(OUT7.parent / ('celltype_composite_benchmark_' + stamp)), 'zip', OUT7)
print('Saved:', OUT7)
from google.colab import files
files.download(archive7)
