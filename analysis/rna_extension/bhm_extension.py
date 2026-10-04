"""Hierarchical RNA models and diagnostics for read-dependent residual-scale analysis."""
from pathlib import Path
import hashlib, json, math, platform, time, traceback, importlib.metadata
from dataclasses import dataclass, replace
import numpy as np
import pandas as pd
import scipy.linalg as la
from scipy.special import logsumexp
from scipy.stats import t as student_t
from numpy.polynomial.hermite import hermgauss
import pymc as pm
import pytensor.tensor as pt
import arviz as az
import xarray as xr

VERSION = 'pdac-extension-v1.1-autorun'
PINS = {'pymc':'5.28.5', 'arviz':'0.23.4', 'nutpie':'0.16.11'}
HASHES = {
 'pca_manifest.json':'a4263b46efb1c0d1d98a01b127ddb143bb7616ce88f8b364959613bbb74dcd30',
 'component_inclusion_tiers.csv':'16494a18ed467223134f0aed2efbf63c084d2f6b09332ca9854a9c403dd67294',
 'pca_input_hashes.csv':'6e2eaedbb5efe7c00d82e82c939d4416a5a3853f14bde881a55471cf874259c9',
 'combined_review_manifest.json':'04b0b1d788fc69e78ded4b44e050d487bfe231dfd5e5bd2b34f6b359c74f4476'}
SPECS = {
 'immune_coarse':dict(compartment='immune_course',tier='tiered_primary_keep',nu='celltype',features=306,celltypes=8,dropped=25,architecture='v2p2_g1_fixed_ctnu',trace_hash='e19078ccd56c5fba623baadcaf86f124e168b73a3950f40906c1e7861015202c'),
 'immune_fine':dict(compartment='immune_fine',tier='tiered_primary_keep',nu='global',features=496,celltypes=14,dropped=58,architecture='v2p1c_g1_fixed_gnu',trace_hash='3cb220f9863fecbc7ee9bce1f01d241d593498e9b8b1767865f6508010f43e98'),
 'nonimmune':dict(compartment='non_immune',tier='strict_qc_keep',nu='global',features=112,celltypes=4,dropped=2,architecture='v2p1c_g1_fixed_gnu',trace_hash='b030acc10bac4f68d2d3dd183b0570adacbcfeaf972e1ab9926736fe6588b3c4')}
PROFILES = {
 'test':dict(chains=2,tune=60,draws=60,target_accept=.90),
 'smoke':dict(chains=2,tune=500,draws=500,target_accept=.97),
 'publication':dict(chains=4,tune=2000,draws=4000,target_accept=.99)}
PRIORS = dict(baseline=.50,celltype_bmi=.20,feature_bmi=.30,celltype_covariate=.25,
 feature_covariate=.25,patient=.50,log_sigma_mean=math.log(.65),log_sigma_sd=.40,
 log_sigma_celltype_sd=.25,log_sigma_feature_sd=.50,nu_rate=.10,
 precision_slope_sd=.25,factor_loading_sd=.35)

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def digest(obj):
    return hashlib.sha256(json.dumps(obj,sort_keys=True,default=str).encode()).hexdigest()

def write_json(obj,path):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,sort_keys=True,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
    tmp.replace(path)

def csv(df,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp');df.to_csv(tmp,index=False);tmp.replace(path)

def require(cond,message):
    if not cond: raise ValueError(message)

def environment():
    return dict(python=platform.python_version(),packages={p:importlib.metadata.version(p)
        for p in ['pymc','arviz','nutpie','pytensor','numpy','pandas','scipy','xarray','h5netcdf']})

def check_environment():
    for p,v in PINS.items(): require(importlib.metadata.version(p)==v,f'{p} must be {v}; restart after installing.')

def check_hash(path,expected):
    require(Path(path).is_file(),f'Missing input: {path}')
    require(sha(path)==expected,f'Input hash mismatch: {path}. Do not bypass; reconcile input version.')

def strict_bool(s):
    out=s.map({True:True,False:False,'True':True,'False':False,'true':True,'false':False})
    require(out.notna().all(),'Invalid inclusion flag');return out.astype(bool)

def load_frozen(base,family,pca_override=None,review_override=None,run_override=None):
    """Reconstruct exactly the published low-read scoring endpoint; do not change source files."""
    base=Path(base); spec=SPECS[family]
    pca=Path(pca_override) if pca_override else base/'bayesprism_initial_issue_verification_v1/publication_pca_tiered_qc_v1/20260906T013040524557Z'
    for n in ['pca_manifest.json','component_inclusion_tiers.csv','pca_input_hashes.csv']:check_hash(pca/n,HASHES[n])
    pmf=json.loads((pca/'pca_manifest.json').read_text())
    review=Path(review_override or pmf['combined_review_directory'])
    check_hash(review/'combined_review_manifest.json',HASHES['combined_review_manifest.json'])
    rm=json.loads((review/'combined_review_manifest.json').read_text())
    run=Path(run_override or rm['selected_runs'][spec['compartment']])
    locked=pd.read_csv(pca/'pca_input_hashes.csv')
    for kind,n in [('cohort','analysis_cohort.csv'),('scores','audit_scores_all_variants.csv.gz')]:
        h=locked.loc[(locked.compartment==spec['compartment'])&(locked.input==kind),'sha256']
        require(len(h)==1,f'Expected one hash for {kind}');check_hash(run/n,h.iloc[0])
    cohort=pd.read_csv(run/'analysis_cohort.csv',dtype={'base_sample_id':str})
    scores=pd.read_csv(run/'audit_scores_all_variants.csv.gz',dtype={'base_sample_id':str})
    qc=pd.read_csv(run/'patient_qc.csv',dtype={'base_sample_id':str})
    tiers=pd.read_csv(pca/'component_inclusion_tiers.csv')
    require(len(cohort)==127 and cohort.base_sample_id.nunique()==127,'Expected 127 unique patients')
    require(cohort.bmi_group.value_counts().to_dict()=={'normal_weight':52,'overweight':57,'obese':18},'BMI groups changed')
    cols=['base_sample_id','bmi','bmi_z','bmi_group','age_z','sex_centered']
    require(set(cols)<=set(cohort),'Cohort schema mismatch')
    require(set(cohort.sex_centered)=={-.5,.5},'Sex coding must remain -0.5/+0.5')
    tt=tiers[tiers.compartment==spec['compartment']]
    keep=tt.loc[strict_bool(tt[spec['tier']]),'CellType'].astype(str)
    f=scores[(scores.variant=='within_component_logcpm_z_cap5')&scores.CellType.isin(keep)].copy()
    require(set(['feature','Signature','score'])<=set(f),'Score schema mismatch')
    f['global_feature']=spec['compartment']+'||'+f.feature.astype(str)
    f['global_celltype']=spec['compartment']+'||'+f.CellType.astype(str)
    f=f.merge(cohort[cols],on='base_sample_id',validate='many_to_one')
    require(not f.duplicated(['base_sample_id','global_feature']).any(),'Duplicate score rows')
    require(f.groupby('global_feature').base_sample_id.nunique().eq(127).all(),'Incomplete full grid')
    require((f.global_feature.nunique(),f.global_celltype.nunique())==(spec['features'],spec['celltypes']),'Feature inventory changed')
    f['score']=pd.to_numeric(f.score,errors='raise');g=f.groupby('global_feature').score
    f['score_model']=(f.score-g.transform('mean'))/g.transform('std')
    f['bmi_standardized']=f.bmi_z
    f=f.merge(qc[['base_sample_id','CellType','component_total','theta']],on=['base_sample_id','CellType'],validate='many_to_one',how='left')
    require(np.isfinite(f[['score_model','bmi','bmi_standardized','age_z','sex_centered','component_total','theta']].to_numpy(float)).all(),'Nonfinite inputs')
    require((f.component_total>=0).all() and f.theta.between(0,1).all(),'Invalid reads/theta')
    f['low_read']=f.component_total<1000
    dropped=f.loc[f.low_read,['base_sample_id','CellType','component_total']].drop_duplicates()
    require(len(dropped)==spec['dropped'],'Low-read counts differ from frozen analysis')
    require(f.loc[~f.low_read].groupby('global_feature').base_sample_id.nunique().min()>=30,'Too few observations')
    files=[pca/n for n in ['pca_manifest.json','component_inclusion_tiers.csv','pca_input_hashes.csv']]+[review/'combined_review_manifest.json']+[run/n for n in ['analysis_cohort.csv','audit_scores_all_variants.csv.gz','patient_qc.csv']]
    meta=dict(family=family,source_hashes={str(p):sha(p) for p in files},source_run=str(run),
        patients=127,features=spec['features'],dropped_profiles=len(dropped),original_scaling='all 127 before removal',
        bmi_sd=float(cohort.bmi.std(ddof=1)),post_hoc=True)
    expected=(cohort.bmi-cohort.bmi.mean())/cohort.bmi.std(ddof=1)
    require(np.allclose(expected,cohort.bmi_z,rtol=1e-7,atol=1e-7),'BMI z coding does not equal sample-SD scaling')
    return f,meta,dropped

@dataclass(frozen=True)
class ModelSpec:
    name:str
    precision:bool=False
    factors:int=0
    patient_mode:str='iid'
    bmi:bool=True
    loading_sd:float=.35
    precision_sd:float=.25

def model_registry(k=2,loading_sd=.35,precision_sd=.25):
    require(1<=k<=5,'Use a small, declared factor rank (1–5)')
    require(np.isfinite([loading_sd,precision_sd]).all() and min(loading_sd,precision_sd)>0,'Prior scales must be positive and finite')
    return {
      'legacy':ModelSpec('legacy',patient_mode='orthogonal'),
      'iid':ModelSpec('iid'),
      'precision':ModelSpec('precision',precision=True,precision_sd=precision_sd),
      'factor':ModelSpec('factor',factors=k,loading_sd=loading_sd),
      'combined':ModelSpec('combined',precision=True,factors=k,loading_sd=loading_sd,precision_sd=precision_sd)}

def precision_transform(frame,train_ids=None):
    """One log-read value per patient/cell type, not weighted by signature counts."""
    q=frame[['base_sample_id','global_celltype','component_total']].drop_duplicates()
    require(not q.duplicated(['base_sample_id','global_celltype']).any(),'Inconsistent read totals')
    train=q[q.component_total>=1000] if train_ids is None else q[q.base_sample_id.isin(train_ids)&(q.component_total>=1000)]
    stats={}
    for c,g in train.groupby('global_celltype'):
        z=np.log1p(g.component_total.to_numpy(float)); sd=float(z.std(ddof=1))
        require(np.isfinite(sd) and sd>1e-8,f'No read-depth variation in {c}')
        stats[c]=(float(z.mean()),sd)
    out=frame.copy();out['read_z']=[(math.log1p(r)-stats[c][0])/stats[c][1] for r,c in zip(out.component_total,out.global_celltype)]
    return out,stats

def make_arrays(frame,patient_frame=None):
    """patient_frame supplies original full design for exact legacy reconstruction."""
    f=frame.sort_values(['global_celltype','global_feature','base_sample_id']).reset_index(drop=True).copy()
    features=sorted(f.global_feature.unique());cts=sorted(f.global_celltype.unique())
    p=frame if patient_frame is None else patient_frame
    predictors=['bmi_standardized','age_z','sex_centered']
    pp=p[['base_sample_id',*predictors]].drop_duplicates()
    require(not pp.base_sample_id.duplicated().any(),'Patient covariates inconsistent')
    pp=pp.sort_values('base_sample_id');patients=pp.base_sample_id.tolist()
    for col,names,out in [('global_feature',features,'fi'),('global_celltype',cts,'ci'),('base_sample_id',patients,'pi')]:
        f[out]=f[col].map(dict(zip(names,range(len(names)))))
    require(f[['fi','ci','pi']].notna().all().all(),'Index mapping failed')
    m=f[['fi','ci']].drop_duplicates().sort_values('fi');require(len(m)==len(features),'Feature-celltype map invalid')
    x=np.column_stack([np.ones(len(pp)),pp[predictors].to_numpy(float)])
    require(np.linalg.matrix_rank(x)==4 and len(pp)>5,'Full-rank intercept/BMI/age/sex design required')
    u,s,_=np.linalg.svd(x,full_matrices=True);rank=int((s>max(x.shape)*np.finfo(float).eps*s.max()).sum())
    q=u[:,rank:]*np.sqrt(len(pp)/(len(pp)-rank))
    require(np.max(np.abs(x.T@q))<1e-9,'Invalid orthogonal basis')
    f2c=m.ci.to_numpy(int)
    # Fixed anchors: round-robin over cell types, sorted feature IDs; never selected by association.
    buckets=[list(np.where(f2c==c)[0]) for c in range(len(cts))]
    order=[b[r] for r in range(max(map(len,buckets))) for b in buckets if r<len(b)]
    return dict(frame=f,features=features,celltypes=cts,patients=patients,f2c=f2c,Q=q,X=x,order=np.array(order),
        fi=f.fi.to_numpy(int),ci=f.ci.to_numpy(int),pi=f.pi.to_numpy(int),
        y=f.score_model.to_numpy(float),bmi=f.bmi_standardized.to_numpy(float),
        cov=f[['age_z','sex_centered']].to_numpy(float),read=f.read_z.to_numpy(float))

def build_model(a,spec,nu_scope='global'):
    j=len(a['features']); n=len(a['patients']);k=spec.factors
    require(k<j,'Factor rank must be less than feature count')
    require(spec.patient_mode in ['iid','orthogonal'],'Unknown patient mode')
    require(not (k and spec.patient_mode=='orthogonal'),'Factor models use explicit iid patient/factor generative baseline')
    coords=dict(feature=a['features'],celltype=a['celltypes'],patient=a['patients'],obs_id=np.arange(len(a['y'])),
                covariate=['age_z','sex_centered'],patient_re_dimension=np.arange(a['Q'].shape[1]))
    if k:coords.update(factor=np.arange(k))
    with pm.Model(coords=coords) as model:
        alpha=pm.Normal('feature_baseline',0,.50,dims='feature')
        if spec.bmi:
            cm=pm.Normal('celltype_bmi_slope',0,.20,dims='celltype')
            bs=pm.HalfNormal('feature_sigma_bmi',.30,dims='celltype')
            br=pm.Normal('feature_bmi_raw',0,1,dims='feature')
            beta=pm.Deterministic('feature_bmi_slope',cm[a['f2c']]+bs[a['f2c']]*br,dims='feature')
        else: beta=pt.zeros(j)
        cc=pm.Normal('celltype_covariate_slope',0,.25,dims=('covariate','celltype'))
        cs=pm.HalfNormal('feature_sigma_covariate',.25,dims=('covariate','celltype'))
        cr=pm.Normal('feature_covariate_raw',0,1,dims=('covariate','feature'))
        gamma=pm.Deterministic('feature_covariate_slope',cc[:,a['f2c']]+cs[:,a['f2c']]*cr,dims=('covariate','feature'))
        ps=pm.HalfNormal('patient_global_sigma',.50)
        if spec.patient_mode=='orthogonal':
            pc=pm.Normal('patient_global_coefficient',0,ps,dims='patient_re_dimension')
            pe=pt.dot(a['Q'],pc)
        else:
            pr=pm.Normal('patient_global_raw',0,1,dims='patient');pe=ps*pr
        pe=pm.Deterministic('patient_global_effect',pe,dims='patient')
        mu=alpha[a['fi']]+beta[a['fi']]*a['bmi']+pe[a['pi']]+pt.sum(gamma[:,a['fi']].T*a['cov'],axis=1)
        loc=pm.Normal('log_sigma_location',math.log(.65),.40)
        sc=pm.HalfNormal('log_sigma_celltype_location_sd',.25)
        rc=pm.Normal('log_sigma_celltype_raw',0,1,dims='celltype')
        lc=pm.Deterministic('log_sigma_celltype_location',loc+sc*rc,dims='celltype')
        sf=pm.HalfNormal('log_sigma_feature_scale',.50,dims='celltype')
        rf=pm.Normal('log_sigma_raw',0,1,dims='feature')
        sigma=pm.Deterministic('sigma_feature',pt.exp(lc[a['f2c']]+sf[a['f2c']]*rf),dims='feature')
        sigma_obs=sigma[a['fi']]
        if spec.precision:
            # Two-sided prior: decreasing variance is tested, not imposed. Scale, not mean, changes.
            eta=pm.Normal('read_log_scale_slope',0,spec.precision_sd,dims='celltype')
            sigma_obs=sigma_obs*pt.exp(eta[a['ci']]*a['read'])
        if k:
            # Triangular anchors fix rotation/sign; no free parameters for structural zeros.
            rr,kk=np.where(np.arange(j)[:,None]>np.arange(k)[None,:])
            sd=spec.loading_sd/np.sqrt(k)
            off=pm.Normal('loading_offdiag',0,sd,shape=len(rr))
            diag=pm.HalfNormal('loading_diag',sd,shape=k)
            ordered=pt.zeros((j,k));ordered=pt.set_subtensor(ordered[rr,kk],off)
            ordered=pt.set_subtensor(ordered[np.arange(k),np.arange(k)],diag)
            load=pm.Deterministic('factor_loading',ordered[np.argsort(a['order'])],dims=('feature','factor'))
            factors=pm.Normal('factor_score',0,1,dims=('patient','factor'))
            mu=mu+pt.sum(factors[a['pi']]*load[a['fi']],axis=1)
        nu=pm.Deterministic('nu',2+pm.Exponential('nu_minus_two',.10,dims='celltype' if nu_scope=='celltype' else None),dims='celltype' if nu_scope=='celltype' else None)
        nu_obs=nu[a['ci']] if nu_scope=='celltype' else nu
        pm.StudentT('y_obs',nu=nu_obs,mu=mu,sigma=sigma_obs,observed=a['y'],dims='obs_id')
    require(np.isfinite(model.compile_logp()(model.initial_point())),'Nonfinite initial log density')
    return model

def posterior_values(idata):
    return {n:np.asarray(v.transpose('chain','draw',*[d for d in v.dims if d not in ('chain','draw')])).reshape((-1,*[v.sizes[d] for d in v.dims if d not in ('chain','draw')])) for n,v in idata.posterior.data_vars.items()}

def components(draw,a,spec,nu_scope,patient_effect=True):
    mu=draw['feature_baseline'][a['fi']]+np.sum(draw['feature_covariate_slope'][:,a['fi']].T*a['cov'],axis=1)
    if spec.bmi:mu=mu+draw['feature_bmi_slope'][a['fi']]*a['bmi']
    if patient_effect:
        mu=mu+draw['patient_global_effect'][a['pi']]
        if spec.factors:mu=mu+np.sum(draw['factor_score'][a['pi']]*draw['factor_loading'][a['fi']],axis=1)
    scale=draw['sigma_feature'][a['fi']]
    if spec.precision:scale=scale*np.exp(draw['read_log_scale_slope'][a['ci']]*a['read'])
    nu=draw['nu'][a['ci']] if nu_scope=='celltype' else draw['nu']
    return mu,scale,nu

def netcdf_save(idata,path):
    for attrs in [idata.attrs]+[idata[g].attrs for g in idata.groups()]:
        for k,v in list(attrs.items()):
            if not isinstance(v,(str,int,float,np.number,np.ndarray)) or isinstance(v,bool):attrs[k]=json.dumps(v,default=str)
    path=Path(path);tmp=path.with_name(path.name+'.tmp');idata.to_netcdf(tmp);tmp.replace(path)

def diagnostics(idata,model,profile):
    # All free RVs; selected non-constant derived estimands; structural loading zeros omitted only.
    names=[v.name for v in model.free_RVs]+[n for n in ['feature_bmi_slope','feature_covariate_slope','patient_global_effect','sigma_feature','nu'] if n in idata.posterior]
    names=list(dict.fromkeys(names));rows=[]
    rh=az.rhat(idata,var_names=names); eb=az.ess(idata,var_names=names,method='bulk');et=az.ess(idata,var_names=names,method='tail')
    for n in names:
        r=np.asarray(rh[n]);b=np.asarray(eb[n]);t=np.asarray(et[n])
        rows.append(dict(parameter=n,max_rhat=float(np.max(r)),min_ess_bulk=float(np.min(b)),min_ess_tail=float(np.min(t)),nonfinite=int((~np.isfinite(r)).sum()+(~np.isfinite(b)).sum()+(~np.isfinite(t)).sum())))
    tab=pd.DataFrame(rows);ss=idata.sample_stats
    depth=ss.maxdepth_reached if 'maxdepth_reached' in ss else (ss.tree_depth>=10)
    rec=dict(profile=profile,chains=int(idata.posterior.sizes['chain']),draws=int(idata.posterior.sizes['draw']),
        divergences=int(ss.diverging.sum()),min_bfmi=float(np.min(az.bfmi(idata))),
        max_depth_fraction=float(np.max(np.asarray(depth).mean(axis=1))),
        max_rhat=float(tab.max_rhat.max()),min_ess_bulk=float(tab.min_ess_bulk.min()),min_ess_tail=float(tab.min_ess_tail.min()),nonfinite=int(tab.nonfinite.sum()))
    rec['geometry_ok']=bool(rec['divergences']==0 and rec['min_bfmi']>=.30 and rec['max_depth_fraction']<.01 and rec['nonfinite']==0)
    rec['publication_gate']=bool(profile=='publication' and rec['chains']>=4 and rec['draws']>=4000 and rec['geometry_ok'] and rec['max_rhat']<1.01 and min(rec['min_ess_bulk'],rec['min_ess_tail'])>=400)
    rec['smoke_gate']=bool(profile=='smoke' and rec['geometry_ok'] and rec['max_rhat']<1.05 and min(rec['min_ess_bulk'],rec['min_ess_tail'])>=100)
    return rec,tab

def fit(a,spec,scope,out,profile='smoke',seed=20260918,cores=2,backend='nutpie'):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);model=build_model(a,spec,scope)
    contract=dict(version=VERSION,engine_sha=sha(__file__),model=spec.__dict__,scope=scope,profile=profile,seed=seed,
        sampler=PROFILES[profile],backend=backend,priors=PRIORS,environment=environment(),
        data_hash=hashlib.sha256(pd.util.hash_pandas_object(a['frame'],index=True).values.tobytes()).hexdigest(),
        basis_hash=hashlib.sha256(a['Q'].tobytes()).hexdigest(),feature_order=a['features'],patient_order=a['patients'])
    lock=out/'fit_contract.json';trace=out/'posterior_trace.nc'
    if lock.exists():require(json.loads(lock.read_text())==contract,f'Existing run contract differs: {out}; use a new STUDY_ID')
    else:write_json(contract,lock)
    if trace.exists():
        receipt=out/'trace_receipt.json';require(receipt.exists(),'Trace lacks receipt; inspect interrupted write')
        check_hash(trace,json.loads(receipt.read_text())['sha256']);idata=az.from_netcdf(trace)
        require(idata.attrs.get('contract_hash')==digest(contract),'Trace/contract mismatch')
    else:
        with model:
            if not (out/'prior_predictive.json').exists():
                prior=pm.sample_prior_predictive(draws=64,random_seed=seed+1,var_names=['y_obs'])
                y=prior.prior_predictive.y_obs.values
                write_json(dict(draws=64,quantiles=np.quantile(y,[.001,.01,.5,.99,.999]).tolist(),fraction_abs_gt_5=float(np.mean(np.abs(y)>5)),fraction_abs_gt_10=float(np.mean(np.abs(y)>10)),finite=bool(np.isfinite(y).all()),review_required=True),out/'prior_predictive.json')
                require(np.isfinite(y).all(),'Prior predictive simulation produced nonfinite values; review model/scales')
                del prior,y
            kw=dict(PROFILES[profile]);kw['cores']=min(cores,kw['chains'])
            if backend=='nutpie':kw.update(nuts_sampler='nutpie',nuts_sampler_kwargs={'maxdepth':10})
            elif backend=='pymc':kw.update(nuts={'max_treedepth':10})
            else:raise ValueError('Unsupported backend')
            print(f'SAMPLING {out} | {kw}',flush=True)
            start=time.time();idata=pm.sample(**kw,random_seed=seed,progressbar=True,return_inferencedata=True,idata_kwargs={'log_likelihood':False})
        idata.attrs.update(contract_hash=digest(contract),seconds=time.time()-start,study_role='post_hoc_model_extension',profile=profile)
        netcdf_save(idata,trace);write_json(dict(sha256=sha(trace)),out/'trace_receipt.json')
    diag,tab=diagnostics(idata,model,profile);write_json(diag,out/'diagnostics.json');csv(tab,out/'diagnostic_parameters.csv')
    return idata,diag

def bfdr(p,q=.05):
    p=np.asarray(p,float);require(np.isfinite(p).all() and ((p>=0)&(p<=1)).all(),'Invalid null probabilities')
    order=np.argsort(p,kind='stable');cum=np.cumsum(p[order])/np.arange(1,len(p)+1)
    n=int(np.sum(cum<=q));selected=np.zeros(len(p),bool)
    if n:selected[order[:n]]=True
    return selected

def effect_table(idata,a,out,publication=False):
    if 'feature_bmi_slope' not in idata.posterior:return None
    beta=idata.posterior.feature_bmi_slope;v=beta.values;s=v.reshape(-1,v.shape[-1]);h=az.hdi(idata,var_names=['feature_bmi_slope'],hdi_prob=.95).feature_bmi_slope.values
    tab=pd.DataFrame(dict(global_feature=a['features'],mean=s.mean(axis=0),sd=s.std(axis=0,ddof=1),hdi_low=h[:,0],hdi_high=h[:,1],prob_positive=(s>0).mean(axis=0)))
    tab['hdi_excludes_zero']=(tab.hdi_low>0)|(tab.hdi_high<0)
    for rope in [.05,.10,.15,.20]:
        ind=(np.abs(v)<=rope).astype(float);p=ind.mean(axis=(0,1));key=f'rope_{rope:.2f}'
        # Constant indicator draws do not establish zero uncertainty for rare events.
        active=np.ptp(ind,axis=(0,1))>0;mc=np.full(len(p),np.nan);ess=np.full(len(p),np.nan)
        if active.any():
            da=xr.DataArray(ind[:,:,active],dims=('chain','draw','feature'))
            ess[active]=az.ess(da,method='mean').to_array().values.ravel()
            mc[active]=az.mcse(da,method='mean').to_array().values.ravel()
        selected=bfdr(p);upper=np.where(np.isfinite(mc),np.minimum(1,p+2*mc),1)
        lower=np.where(np.isfinite(mc),np.maximum(0,p-2*mc),0)
        tab[key+'_null']=p;tab[key+'_indicator_ess']=ess;tab[key+'_mcse']=mc
        tab[key+'_selected']=selected;tab[key+'_selected_upper_mc']=bfdr(upper);tab[key+'_selected_lower_mc']=bfdr(lower)
        loo=[]
        for c in range(v.shape[0]):loo.append(bfdr(np.delete(ind,c,axis=0).mean(axis=(0,1))))
        tab[key+'_chain_stable']=np.all(np.array(loo)==selected[None,:],axis=0)
    tab['publication_sampling_gate']=publication
    tab['interpretation']='conditional estimates; predictive adequacy and robustness require separate assessment'
    csv(tab,Path(out)/'feature_effects.csv');return tab

def fixed_pairs(a,max_pairs=512,seed=713):
    # Deterministic uniform sample of within-celltype pairs; fixed before PPC values inspected.
    pairs=[(i,j) for i in range(len(a['features'])) for j in range(i+1,len(a['features'])) if a['f2c'][i]==a['f2c'][j]]
    require(len(pairs)>0,'PPC needs at least two features in a cell type')
    rng=np.random.default_rng(seed)
    take=np.sort(rng.choice(len(pairs),min(max_pairs,len(pairs)),replace=False))
    return np.asarray(pairs,dtype=int)[take]

def pair_corr(matrix,pairs):
    x=matrix[:,pairs[:,0]];y=matrix[:,pairs[:,1]];mask=np.isfinite(x)&np.isfinite(y);n=mask.sum(axis=0)
    xx=np.where(mask,x,0);yy=np.where(mask,y,0)
    xx=np.where(mask,xx-xx.sum(axis=0)/np.maximum(n,1),0)
    yy=np.where(mask,yy-yy.sum(axis=0)/np.maximum(n,1),0)
    den=np.sqrt(np.sum(xx*xx,axis=0)*np.sum(yy*yy,axis=0))
    return np.divide(np.sum(xx*yy,axis=0),den,out=np.full(len(pairs),np.nan),where=(den>0)&(n>=4))

def discrepancy(y,mu,scale,a,pairs):
    mat=np.full((len(a['patients']),len(a['features'])),np.nan);mat[a['pi'],a['fi']]=(y-mu)/scale
    centered=mat-np.nanmean(mat,axis=0);sd=np.nanstd(mat,axis=0,ddof=1)
    tail=np.nanmax(np.abs(centered),axis=0)/np.maximum(sd,1e-12)
    cor=pair_corr(mat,pairs);require(np.isfinite(cor).all(),'Undefined pairwise residual correlations')
    return sd,tail,np.array([np.mean(np.abs(cor)),np.sum(np.abs(cor)>=.90)])

def ppc(idata,a,spec,scope,out,draws=4000,seed=251):
    out=Path(out);post=posterior_values(idata);total=len(next(iter(post.values())));require(total>=draws and draws>=20 and draws%2==0,'PPC needs an even, available number of draws >=20')
    rng=np.random.default_rng(seed);ix=rng.choice(total,draws,replace=False);pairs=fixed_pairs(a)
    obs=[];rep=[];od=[];rd=[];pred=[]
    for t in ix:
        d={n:v[t] for n,v in post.items()};mu,scale,nu=components(d,a,spec,scope)
        yy=mu+scale*rng.standard_t(nu,size=len(mu))
        so,to,co=discrepancy(a['y'],mu,scale,a,pairs);sr,tr,cr=discrepancy(yy,mu,scale,a,pairs)
        obs.append(np.stack([so,to]));rep.append(np.stack([sr,tr]));od.append(co);rd.append(cr)
        # Store only a modest evenly spaced subset for in-sample predictive intervals.
        if len(pred)<200 and len(obs)%max(1,draws//200)==0:pred.append(yy)
    obs=np.asarray(obs);rep=np.asarray(rep);od=np.asarray(od);rd=np.asarray(rd);half=draws//2
    lo,hi=np.quantile(rep[:half],[.025,.975],axis=0)
    oc=((obs[half:]<lo)|(obs[half:]>hi)).sum(axis=2)
    rc=((rep[half:]<lo)|(rep[half:]>hi)).sum(axis=2)
    O=np.column_stack([oc,od[half:]]);R=np.column_stack([rc,rd[half:]])
    names=['standardized_residual_dispersion_flag_count','standardized_residual_tail_flag_count','mean_absolute_residual_correlation','residual_pairs_abs_r_ge_0p90']
    rows=[]
    for i,n in enumerate(names):
        up=(1+np.sum(R[:,i]>=O[:,i]))/(half+1);low=(1+np.sum(R[:,i]<=O[:,i]))/(half+1)
        rows.append(dict(statistic=n,observed_median=float(np.median(O[:,i])),replicated_median=float(np.median(R[:,i])),replicated_025=float(np.quantile(R[:,i],.025)),replicated_975=float(np.quantile(R[:,i],.975)),upper_p=up,lower_p=low,flag=bool(min(up,low)<.025)))
    tab=pd.DataFrame(rows);csv(tab,out/'ppc_summary.csv')
    ftab=pd.DataFrame({'global_feature':a['features']})
    for z,name in enumerate(['dispersion','tail']):
        ftab[name+'_observed_median']=np.median(obs[half:,z],axis=0)
        ftab[name+'_replicated_median']=np.median(rep[half:,z],axis=0)
        ftab[name+'_calibration_low']=lo[z];ftab[name+'_calibration_high']=hi[z]
    pi=np.quantile(np.asarray(pred),[.05,.95],axis=0);covered=(a['y']>=pi[0])&(a['y']<=pi[1])
    ftab['in_sample_90ppi_coverage']=np.bincount(a['fi'],weights=covered,minlength=len(ftab))/np.bincount(a['fi'],minlength=len(ftab))
    csv(ftab,out/'ppc_features.csv');csv(pd.DataFrame(pairs,columns=['feature_index_1','feature_index_2']),out/'ppc_fixed_pairs.csv')
    write_json(dict(version='extension-pairwise-common-patients-v1',draws=draws,calibration=half,evaluation=half,pair_count=len(pairs),pair_selection='uniform within-celltype; seed713; independent of scores',coverage='in-sample, not validation',notes='Observed and replicated residuals use matched draws. Correlations use common observed patients; no zero-fill. These PPC summaries are not identical to historical PPC v2.'),out/'ppc_contract.json')
    return tab

def inspect_legacy(base,family,a,out):
    spec=SPECS[family];path=Path(base)/'normalized_state_bhm_lowread_sensitivity_v1'/f"{spec['architecture']}__{family}__lowread1000/publication/posterior_trace.nc"
    check_hash(path,spec['trace_hash']);idata=az.from_netcdf(path)
    for dim,values in [('feature',a['features']),('celltype',a['celltypes']),('patient',a['patients'])]:
        require(idata.posterior[dim].values.astype(str).tolist()==values,f'Legacy trace coordinate mismatch: {dim}')
    require('observed_data' in idata.groups() and 'y_obs' in idata.observed_data,'Legacy trace lacks observed scores needed for endpoint verification')
    require(np.asarray(idata.observed_data.y_obs).shape==a['y'].shape and np.allclose(idata.observed_data.y_obs,a['y'],atol=1e-10,rtol=1e-10),'Legacy trace observed scores or row ordering differ')
    post=posterior_values(idata);v=post['patient_global_effect']
    require(np.max(np.abs(v@a['X']))<1e-6,'Legacy trace is not orthogonal to current patient design')
    if 'patient_global_coefficient' in post:
        require(np.allclose(v,post['patient_global_coefficient']@a['Q'].T,atol=1e-8),'Legacy basis/sign mismatch')
    model=build_model(a,model_registry()['legacy'],spec['nu']);diag,dt=diagnostics(idata,model,'publication')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    write_json(dict(source=str(path),sha256=sha(path),read_only=True,diagnostics=diag),out/'legacy_audit.json');csv(dt,out/'diagnostic_parameters.csv')
    effect_table(idata,a,out,diag['publication_gate'])
    return idata,diag

def score_from_genes(full,train_ids,manifest_path):
    """Fold-wise reconstruction. NPZs hold per-sample normalized expression BEFORE cohort gene z-scoring.
    Manifest: {components: {global_celltype: path}, signature_membership: path}.
    NPZ keys: sample_ids (unicode), genes (unicode), log1p_component_cpm (float [patient,gene]).
    Membership CSV: global_feature,gene,sign; sign is +1 or -1. A gene occurs at most once/feature.
    """
    manifest_path=Path(manifest_path);m=json.loads(manifest_path.read_text());root=manifest_path.parent
    def resolve(p): return Path(p) if Path(p).is_absolute() else root/p
    membership=pd.read_csv(resolve(m['signature_membership']),dtype={'global_feature':str,'gene':str})
    require(set(['global_feature','gene','sign'])<=set(membership),'Membership schema mismatch')
    require(not membership.duplicated(['global_feature','gene']).any(),'Duplicate gene/signature membership')
    require(membership.sign.isin([-1,1]).all(),'Signs must be -1 or +1')
    result=full.copy();result['score']=np.nan
    for cell,g in full.groupby('global_celltype'):
        require(cell in m['components'],f'No gene matrix for {cell}')
        with np.load(resolve(m['components'][cell]),allow_pickle=False) as z:
            ids=z['sample_ids'].astype(str);genes=z['genes'].astype(str);x=np.asarray(z['log1p_component_cpm'],float)
        require(x.shape==(len(ids),len(genes)) and np.isfinite(x).all() and (x>=0).all(),'Invalid log1p CPM matrix')
        require(len(set(ids))==len(ids) and len(set(genes))==len(genes),'Duplicate matrix IDs')
        require(set(g.base_sample_id)<=set(ids),'Gene matrix missing patients')
        is_train=np.isin(ids,list(train_ids));require(is_train.sum()==len(train_ids),'Training IDs mismatch')
        mean=x[is_train].mean(axis=0);sd=x[is_train].std(axis=0,ddof=1)
        # Zero-variance genes contribute zero rather than inventing a new scale.
        scaled=np.divide(x-mean,sd,out=np.zeros_like(x),where=sd>1e-12);scaled=np.clip(scaled,-5,5)
        gi={v:i for i,v in enumerate(genes)};pi={v:i for i,v in enumerate(ids)}
        for feature,fg in g.groupby('global_feature'):
            mm=membership[membership.global_feature==feature]
            require(len(mm)>0 and set(mm.gene)<=set(genes),f'Incomplete fixed gene membership: {feature}')
            value=np.zeros(len(ids))
            for sign in [1,-1]:
                names=mm.loc[mm.sign==sign,'gene']
                if len(names):value+=sign*scaled[:,[gi[v] for v in names]].mean(axis=1)
            result.loc[fg.index,'score']=[value[pi[p]] for p in fg.base_sample_id]
    require(np.isfinite(result.score).all(),'Unscored feature')
    for feature,g in result.groupby('global_feature'):
        tr=g[g.base_sample_id.isin(train_ids)].score
        sd=float(tr.std(ddof=1));require(sd>1e-10,f'Zero training signature SD: {feature}')
        result.loc[g.index,'score_model']=(g.score-tr.mean())/sd
    # Re-standardizing an existing affine z-score on training patients removes its original center/scale.
    cp=result[['base_sample_id','bmi','age_z']].drop_duplicates()
    cp=cp[cp.base_sample_id.isin(train_ids)]
    for source,target in [('bmi','bmi_standardized'),('age_z','age_z')]:
        mean=float(cp[source].mean());sd=float(cp[source].std(ddof=1));require(sd>0,'Zero covariate SD')
        result[target]=(result[source]-mean)/sd
    return result

def verify_gene_endpoint(full,manifest_path,tolerance=1e-6):
    reconstructed=score_from_genes(full,set(full.base_sample_id),manifest_path)
    err=float(np.max(np.abs(reconstructed.score.to_numpy()-full.score.to_numpy())))
    require(err<tolerance,f'Gene-to-signature endpoint differs from original score (max error {err:g}). Reconcile gene universe, membership, ddof and zero-variance handling before CV.')
    return dict(max_score_error=err,tolerance=tolerance)

def assign_folds(full,k=5,seed=129):
    require(2<=k<=10,'Use 2–10 patient folds')
    p=full[['base_sample_id','bmi_group']].drop_duplicates().sort_values('base_sample_id')
    rng=np.random.default_rng(seed);rows=[]
    for group,g in p.groupby('bmi_group'):
        ids=rng.permutation(g.base_sample_id.to_numpy())
        rows.extend(dict(base_sample_id=v,bmi_group=group,fold=i%k) for i,v in enumerate(ids))
    return pd.DataFrame(rows).sort_values('base_sample_id')

def marginal_scores(idata,test_a,spec,scope,draws=1024,nodes=64,seed=351):
    """Univariate predictive densities marginalized over NEW iid patient and factor effects.
    Mean of marginal feature log scores is a composite score, NOT a joint patient likelihood.
    Normal nuisance variance = patient_sigma^2 + sum(loadings^2). Student-t scale stays separate.
    """
    require(spec.patient_mode=='iid','Legacy projected patient basis has no implemented new-patient generative extension')
    p=posterior_values(idata);n=len(next(iter(p.values())));rng=np.random.default_rng(seed)
    ii=rng.choice(n,min(draws,n),replace=False);gh,gw=hermgauss(nodes);lw=np.log(gw)-.5*np.log(np.pi)
    logmix=np.full(len(test_a['y']),-np.inf);means=np.zeros(len(logmix))
    for i in ii:
        d={name:v[i] for name,v in p.items()};mu,scale,nu=components(d,test_a,spec,scope,patient_effect=False)
        vv=np.full(len(mu),float(d['patient_global_sigma'])**2)
        if spec.factors:vv+=np.sum(d['factor_loading'][test_a['fi']]**2,axis=1)
        for start in range(0,len(mu),256):
            sl=slice(start,min(start+256,len(mu)))
            loc=mu[sl,None]+np.sqrt(2*vv[sl,None])*gh[None,:]
            df=np.broadcast_to(nu,len(mu))[sl,None]
            logp=student_t.logpdf(test_a['y'][sl,None],df=df,loc=loc,scale=scale[sl,None])
            logmix[sl]=np.logaddexp(logmix[sl],logsumexp(logp+lw,axis=1))
        means+=mu
    logmix-=np.log(len(ii));means/=len(ii)
    out=test_a['frame'][['base_sample_id','global_feature']].copy()
    out['marginal_log_score']=logmix;out['squared_error']=(test_a['y']-means)**2
    return out

def cv_run(full,scope,registry,config,root):
    path=config.get('gene_manifest')
    require(path and Path(path).exists(), 'CV requires gene_manifest: per-component pre-z log1p CPM NPZs plus exact signed membership. Frozen full-cohort score CSVs alone cannot provide leakage-free CV. See notebook input contract.')
    require(config['profile']=='publication','CV report requires publication sampling profile')
    root=Path(root)/'cv';root.mkdir(parents=True,exist_ok=True)
    write_json(verify_gene_endpoint(full,path),root/'endpoint_reconstruction.json')
    manifest=json.loads(Path(path).read_text());r=Path(path).parent
    paths=[Path(path)]+[(Path(p) if Path(p).is_absolute() else r/p) for p in [manifest['signature_membership'],*manifest['components'].values()]]
    genes_hash={str(p):sha(p) for p in paths};lock=root/'gene_input_hashes.json'
    if lock.exists():require(json.loads(lock.read_text())==genes_hash,'CV gene inputs changed')
    else:write_json(genes_hash,lock)
    folds=assign_folds(full,config.get('cv_folds',5));csv(folds,root/'patient_folds.csv');rows=[];status=[]
    for fold in sorted(folds.fold.unique()):
        test=set(folds.loc[folds.fold==fold,'base_sample_id']);train=set(folds.base_sample_id)-test
        f=score_from_genes(full,train,path);f,stats=precision_transform(f,train)
        tr=f[(~f.low_read)&f.base_sample_id.isin(train)].copy();te=f[(~f.low_read)&f.base_sample_id.isin(test)].copy()
        a=make_arrays(tr);b=make_arrays(te)
        require(a['features']==b['features'],'CV fold lost a feature')
        require(set(a['patients']).isdisjoint(b['patients']),'Patient leakage')
        for name in config.get('cv_models',['iid','combined']):
            require(name!='legacy','Legacy model excluded from new-patient CV')
            for use_bmi in [True,False]:
                spec=replace(registry[name],bmi=use_bmi,name=name+('' if use_bmi else '_no_bmi'))
                out=root/f'fold_{fold}'/spec.name
                try:
                    idata,diag=fit(a,spec,scope,out,'publication',seed=9400+int(fold),cores=config['cores'])
                    require(diag['publication_gate'],'CV sampling diagnostics failed')
                    low=marginal_scores(idata,b,spec,scope,config.get('cv_draws',1024),32)
                    high=marginal_scores(idata,b,spec,scope,config.get('cv_draws',1024),64)
                    error=high.marginal_log_score-low.marginal_log_score
                    # Conservative numerical check at observation AND patient-average levels.
                    numerical=bool(np.max(np.abs(error))<.05)
                    require(numerical,f'Quadrature 32/64 disagreement {np.max(np.abs(error)):.3g}; increase integration accuracy before interpreting CV')
                    high['fold']=fold;high['model']=spec.name;rows.append(high)
                    csv(high,out/'heldout_marginal_scores.csv');status.append(dict(fold=int(fold),model=spec.name,passed=True))
                    del idata
                except Exception as e:
                    status.append(dict(fold=int(fold),model=spec.name,passed=False,error=str(e)))
                csv(pd.DataFrame(status),root/'cv_status.csv')
    write_json(dict(target='new patient, same predefined feature panel and eligibility rules',score='mean univariate marginal log score per patient; not joint ELPD',conditional_on='observed allocated read support; fixed original panel/QC policy',foldwise='gene means/SD, cap5 scores, feature mean/SD, BMI/age scaling, read-depth scaling',scope='no claim of external validation or validation of panel discovery',all_passed=all(s['passed'] for s in status)),root/'cv_contract.json')
    require(all(s['passed'] for s in status),'Incomplete/failed CV: status saved; no selected-fold comparison is valid')
    allrows=pd.concat(rows,ignore_index=True);csv(allrows,root/'heldout_scores.csv')
    patient=allrows.groupby(['model','base_sample_id'],as_index=False)[['marginal_log_score','squared_error']].mean()
    csv(patient,root/'patient_scores.csv');comparisons=[]
    for name in config.get('cv_models',['iid','combined']):
        pivot=patient[patient.model.isin([name,name+'_no_bmi'])].pivot(index='base_sample_id',columns='model',values='marginal_log_score').dropna()
        d=pivot[name]-pivot[name+'_no_bmi']
        comparisons.append(dict(comparison=name+' minus '+name+'_no_bmi',n_patients=len(d),mean_log_score_gain=float(d.mean()),descriptive_patient_se=float(d.std(ddof=1)/np.sqrt(len(d))),note='CV fold dependence makes this SE approximate; not a p-value; shared across correlated feature outcomes'))
    csv(pd.DataFrame(comparisons),root/'bmi_predictive_gain.csv')
    return pd.DataFrame(comparisons)

def simulate_frame(frame,scenario,seed):
    """Fixed-design stress experiments; not SBC and not independent biological evidence."""
    rng=np.random.default_rng(seed);a=make_arrays(frame);j=len(a['features']);n=len(a['patients']);c=len(a['celltypes'])
    settings={
      'null_independent':(0.,False,False,8.),
      'null_correlated':(0.,True,False,8.),
      'weak_correlated':(.05,True,False,8.),
      'moderate_heteroscedastic':(.15,False,True,8.),
      'combined':(.15,True,True,8.),
      'combined_heavy_tails':(.15,True,True,3.)}
    size,corr,hetero,df=settings[scenario]
    active=rng.choice(j,max(1,j//4),replace=False);beta=np.zeros(j);beta[active]=size*rng.choice([-1,1],len(active))
    alpha=rng.normal(0,.15,j);gamma=rng.normal(0,.08,(2,j));patient=rng.normal(0,.25,n)
    mu=alpha[a['fi']]+beta[a['fi']]*a['bmi']+np.sum(gamma[:,a['fi']].T*a['cov'],axis=1)+patient[a['pi']]
    if corr:
        # Block structure is deliberately not exactly the candidate global two-factor model.
        z=rng.normal(size=(n,c));load=rng.uniform(.35,.65,j)*rng.choice([-1,1],j)
        mu+=z[a['pi'],a['ci']]*load[a['fi']]
    sigma=rng.uniform(.4,.7,j)[a['fi']]
    if hetero:sigma*=np.exp(-.25*a['read'])
    f=a['frame'].copy();f['score_model']=mu+sigma*rng.standard_t(df,size=len(mu))
    return f,beta,dict(scenario=scenario,seed=seed,df=df,active_fraction=len(active)/j,true_abs_active_slope=size,correlated=corr,heteroscedastic=hetero)

def simulations(frame,scope,registry,config,root):
    root=Path(root)/'simulation'/config['profile'];root.mkdir(parents=True,exist_ok=True);rows=[]
    scenarios=config.get('sim_scenarios',['null_independent','null_correlated','weak_correlated','moderate_heteroscedastic','combined','combined_heavy_tails'])
    reps=config.get('sim_replicates',50);shard=config.get('sim_shard',0);shards=config.get('sim_shards',1)
    require(0<=shard<shards,'Invalid shard')
    plan=dict(scenarios=scenarios,replicates=reps,models=config.get('sim_models',['iid','precision','factor','combined']),shards=shards)
    if (root/'plan.json').exists():require(json.loads((root/'plan.json').read_text())==plan,'Simulation plan changed')
    else:write_json(plan,root/'plan.json')
    for si,scenario in enumerate(scenarios):
        for rep in range(reps):
            if (si*reps+rep)%shards!=shard:continue
            data,truth,settings=simulate_frame(frame,scenario,640000+10000*si+rep);a=make_arrays(data)
            for name in config.get('sim_models',['iid','precision','factor','combined']):
                print(f'SIMULATION {scope}: {scenario} replicate {rep+1}/{reps} model {name}',flush=True)
                out=root/scenario/f'rep_{rep:03d}'/name;row=dict(scenario=scenario,replicate=rep,model=name,passed=False)
                try:
                    cache=out/'simulation_metric_receipt.json'
                    token=digest(dict(engine=sha(__file__),environment=environment(),scope=scope,
                        model=registry[name].__dict__,sampler=PROFILES[config['profile']],settings=settings,
                        data=hashlib.sha256(pd.util.hash_pandas_object(a['frame'],index=True).values.tobytes()).hexdigest()))
                    if cache.exists():
                        saved=json.loads(cache.read_text());require(saved['token']==token,'Simulation checkpoint contract changed')
                        for rel,h in saved['files'].items():check_hash(out/rel,h)
                        rows.append(saved['row']);csv(pd.DataFrame(rows),root/f'shard_{shard}_metrics.csv');continue
                    out.mkdir(parents=True,exist_ok=True);write_json(settings,out/'truth_settings.json');csv(pd.DataFrame({'global_feature':a['features'],'true_beta':truth}),out/'truth.csv')
                    idata,diag=fit(a,registry[name],scope,out,config['profile'],seed=720000+10000*si+rep,cores=config['cores'])
                    tab=effect_table(idata,a,out,diag['publication_gate'])
                    post=idata.posterior.feature_bmi_slope.values.reshape(-1,len(truth));p=(np.abs(post)<=.10).mean(axis=0);sel=bfdr(p)
                    real=np.abs(truth)>.10
                    row.update(passed=diag['publication_gate'],bias=float(np.mean(tab['mean']-truth)),rmse=float(np.sqrt(np.mean((tab['mean']-truth)**2))),coverage=float(np.mean((truth>=tab.hdi_low)&(truth<=tab.hdi_high))),mean_interval_width=float(np.mean(tab.hdi_high-tab.hdi_low)),selected=int(sel.sum()),fdp=float(np.sum(sel&~real)/max(1,sel.sum())),power=float(np.sum(sel&real)/real.sum()) if real.any() else np.nan,geometry_ok=diag['geometry_ok'])
                    del idata
                    files={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='simulation_metric_receipt.json' and not p.name.endswith('.tmp')}
                    write_json(dict(token=token,row=row,files=files),cache)
                except Exception as e:row['error']=str(e)
                rows.append(row);csv(pd.DataFrame(rows),root/f'shard_{shard}_metrics.csv')
    return pd.DataFrame(rows)

def aggregate_simulations(root):
    simroot=Path(root)/'simulation/publication'
    paths=sorted(simroot.glob('shard_*_metrics.csv'));require(paths,'No publication-profile simulation shards')
    plan=json.loads((simroot/'plan.json').read_text())
    d=pd.concat([pd.read_csv(p) for p in paths],ignore_index=True)
    require(not d.duplicated(['scenario','replicate','model']).any(),'Duplicate simulation replicate/model rows')
    rows=[]
    for (scenario,model),g in d.groupby(['scenario','model']):
        ok=g[g.passed==True];r=dict(scenario=scenario,model=model,attempted=len(g),passed=len(ok),failed=len(g)-len(ok),failure_fraction=1-len(ok)/len(g))
        for metric in ['bias','rmse','coverage','mean_interval_width','fdp','power']:
            vals=pd.to_numeric(ok.get(metric,pd.Series(dtype=float)),errors='coerce').dropna()
            r[metric+'_mean']=float(vals.mean()) if len(vals) else None
            r[metric+'_mcse']=float(vals.std(ddof=1)/np.sqrt(len(vals))) if len(vals)>1 else None
        rows.append(r)
    expected=len(plan['scenarios'])*plan['replicates']*len(plan['models'])
    write_json(dict(expected=expected,attempted=len(d),missing=expected-len(d),complete=len(d)==expected,conditional_summary='metrics summarize passed fits only; failed/attempted counts remain explicit'),simroot/'completion.json')
    paired=[]
    for scenario,g in d.groupby('scenario'):
        good=g[g.passed==True]
        for model in [m for m in plan['models'] if m!='iid']:
            for metric in ['rmse','coverage','mean_interval_width','fdp','power']:
                if metric not in good:continue
                pivot=good.pivot(index='replicate',columns='model',values=metric)
                if not {'iid',model}<=set(pivot):continue
                pair=pivot[['iid',model]].dropna();delta=pair[model]-pair['iid']
                paired.append(dict(scenario=scenario,model=model,metric=metric,paired_passed_replicates=len(delta),planned_replicates=plan['replicates'],mean_difference=float(delta.mean()) if len(delta) else None,mcse=float(delta.std(ddof=1)/np.sqrt(len(delta))) if len(delta)>1 else None))
    csv(pd.DataFrame(paired),simroot/'paired_differences_vs_iid.csv')
    out=pd.DataFrame(rows);csv(out,simroot/'aggregate.csv');return out

def plot_reports(tables,ppcs,out):
    import matplotlib.pyplot as plt
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if 'iid' in tables:
        ref=tables['iid'].set_index('global_feature')
        sensitivity=[]
        for name,table in tables.items():
            if name in ['iid','legacy']:continue
            t=table.set_index('global_feature').loc[ref.index]
            sensitivity.append(dict(model=name,reference='iid',mean_absolute_slope_shift=float(np.mean(np.abs(t['mean']-ref['mean']))),max_absolute_slope_shift=float(np.max(np.abs(t['mean']-ref['mean']))),hdi_status_changes=int(np.sum(t.hdi_excludes_zero!=ref.hdi_excludes_zero)),selected_rope_010=int(t['rope_0.10_selected'].sum()),median_interval_width=float(np.median(t.hdi_high-t.hdi_low))))
            fig,ax=plt.subplots(figsize=(5,5));ax.scatter(ref['mean'],t['mean'],s=10,alpha=.55)
            lim=np.max(np.abs(np.r_[ref['mean'],t['mean']]))*1.1+.001;ax.plot([-lim,lim],[-lim,lim],color='grey',lw=1)
            ax.set(xlabel='IID baseline BMI slope',ylabel=f'{name}: BMI slope',title='Model sensitivity; not independent replication');fig.tight_layout()
            fig.savefig(out/f'effect_comparison_{name}.pdf');plt.close(fig)
        csv(pd.DataFrame(sensitivity),out/'model_sensitivity_summary.csv')
    if ppcs:
        names=list(ppcs);fig,axs=plt.subplots(1,2,figsize=(9,3.8))
        for ax,stat in zip(axs,['standardized_residual_tail_flag_count','mean_absolute_residual_correlation']):
            for i,name in enumerate(names):
                r=ppcs[name].set_index('statistic').loc[stat]
                ax.plot([r.replicated_025,r.replicated_975],[i,i],lw=3,color='lightgrey')
                ax.scatter(r.observed_median,i,color='black',s=22)
            ax.set_yticks(range(len(names)),names);ax.set_xlabel(stat.replace('_',' '));ax.set_title('Observed dot; replicated 95% range')
        fig.tight_layout();fig.savefig(out/'predictive_discrepancies.pdf');plt.close(fig)

def study(config):
    """Dispatch stages; no writes below the historical model directories."""
    check_environment();family=config['family'];base=Path(config['base']);scope=SPECS[family]['nu']
    root=base/'normalized_state_bhm_extensions_v1'/config['study_id']/family
    full,meta,dropped=load_frozen(base,family,config.get('pca_override'),config.get('review_override'),config.get('run_override'))
    f,depth=precision_transform(full.loc[~full.low_read].copy())
    a=make_arrays(f,full);registry=model_registry(config.get('factor_rank',2),config.get('factor_loading_sd',.35),config.get('precision_slope_sd',.25))
    root.mkdir(parents=True,exist_ok=True)
    scientific={k:v for k,v in config.items() if k not in ['stage','profile','cores','sim_shard','gene_manifest']}
    lock=dict(config=scientific,inputs=meta,engine_sha=sha(__file__),priors=PRIORS,models={k:v.__dict__ for k,v in registry.items()},post_hoc=True)
    lp=root/'study_lock.json'
    if lp.exists():require(json.loads(lp.read_text())==lock,'Study lock changed; use a new study_id and explain amendment')
    else:write_json(lock,lp)
    csv(dropped,root/'excluded_profiles.csv');csv(f[['global_feature','global_celltype','CellType','Signature']].drop_duplicates(),root/'feature_inventory.csv')
    write_json(depth,root/'read_depth_scaling.json');write_json(environment(),root/'environment.json')
    stage=config['stage'];profile=config['profile']
    print(f'{family}: {len(a["patients"])} patients, {len(a["features"])} features, {len(a["y"])} observations. Output: {root}',flush=True)
    if stage=='preflight':
        rows=[]
        for name,spec in registry.items():
            model=build_model(a,spec,scope);rows.append(dict(model=name,initial_logp=float(model.compile_logp()(model.initial_point())),free_variables=[v.name for v in model.free_RVs]))
        write_json(rows,root/'model_preflight.json');return pd.DataFrame(rows)
    if stage=='legacy':
        idata,diag=inspect_legacy(base,family,a,root/'legacy_readonly')
        ppc(idata,a,registry['legacy'],scope,root/'legacy_readonly',config.get('ppc_draws',4000));return diag
    if stage=='cv':return cv_run(full,scope,registry,config,root)
    if stage=='simulation':return simulations(f,scope,registry,config,root)
    if stage=='aggregate_simulations':return aggregate_simulations(root)
    require(stage in ['fit','report'],'Unknown stage')
    rows=[];tables={};ppcs={}
    for i,name in enumerate(config.get('models',['iid','precision','factor','combined'])):
        require(name!='legacy','Use legacy stage to read the original trace; do not overwrite/refit it here')
        spec=registry[name];out=root/name/profile
        try:
            if stage=='fit':
                idata,diag=fit(a,spec,scope,out,profile,seed=20260918+i,cores=config['cores'])
            else:
                require((out/'posterior_trace.nc').exists(),'Missing fitted trace')
                check_hash(out/'posterior_trace.nc',json.loads((out/'trace_receipt.json').read_text())['sha256'])
                idata=az.from_netcdf(out/'posterior_trace.nc');model=build_model(a,spec,scope)
                contract=json.loads((out/'fit_contract.json').read_text())
                require(contract['engine_sha']==sha(__file__) and contract['model']==spec.__dict__,'Report code/model differs from fitted contract')
                require(idata.attrs.get('contract_hash')==digest(contract),'Trace contract mismatch')
                diag,dt=diagnostics(idata,model,profile)
            tab=effect_table(idata,a,out,diag['publication_gate']);tables[name]=tab
            count=200 if profile!='publication' else config.get('ppc_draws',4000)
            if idata.posterior.sizes['chain']*idata.posterior.sizes['draw']>=count:ppcs[name]=ppc(idata,a,spec,scope,out,count)
            rows.append(dict(model=name,status='completed',**diag));del idata
        except Exception as e:
            rows.append(dict(model=name,status='failed',error=str(e)));print(f'{name}: {e}',flush=True)
            out.mkdir(parents=True,exist_ok=True);(out/'error.txt').write_text(traceback.format_exc())
        csv(pd.DataFrame(rows),root/f'{profile}_status.csv')
    plot_reports(tables,ppcs,root/f'{profile}_figures')
    write_json(dict(models_attempted=len(rows),completed=sum(r['status']=='completed' for r in rows),all_sampling_pass=all(r.get('publication_gate',False) for r in rows),automatic_winner=False,interpretation='Compare adequacy and estimation stability. Zero discoveries is an admissible result. Computational gates are not model validation.'),root/f'{profile}_review.json')
    return pd.DataFrame(rows)
