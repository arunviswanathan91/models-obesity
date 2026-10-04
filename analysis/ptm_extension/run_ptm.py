from pathlib import Path
import sys,json,gc,traceback
import numpy as np
import pandas as pd
import pymc as pm
import ptm_core as p
import run_common as c

VARIANTS=['fixed_plex','stage_subset_base','stage_adjusted','wes_subset_base','wes_adjusted',
          'mutation_adjusted','exclude_recorded_weight_loss','exclude_adenosquamous']

def load(cfg):
    work=Path(cfg['work'])/'ptm_inputs';work.mkdir(parents=True,exist_ok=True)
    source,ih=p.unpack_inputs(Path(cfg['input_zip']),work)
    if ih!=cfg['input_sha256']:raise ValueError('Multiomics input archive changed')
    links,payload,proof=p.load_signature_links(Path(cfg['signature_zip']),ih)
    data=p.load_data(source,payload);features,eligibility,arrays=p.prepare_features(data,links)
    return data,features,arrays

def sensitivity(cfg,data,features,arrays):
    root=Path(cfg['output'])/'ptm_sensitivity';root.mkdir(parents=True,exist_ok=True)
    primary=pd.read_csv(Path(cfg['report'])/'all_effects_95HDI.csv')
    if c.sha(Path(cfg['report'])/'all_effects_95HDI.csv')!=cfg['effects_sha256']:raise ValueError('Original effect table changed')
    selected=features[features.key.isin(cfg['target_keys'])].sort_values('key')
    if set(selected.key)!=set(cfg['target_keys']):raise ValueError('Some locked target features are no longer eligible')
    c.write_csv(root/'targets.csv',selected);rows=[]
    old=primary[primary.key.isin(cfg['target_keys'])&primary.variant.eq('primary')].copy()
    old['source']='historical_primary_reused';c.write_csv(root/'historical_primary_95HDI.csv',old)
    for f in selected.to_dict('records'):
        d=p.feature_cohort(data,f['assay']);a=arrays[f['key']]
        for variant in VARIANTS:
            folder=root/'fits'/f['key']/variant
            try:
                mask,adj,extra,cat=p.variant_spec(d,a['mask'],variant)
                des=p.make_design(d,mask,adj,extra,cat);Y=a['Y'][mask]
                contract=dict(feature=f,variant=variant,patient_ids=des['data'].base_sample_id.tolist(),X=c.array_hash(des['X']),Y=c.array_hash(Y),
                    terms=des['names'],group=des['group'].tolist(),sampling=cfg['sampling'],priors=cfg['ptm_model_config'],
                    input=cfg['input_sha256'],signature=cfg['signature_sha256'],source_hash=cfg['source_hash'],versions=c.versions(),seed=cfg['seed'])
                c.lock(folder,contract)
                if c.complete(folder):rec=json.loads((folder/'diagnostics.json').read_text());state='reused'
                else:
                    model=p.build_model(des,Y,cfg['ptm_model_config']);seed=int(c.digest([cfg['seed'],f['key'],variant])[:8],16)%(2**31-1)
                    z,rec=c.fit(model,folder,cfg['sampling'],seed,'pymc',['focal_effect','ptm_conditional_effect','protein_to_ptm_coupling'])
                    t=p.effect_summary(z,f,variant,des,rec['numerical_pass']);c.write_csv(folder/'effects_95HDI.csv',t)
                    pc=p.compute_ppc(z,des,Y,cfg['ptm_model_config'],folder,seed)
                    c.write_csv(folder/'design.csv',pd.DataFrame(des['X'],columns=des['names']).assign(base_sample_id=des['data'].base_sample_id))
                    c.write_json(folder/'design_notes.json',dict(n=len(Y),terms=des['names'],dropped_constant=des['dropped_constant_terms'],
                        scaling='frozen original primary paired-cohort means/SD',post_hoc=True,predictive_flagged_checks=pc))
                    c.finish(folder,['effects_95HDI.csv','predictive_checks.csv','diagnostics.json','design.csv','design_notes.json'],rec['numerical_pass'])
                    state='completed';del z,model;gc.collect()
                rows.append(dict(key=f['key'],gene=f['gene'],variant=variant,state=state,**rec))
            except Exception as ex:
                c.write_json(folder/'error.json',dict(error=str(ex),traceback=traceback.format_exc()))
                rows.append(dict(key=f['key'],gene=f['gene'],variant=variant,state='blocked_or_failed',error=str(ex)))
            c.write_csv(root/'status.csv',pd.DataFrame(rows))
    tables=[old]+[pd.read_csv(f) for f in (root/'fits').glob('*/*/effects_95HDI.csv') if c.complete(f.parent)]
    all_effects=pd.concat(tables,ignore_index=True);c.write_csv(root/'all_effects_95HDI.csv',all_effects)
    pairs=[('stage_subset_base','stage_adjusted'),('wes_subset_base','wes_adjusted')]+[('primary',v) for v in VARIANTS if v not in ['stage_adjusted','wes_adjusted']]
    out=[]
    for base,alt in pairs:
        a=all_effects[all_effects.variant.eq(base)];b=all_effects[all_effects.variant.eq(alt)]
        m=a.merge(b,on=['key','quantity','term'],suffixes=('_base','_comparison'),validate='one_to_one')
        if len(m):
            m['mean_change']=m.posterior_mean_comparison-m.posterior_mean_base
            m['interval_class_changed']=m.hdi_excludes_zero_base!=m.hdi_excludes_zero_comparison
            m['comparison']=base+' -> '+alt;out.append(m)
    if out:c.write_csv(root/'matched_comparisons.csv',pd.concat(out,ignore_index=True))

def extract_quantities(z):
    b=c.flat(z,'focal_effect')[:,0,:]
    return {'beta_protein':b[:,0],'beta_ptm':b[:,1],'beta_conditional':c.flat(z,'ptm_conditional_effect')[:,0],
            'rho':c.flat(z,'residual_correlation').ravel(),'coupling':c.flat(z,'protein_to_ptm_coupling').ravel(),
            'nu':c.flat(z,'nu').ravel()}

def generate(des,observed,cfg,scenario,seed):
    rng=np.random.default_rng(seed)
    if scenario=='sbc':
        model=p.build_model(des,observed,cfg['ptm_model_config'])
        with model:
            prior=pm.sample_prior_predictive(samples=1,random_seed=seed,
                var_names=['paired_score','focal_effect','ptm_conditional_effect','residual_correlation','protein_to_ptm_coupling','nu'])
        q={n:prior.prior[n].values.reshape((-1,*prior.prior[n].shape[2:]))[0] for n in prior.prior.data_vars}
        Y=prior.prior_predictive.paired_score.values[0,0]
        truth=dict(beta_protein=float(q['focal_effect'][0,0]),beta_ptm=float(q['focal_effect'][0,1]),
            beta_conditional=float(q['ptm_conditional_effect'][0]),rho=float(q['residual_correlation']),
            coupling=float(q['protein_to_ptm_coupling']),nu=float(q['nu']))
        return Y,truth
    sigma=np.array([.8,.9]);rho=.7 if scenario in ['protein_only','heavy_tails'] else .5
    coupling=rho*sigma[1]/sigma[0];nu=3. if scenario=='heavy_tails' else 8.
    bp=0. if scenario=='null' else .25
    bc=.2 if scenario in ['conditional_signal','heavy_tails'] else 0.
    bm=bc+coupling*bp
    coef=rng.normal(0,.03,size=(len(des['names']),2));focal=des['names'].index('bmi_5_within');coef[focal]=[bp,bm]
    # Nuisance effects vary between datasets; focal truth is fixed per scenario.
    ge=rng.normal(0,.15,size=(len(des['levels']),2));mu=des['X']@coef+ge[des['group']]
    covariance=np.outer(sigma,sigma)*np.array([[1,rho],[rho,1]])
    eps=rng.normal(size=observed.shape)@np.linalg.cholesky(covariance).T
    eps/=np.sqrt(rng.chisquare(nu,size=len(mu))/nu)[:,None]
    return mu+eps,dict(beta_protein=bp,beta_ptm=bm,beta_conditional=bc,rho=rho,coupling=coupling,nu=nu)

def calibration(cfg,data,features,arrays):
    root=Path(cfg['output'])/'ptm_calibration';root.mkdir(parents=True,exist_ok=True)
    # Coverage templates chosen by sample coverage, never by estimated BMI effect.
    ordered=features.sort_values(['n_paired','key']);templates=[ordered.iloc[0].to_dict(),ordered.iloc[-1].to_dict()]
    c.write_csv(root/'templates.csv',pd.DataFrame(templates))
    rows=[];scenarios=['null','protein_only','conditional_signal','heavy_tails','sbc']
    for scenario in scenarios:
        reps=cfg['sbc_replicates'] if scenario=='sbc' else cfg['calibration_replicates']
        for r in range(reps):
            f=templates[r%2];d=p.feature_cohort(data,f['assay']);a=arrays[f['key']]
            des=p.make_design(d,a['mask']);observed=a['Y'][a['mask']]
            folder=root/scenario/f'replicate_{r:04d}';seed=int(c.digest([cfg['seed'],scenario,r])[:8],16)%(2**31-1)
            try:
                c.lock(folder,dict(scenario=scenario,replicate=r,template=f['key'],X=c.array_hash(des['X']),group=des['group'].tolist(),
                    sampling=cfg['sampling'],priors=cfg['ptm_model_config'],source_hash=cfg['source_hash'],versions=c.versions(),seed=seed,
                    interpretation='conditional on fixed score scale and observed missingness mask; no end-to-end deconvolution or MNAR simulation'))
                if c.complete(folder):
                    rec=json.loads((folder/'diagnostics.json').read_text());state='reused'
                else:
                    yp=folder/'generated_Y.npy';tp=folder/'truth.json'
                    if yp.exists() and tp.exists():
                        truth=json.loads(tp.read_text());Y=np.load(yp,allow_pickle=False)
                        if c.array_hash(Y)!=truth.pop('_Y_hash'):raise ValueError('Generated-data checksum failure')
                    else:
                        Y,truth=generate(des,observed,cfg,scenario,seed)
                        temp=folder/'generated_Y.partial.npy';np.save(temp,Y);temp.replace(yp)
                        c.write_json(tp,{**truth,'_Y_hash':c.array_hash(Y)})
                    model=p.build_model(des,Y,cfg['ptm_model_config'])
                    z,rec=c.fit(model,folder,cfg['sampling'],seed+211,'pymc',['focal_effect','ptm_conditional_effect','protein_to_ptm_coupling'])
                    metrics=[]
                    for name,v in extract_quantities(z).items():
                        lo,hi=c.hdi(v);mu=float(v.mean());t=truth[name]
                        # Equally spaced draws within every chain, limited to 64 per chain.
                        # flat() stacks chain then draw.
                        nchain=z.posterior.sizes['chain'];ndraw=z.posterior.sizes['draw']
                        vv=v.reshape(nchain,ndraw);chosen=vv[:,np.linspace(0,ndraw-1,64,dtype=int)].ravel()
                        metrics.append(dict(scenario=scenario,replicate=r,template=f['key'],n_patients=len(Y),quantity=name,truth=t,
                            estimate=mu,bias=mu-t,squared_error=(mu-t)**2,hdi_low=lo,hdi_high=hi,covered=bool(lo<=t<=hi),
                            numerical_pass=rec['numerical_pass'],rank_less=int((chosen<t).sum()),rank_equal=int((chosen==t).sum()),rank_draws=len(chosen)))
                    c.write_csv(folder/'metrics.csv',pd.DataFrame(metrics));c.finish(folder,['metrics.csv','diagnostics.json','truth.json','generated_Y.npy'],rec['numerical_pass'])
                    state='completed';del z,model;gc.collect()
                rows.append(dict(scenario=scenario,replicate=r,state=state,**rec))
            except Exception as ex:
                c.write_json(folder/'error.json',dict(error=str(ex),traceback=traceback.format_exc()));rows.append(dict(scenario=scenario,replicate=r,state='failed',error=str(ex)))
            c.write_csv(root/'status.csv',pd.DataFrame(rows))
            aggregate_calibration(root,cfg)

def aggregate_calibration(root,cfg):
    files=[f for f in root.glob('*/replicate_*/metrics.csv') if c.complete(f.parent)]
    if not files:return
    d=pd.concat([pd.read_csv(f) for f in files],ignore_index=True);c.write_csv(root/'all_metrics.csv',d)
    rows=[]
    for scope,df in [('all_available',d),('numerically_adequate',d[d.numerical_pass])]:
        for (scenario,q),g in df.groupby(['scenario','quantity']):
            n=len(g);rate=float(g.covered.mean())
            rows.append(dict(scope=scope,scenario=scenario,quantity=q,n_datasets=n,
                planned=cfg['sbc_replicates'] if scenario=='sbc' else cfg['calibration_replicates'],
                coverage=rate,coverage_mcse=float(g.covered.astype(float).std(ddof=1)/np.sqrt(n)) if n>1 else None,
                bias=float(g.bias.mean()),bias_mcse=float(g.bias.std(ddof=1)/np.sqrt(n)) if n>1 else None,
                rmse=float(np.sqrt(g.squared_error.mean())),interpretation='SBC prior-predictive coverage' if scenario=='sbc' else 'fixed-focal-truth coverage'))
    c.write_csv(root/'calibration_summary.csv',pd.DataFrame(rows))

if __name__=='__main__':
    cfg=json.loads(Path(sys.argv[1]).read_text());data,features,arrays=load(cfg)
    if sys.argv[2]=='sensitivity':sensitivity(cfg,data,features,arrays)
    elif sys.argv[2]=='calibration':calibration(cfg,data,features,arrays)
    else:raise ValueError('Unknown stage')
