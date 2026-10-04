"""New-output-only, per-chain checkpointing. Never writes to historical runs."""
from pathlib import Path
import hashlib,json,time,os,importlib.metadata,gc
import numpy as np
import pandas as pd
import pymc as pm
import arviz as az

def digest(x):
    return hashlib.sha256(json.dumps(x,sort_keys=True,default=str).encode()).hexdigest()
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def array_hash(x):
    x=np.ascontiguousarray(x)
    return digest([str(x.dtype),list(x.shape),hashlib.sha256(x.tobytes()).hexdigest()])
def write_json(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_name(p.name+'.partial');q.write_text(json.dumps(x,indent=2,default=lambda z:z.item() if isinstance(z,np.generic) else str(z)));os.replace(q,p)
def write_csv(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_name(p.name+'.partial');d.to_csv(q,index=False);os.replace(q,p)
def versions():
    out={}
    for n in ['pymc','arviz','numpy','pytensor','pandas','scipy','xarray','nutpie','h5netcdf']:
        try:out[n]=importlib.metadata.version(n)
        except importlib.metadata.PackageNotFoundError:pass
    return out
def lock(folder,contract):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    p=folder/'contract.json'
    contract=json.loads(json.dumps(contract,default=str))
    if p.exists():
        if json.loads(p.read_text())!=contract:
            raise ValueError(f'Configuration/data/code changed at {folder}; use a new RUN_LABEL. Do not edit old contracts.')
    else:write_json(p,contract)
def complete(folder):
    p=Path(folder)/'complete.json'
    if not p.exists():return False
    m=json.loads(p.read_text())
    for name,h in m['files'].items():
        f=Path(folder)/name
        if not f.is_file() or sha(f)!=h:raise ValueError(f'Completed output missing/changed: {f}')
    return True
def finish(folder,files,passed):
    write_json(Path(folder)/'complete.json',dict(numerical_pass=bool(passed),files={n:sha(Path(folder)/n) for n in files}))
def hdi(v):
    v=np.sort(np.asarray(v).ravel());k=int(np.floor(.95*len(v)))
    if len(v)<20 or not np.isfinite(v).all():raise ValueError('Incomplete or nonfinite posterior')
    i=int(np.argmin(v[k:]-v[:-k]));return float(v[i]),float(v[i+k])
def flat(idata,name):
    a=idata.posterior[name];dims=[x for x in a.dims if x not in ('chain','draw')]
    return a.stack(sample=('chain','draw')).transpose('sample',*dims).values
def load_trace(p,chains,draws):
    p=Path(p);rec=p.with_suffix(p.suffix+'.json')
    if not rec.exists() or sha(p)!=json.loads(rec.read_text())['sha256']:raise ValueError(f'Trace checksum failure: {p}')
    z=az.from_netcdf(p);z.load();z.close()
    if z.posterior.sizes.get('chain')!=chains or z.posterior.sizes.get('draw')!=draws:
        raise ValueError(f'Incomplete trace, including zero-draw traces, is NOT reusable: {p}')
    if not all(np.isfinite(v.values).all() for v in z.posterior.data_vars.values()):raise ValueError(f'Nonfinite trace: {p}')
    return z
def save_trace(p,z):
    # nutpie 0.16 stores nested metadata on InferenceData; netCDF attrs must be scalar/array.
    for attrs in [z.attrs]+[z[g].attrs for g in z.groups()]:
        for key,value in list(attrs.items()):
            if isinstance(value,(bool,np.bool_)) or not isinstance(value,(str,bytes,int,float,np.number,np.ndarray)):
                attrs[key]=json.dumps(value,default=str)
    p=Path(p);tmp=p.with_name(p.name+'.partial');z.to_netcdf(tmp);os.replace(tmp,p)
    write_json(p.with_suffix(p.suffix+'.json'),{'sha256':sha(p)})
def sample_chains(model,folder,cfg,seed,attempt,backend):
    folder=Path(folder)/f'attempt_{attempt+1}';folder.mkdir(parents=True,exist_ok=True)
    draws=cfg['draws']*2**attempt;tune=cfg['tune']*2**attempt
    lock(folder,dict(seed=seed,draws=draws,tune=tune,chains=cfg['chains'],backend=backend,target_accept=cfg['target_accept'],maxdepth=cfg['maxdepth']))
    pieces=[]
    for chain in range(cfg['chains']):
        dest=folder/f'chain_{chain:02d}.nc'
        if dest.exists():
            print('REUSE chain',dest,flush=True);part=load_trace(dest,1,draws)
        else:
            print(f'SAMPLE {folder.parent.name} attempt={attempt+1} chain={chain+1}/{cfg["chains"]} tune={tune} draws={draws}',flush=True)
            with model:
                extra=({'nuts_sampler_kwargs':{'maxdepth':cfg['maxdepth']}} if backend=='nutpie'
                       else {'nuts':{'max_treedepth':cfg['maxdepth']}})
                part=pm.sample(chains=1,cores=1,tune=tune,draws=draws,random_seed=(seed+chain*100003+attempt*1009)%(2**31-1),
                    target_accept=cfg['target_accept'],nuts_sampler=backend,progressbar=True,
                    return_inferencedata=True,compute_convergence_checks=False,idata_kwargs={'log_likelihood':False},**extra)
            if part.posterior.sizes.get('draw')!=draws or part.posterior.sizes.get('chain')!=1:
                raise RuntimeError('Sampling interrupted: expected full posterior draws. Partial trace was not accepted.')
            if not all(np.isfinite(v.values).all() for v in part.posterior.data_vars.values()):raise RuntimeError('Nonfinite samples')
            save_trace(dest,part)
        pieces.append(part)
    return az.concat(*pieces,dim='chain',reset_dim=True) if len(pieces)>1 else pieces[0]
def diagnostics(z,model,cfg,derived=()):
    names=list(dict.fromkeys([v.name for v in model.free_RVs]+list(derived)))
    tab=az.summary(z,var_names=names,kind='diagnostics',round_to='none')
    vals=tab[['r_hat','ess_bulk','ess_tail']].to_numpy();finite=bool(np.isfinite(vals).all())
    ss=z.sample_stats
    if 'maxdepth_reached' in ss:depth=int(ss.maxdepth_reached.sum())
    elif 'reached_max_treedepth' in ss:depth=int(ss.reached_max_treedepth.sum())
    elif 'tree_depth' in ss:depth=int((ss.tree_depth>=cfg['maxdepth']).sum())
    elif 'depth' in ss:depth=int((ss.depth>=cfg['maxdepth']).sum())
    else:raise ValueError('No sampler tree-depth diagnostic')
    rec=dict(max_rhat=float(tab.r_hat.max()),min_ess_bulk=float(tab.ess_bulk.min()),min_ess_tail=float(tab.ess_tail.min()),
             divergences=int(ss.diverging.sum()),min_bfmi=float(np.min(az.bfmi(z))),depth_hits=depth,finite=finite)
    rec['numerical_pass']=bool(finite and rec['max_rhat']<=1.01 and min(rec['min_ess_bulk'],rec['min_ess_tail'])>=400
                               and rec['divergences']==0 and rec['min_bfmi']>=.3 and depth==0)
    return rec,tab.rename_axis('parameter').reset_index()
def fit(model,folder,cfg,seed,backend,derived=()):
    for attempt in range(2 if cfg['retry_once'] else 1):
        z=sample_chains(model,folder,cfg,seed,attempt,backend)
        rec,tab=diagnostics(z,model,cfg,derived);rec['attempt']=attempt+1
        write_csv(Path(folder)/f'diagnostics_attempt_{attempt+1}.csv',tab)
        write_json(Path(folder)/f'diagnostics_attempt_{attempt+1}.json',rec)
        if rec['numerical_pass'] or attempt==int(cfg['retry_once']):break
        del z;gc.collect()
    write_json(Path(folder)/'diagnostics.json',rec)
    return z,rec
