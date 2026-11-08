-- Regenerates web/explainer-data.json from atlas_web_rows.
-- Run each block against the models-obesity Postgres and paste the JSON into the
-- matching key. Booleans are JSON booleans in the RNA tables and strings in the
-- PTM comparison table, hence the mixed comparisons below.

-- ptm: per-feature protein, PTM, conditional slopes with kappa and rho
with e as (
  select data->>'feature' fid, data->>'component' comp, (data->>'estimate')::float est
  from atlas_web_rows where dataset_id='ptm/effects'
    and data->>'contrast' in ('bmi_5_within','residual')
), p as (
  select fid,
    max(est) filter (where comp='Protein') bp,
    max(est) filter (where comp='PTM') bm,
    max(est) filter (where comp='PTM_given_Protein') cond,
    max(est) filter (where comp='protein_to_ptm_coupling') kappa,
    max(est) filter (where comp='residual_Protein_PTM_correlation') rho
  from e group by fid
), c as (
  select data->>'feature_id' fid, data->>'assay' assay, data->>'interval_transition' tr,
         data->>'candidate' cand, data->>'gene' gene,
         (data->>'n_patients')::int np, (data->>'n_plexes')::int npx,
         row_number() over (partition by data->>'interval_transition'
           order by abs((data->>'posterior_mean_conditional')::float) desc) rk
  from atlas_web_rows where dataset_id='ptm/marginal_conditional_feature_comparison'
    and data->>'numerical_checks_pass'='True'
)
select jsonb_agg(jsonb_build_array(c.gene,c.assay,c.tr,c.cand='True',
  round(p.bm::numeric,4),round(p.bp::numeric,4),round(p.cond::numeric,4),
  round(p.kappa::numeric,4),round(p.rho::numeric,4),c.np,c.npx)
  order by c.tr, abs(p.cond) desc)
from p join c using(fid) where c.cand='True' or c.rk<=8;

-- funnel and transitions
select jsonb_agg(data order by row_number) from atlas_web_rows where dataset_id='ptm/ptm_coverage_counts';
select jsonb_agg(data order by row_number) from atlas_web_rows where dataset_id='ptm/interval_transition_counts';
select jsonb_agg(data order by row_number) from atlas_web_rows where dataset_id='ptm/scope_summary';

-- widening: fmax is the largest widening factor an interval survives
with s as (
  select data->>'cell_type' ct,
         replace(split_part(data->>'feature','||',3),'_Signature','') sig,
         (data->>'estimate')::float est,
         (data->>'hdi_lower')::float lo,(data->>'hdi_upper')::float hi
  from atlas_web_rows
  where dataset_id in ('immune_coarse/effects','immune_fine/effects','nonimmune/effects')
    and data->>'entity_level'='feature'
), f as (
  select ct,sig,est,lo,hi,
    case when est>0 then est/nullif(est-lo,0) else (-est)/nullif(hi-est,0) end fmax
  from s where (lo>0 or hi<0)
)
select jsonb_build_object(
  'total',(select count(*) from s),'excludes_zero',(select count(*) from f),
  'survive_137',(select count(*) from f where fmax>1.3745),
  'rows',(select jsonb_agg(jsonb_build_array(ct,sig,round(est::numeric,4),
            round(lo::numeric,4),round(hi::numeric,4),round(fmax::numeric,3))
          order by fmax desc) from f where fmax>1.20));

-- support metrics per deconvolution component
select jsonb_agg(jsonb_build_array(data->>'CellType',data->>'compartment',data->>'analysis_tier',
  round(nullif(data->>'median_component_total','')::numeric,0),
  round(nullif(data->>'minimum_component_total','')::numeric,0),
  round(nullif(data->>'median_theta','')::numeric,5),
  round(nullif(data->>'median_theta_cv','')::numeric,3),
  round(nullif(data->>'fraction_theta_below_0_001','')::numeric,3),
  (data->>'strict_qc_keep')='true',(data->>'tiered_primary_keep')='true')
  order by data->>'compartment', nullif(data->>'median_component_total','')::numeric)
from atlas_web_rows where dataset_id='reference/component_inclusion_tiers';

-- simulation operating characteristics, Bayes against OLS
with t as (
  select data->>'scenario' sc, data->>'kind' kd, data->>'contrast' ct,
    case data->>'rule' when 'Bayes_estimation' then 'B' else 'O' end rule,
    avg(nullif(data->>'rmse','')::float) rmse, avg(nullif(data->>'coverage','')::float) cov,
    avg(nullif(data->>'mean_interval_width','')::float) wid,
    avg(nullif(data->>'bias','')::float) bias,
    avg(nullif(data->>'practical_power','')::float) pw,
    avg(nullif(data->>'practical_fdp','')::float) fdp
  from atlas_web_rows where dataset_id='simulation/simulation_task_estimation_metrics'
    and data->>'geometry'='generic_K8' and data->>'contrast' in ('bmi','obese_vs_normal')
  group by 1,2,3,4
)
select jsonb_agg(jsonb_build_array(sc,kd,ct,rule,round(rmse::numeric,4),round(cov::numeric,4),
  round(wid::numeric,4),round(bias::numeric,4),round(pw::numeric,4),round(fdp::numeric,4))
  order by sc,ct,rule) from t;

-- rna-protein: within-plex slopes, their difference and the residual correlation
with r as (
  select data->>'feature' f,max(data->>'label') lab,
    max((data->>'estimate')::float) filter (where data->>'contrast'='RNA_BMI_per5') rna,
    max((data->>'estimate')::float) filter (where data->>'contrast'='Protein_BMI_per5') pro,
    max((data->>'estimate')::float) filter (where data->>'contrast'='Protein_minus_RNA_BMI_per5') d,
    max((data->>'estimate')::float) filter (where data->>'contrast'='residual_RNA_Protein_correlation') rho
  from atlas_web_rows where dataset_id='rna_protein/effects' group by 1
)
select jsonb_agg(jsonb_build_array(lab,round(rna::numeric,4),round(pro::numeric,4),
  round(d::numeric,4),round(rho::numeric,4)) order by abs(d) desc)
from r where rna is not null and pro is not null and rho is not null;
