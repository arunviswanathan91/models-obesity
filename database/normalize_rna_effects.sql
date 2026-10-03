insert into public.atlas_effects(source_key,row_number,run_key,entity_level,feature_key,cell_type,signature,contrast,estimate,hdi_lower,hdi_upper,units,interval_excludes_zero)
select r.source_key,r.row_number,s.run_key,
case when s.source_path like '%celltype_%' then 'celltype' else 'feature' end,
coalesce(r.payload->>'global_feature',r.payload->>'global_celltype'),
r.payload->>'celltype_label',r.payload->>'signature_label',
case when u.exposure_model='continuous' then 'bmi_per_5kg_m2' else r.payload->>'contrast' end,
(case when u.exposure_model='continuous' then r.payload->>'posterior_mean_per_5kg' else r.payload->>'posterior_mean' end)::double precision,
(case when u.exposure_model='continuous' then r.payload->>'hdi_95_low_per_5kg' else r.payload->>'hdi_95_low' end)::double precision,
(case when u.exposure_model='continuous' then r.payload->>'hdi_95_high_per_5kg' else r.payload->>'hdi_95_high' end)::double precision,
case when u.exposure_model='continuous' then 'score SD per 5 kg/m2 BMI' else 'score SD between BMI groups' end,
(r.payload->>'hdi_95_excludes_zero')::boolean
from public.atlas_records r join public.atlas_sources s using(source_key) join public.atlas_runs u using(run_key)
where s.import_status='verified' and (s.source_path like '%feature_bmi_posterior.csv' or s.source_path like '%celltype_bmi_posterior.csv' or s.source_path like '%feature_contrasts.csv' or s.source_path like '%celltype_contrasts.csv')
on conflict do nothing;
insert into public.atlas_diagnostics(source_key,row_number,run_key,diagnostic_kind,payload)
select r.source_key,r.row_number,s.run_key,'report_ppc_domain_checks',r.payload
from public.atlas_records r join public.atlas_sources s using(source_key)
where s.import_status='verified' and s.source_path like 'BiB_figure_report_v4/%ppc_domain_checks.csv'
on conflict do nothing;
update public.atlas_runs set source_manifest=source_manifest || '{"continuous_selection":"lowread1000","diagnostics_linkage":"exact_trace_sha256"}'::jsonb where exposure_model='continuous';
