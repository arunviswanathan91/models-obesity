insert into atlas_effects(source_key,row_number,run_key,entity_level,feature_key,cell_type,signature,gene,contrast,estimate,hdi_lower,hdi_upper,units,interval_excludes_zero,effect_component,variant,assay,numerical_checks_pass)
select r.source_key,r.row_number,s.run_key,
case when s.run_key='rna_protein' then 'programme' else 'feature' end,
case when s.run_key='rna_protein' then r.payload->>'programme' else r.payload->>'feature_id' end,
null,null,r.payload->>'gene',
case when s.run_key='rna_protein' then r.payload->>'quantity' else r.payload->>'term' end,
(r.payload->>'posterior_mean')::double precision,
(r.payload->>'hdi_95_lower')::double precision,
(r.payload->>'hdi_95_upper')::double precision,
r.payload->>'effect_units',
((r.payload->>'hdi_95_lower')::double precision>0 or (r.payload->>'hdi_95_upper')::double precision<0),
r.payload->>'quantity',r.payload->>'variant',r.payload->>'assay',
(r.payload->>'numerical_checks_pass')::boolean
from atlas_records r join atlas_sources s using(source_key)
where s.import_status='verified' and (s.source_path like '%unique_bulk_programme_effects_95HDI.csv' or s.source_path like '%sensitivity_effects_95HDI.csv' or s.source_path='aggregate/all_effects_95HDI.csv')
on conflict do nothing;
