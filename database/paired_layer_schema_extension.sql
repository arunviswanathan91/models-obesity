alter table public.atlas_effects add column effect_component text not null default 'estimate', add column variant text, add column assay text, add column numerical_checks_pass boolean;
alter table public.atlas_effects drop constraint atlas_effects_pkey;
alter table public.atlas_effects add primary key(source_key,row_number,effect_component);
create index atlas_effects_layer_idx on public.atlas_effects(run_key,assay,variant,effect_component);