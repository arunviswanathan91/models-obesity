
create table public.atlas_runs (
 run_key text primary key,
 modality text not null,
 compartment text,
 exposure_model text not null,
 source_manifest jsonb not null default '{}'::jsonb,
 trace_sha256 text,
 provisional boolean not null default false,
 publication_status text not null default 'staged' check(publication_status in ('staged','published','excluded')),
 created_at timestamptz not null default now()
);
create table public.atlas_sources (
 source_key text primary key,
 run_key text not null references public.atlas_runs(run_key),
 drive_file_id text,
 source_path text not null,
 sha256 text not null check(length(sha256)=64),
 expected_rows integer not null check(expected_rows>=0),
 columns_json jsonb not null,
 import_status text not null default 'pending' check(import_status in ('pending','verified','failed')),
 imported_at timestamptz not null default now(),
 unique(run_key,source_path,sha256)
);
create index atlas_sources_run_idx on public.atlas_sources(run_key);
create table public.atlas_records (
 source_key text not null references public.atlas_sources(source_key),
 row_number integer not null check(row_number>0),
 payload jsonb not null check(jsonb_typeof(payload)='object'),
 primary key(source_key,row_number)
);
create table public.atlas_effects (
 source_key text not null,
 row_number integer not null,
 run_key text not null references public.atlas_runs(run_key),
 entity_level text not null check(entity_level in ('feature','celltype','programme')),
 feature_key text not null,
 cell_type text,
 signature text,
 gene text,
 contrast text not null,
 estimate double precision not null,
 hdi_lower double precision not null,
 hdi_upper double precision not null,
 units text not null,
 interval_probability double precision not null default 0.95 check(interval_probability>0 and interval_probability<1),
 interval_excludes_zero boolean not null,
 primary key(source_key,row_number),
 foreign key(source_key,row_number) references public.atlas_records(source_key,row_number),
 check(hdi_lower<=hdi_upper)
);
create index atlas_effects_filter_idx on public.atlas_effects(run_key,entity_level,contrast,cell_type);
create index atlas_effects_feature_idx on public.atlas_effects(feature_key);
create table public.atlas_diagnostics (
 source_key text not null,
 row_number integer not null,
 run_key text not null references public.atlas_runs(run_key),
 diagnostic_kind text not null,
 payload jsonb not null,
 primary key(source_key,row_number),
 foreign key(source_key,row_number) references public.atlas_records(source_key,row_number)
);
create index atlas_diagnostics_run_idx on public.atlas_diagnostics(run_key,diagnostic_kind);
create table public.atlas_figure_specs (
 figure_key text primary key,
 run_key text references public.atlas_runs(run_key),
 title text not null,
 renderer_key text not null,
 renderer_version text not null,
 parameters jsonb not null default '{}'::jsonb,
 placement text not null default 'extended' check(placement in ('paper_main','paper_supplement','extended')),
 enabled boolean not null default false
);
create index atlas_figure_specs_run_idx on public.atlas_figure_specs(run_key);
create table public.atlas_import_batches (
 batch_key text primary key,
 started_at timestamptz not null default now(),
 completed_at timestamptz,
 status text not null check(status in ('running','verified','failed')),
 report jsonb not null default '{}'::jsonb
);
alter table public.atlas_runs enable row level security;
revoke all on public.atlas_runs from anon, authenticated;
grant all on public.atlas_runs to service_role;
create policy service_access on public.atlas_runs for all to service_role using(true) with check(true);
alter table public.atlas_sources enable row level security;
revoke all on public.atlas_sources from anon, authenticated;
grant all on public.atlas_sources to service_role;
create policy service_access on public.atlas_sources for all to service_role using(true) with check(true);
alter table public.atlas_records enable row level security;
revoke all on public.atlas_records from anon, authenticated;
grant all on public.atlas_records to service_role;
create policy service_access on public.atlas_records for all to service_role using(true) with check(true);
alter table public.atlas_effects enable row level security;
revoke all on public.atlas_effects from anon, authenticated;
grant all on public.atlas_effects to service_role;
create policy service_access on public.atlas_effects for all to service_role using(true) with check(true);
alter table public.atlas_diagnostics enable row level security;
revoke all on public.atlas_diagnostics from anon, authenticated;
grant all on public.atlas_diagnostics to service_role;
create policy service_access on public.atlas_diagnostics for all to service_role using(true) with check(true);
alter table public.atlas_figure_specs enable row level security;
revoke all on public.atlas_figure_specs from anon, authenticated;
grant all on public.atlas_figure_specs to service_role;
create policy service_access on public.atlas_figure_specs for all to service_role using(true) with check(true);
alter table public.atlas_import_batches enable row level security;
revoke all on public.atlas_import_batches from anon, authenticated;
grant all on public.atlas_import_batches to service_role;
create policy service_access on public.atlas_import_batches for all to service_role using(true) with check(true);