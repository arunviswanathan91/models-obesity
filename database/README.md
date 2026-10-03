# PDAC obesity atlas database

Project: models-obesity (lab-owned Supabase).

## Data layers
- atlas_runs: immutable source report manifests, trace hashes, model and compartment; publication status defaults to staged.
- atlas_sources: source file identity, SHA-256, expected row count and import status.
- atlas_records: original parsed rows. CSV values remain strings to preserve source representation.
- atlas_effects: typed plotting estimates with explicit units and interval bounds, linked to original rows.
- atlas_diagnostics: archived checks linked to the same source trace.
- atlas_figure_specs: renderer version and display parameters, separately tagged as paper_main, paper_supplement or extended. Disabled by default.
- atlas_import_batches: import completion and verification receipts.

## Access
All seven tables have RLS enabled. anon and authenticated have no privileges. Backend service access only; never place a privileged key in frontend code.
No data is published until its run and rendering recipes are reviewed.

## Scientific rules
Continuous report estimates come from lowread1000 runs. Earlier unfiltered continuous fits are different runs and must never share the same diagnostic linkage.
Categorical estimates use the archived publication/publication_long profiles.
Keep the fine-immune categorical provisional flag and all predictive-check caveats.
Continuous BMI units are score SD per 5 kg/m2; categorical contrasts are score SD between groups.
95% interval exclusion is descriptive and is not a multiplicity-controlled discovery declaration.
No new fitting. Exclude unfinished remaining-analysis shards.
Import full result tables even when some rows appear in the paper, but restrict the web figure catalogue to additional displays.
Keep patient-level source files and full posterior traces out of public tables.

## Imports
The first batch targets the six RNA report families and matched diagnostics. Extension, simulation, full parameter/chain diagnostics and browser rendering recipes are subsequent imports, not automatically inferred from the first batch.
SQL files here define the schema and typed RNA normalization. Source data and Drive identifiers are retained privately in Supabase, not committed to this public repository.
