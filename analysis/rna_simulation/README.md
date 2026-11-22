# Cellular RNA simulation

[cellular_rna_simulation.ipynb](cellular_rna_simulation.ipynb) generates synthetic scores, fits continuous and categorical BMI models, and saves posterior traces, diagnostics and performance summaries for Figures 2 and S1.

The manuscript used 20 conditions × 10 datasets × 2 models = 400 assessments. The reference design has 127 synthetic patients, with BMI groups of 52, 57 and 18, and 48 features across eight cell types. Other conditions vary effect size and shape, sample size, feature layout, residual dependence, heavy tails, contamination, sparsity and patient effects.

## Run

1. Open the notebook in Colab and run the dependency cell.
2. Set `OUTPUT_BASE` to a new output directory. Use four runtimes with `SHARD_ID=0`, `1`, `2` and `3`.
3. Keep `REPLICATES_TO_RUN=10`, set `RUN_SIMULATION=True` and run the cells in order. Each shard has 100 tasks and resumes completed work.
4. After all shards finish, use one runtime with `RUN_SIMULATION=False` and `AGGREGATE=True`. Set `DOWNLOAD_REPORTS=True` to download the summaries.

Keep `PLANNED_REPLICATES=50`: it belongs to the original experiment contract, although only ten repetitions were run. Sampling uses four chains, 2,000 warmup iterations and 4,000 draws, with an 8,000-draw escalation when required. No patient-level input files are needed for this synthetic study.

Outputs include `experiment.json`, per-task synthetic data, posterior traces, diagnostics, posterior means and 95% HDIs, and reports under `reports/shard_<ID>/`. The manuscript's pointwise HDI summaries use the 399 numerically adequate assessments; the notebook also retains its separate ROPE/BFDR precision checks. These decision rules are distinct.

The fitting code matches the saved v5 shard notebooks for `focused_0f5aa4b532be5d0d`. Source hashes and checks are recorded in [provenance.json](provenance.json). Package versions contribute to the experiment ID; a different environment produces a different ID.
