# PTM sensitivity and calibration

This directory contains the final analysis code for 24 post hoc PTM candidates across eight sensitivity specifications (192 fits), and 250 simulated datasets across five calibration scenarios. It uses the archived source from run `ee305d82dbea549b`. The primary paired PTM estimates are read from the frozen primary analysis; they are not refitted here.

## Code

- [run.py](run.py): configurable entry point for sensitivity, calibration, or both.
- [run_ptm.py](run_ptm.py): candidate selection, matched sensitivity comparisons, simulation generation and calibration summaries.
- [ptm_core.py](ptm_core.py): frozen paired protein/PTM preparation, design matrices, model and posterior summaries.
- [run_common.py](run_common.py): per-chain checkpoints, sampling, numerical diagnostics and integrity checks.
- [settings.json](settings.json): fixed candidate keys, priors, seeds, sampling settings and input checksums.
- [requirements.txt](requirements.txt): package versions recorded in the completed fits.
- [provenance.json](provenance.json): source fingerprints and curation validation.

## Run

Install dependencies in an isolated Python environment. In Colab, clone the repository and run these commands from its root after mounting the storage containing the frozen inputs.

```bash
python -m pip install -r analysis/ptm_extension/requirements.txt
python analysis/ptm_extension/run.py \
  --stage both \
  --input-zip /path/to/PDAC_main_cohort_multiomics_inputs.zip \
  --signature-zip /path/to/PDAC_original_signature_inputs.zip \
  --primary-report /path/to/primary/aggregate \
  --output /path/to/new/ptm_extension_results \
  --work /path/to/local/ptm_extension_work
```

`--primary-report` must contain the original `all_effects_95HDI.csv`. Input archives, patient-level data and posterior traces are not distributed in this repository. Their frozen hashes are checked before fitting. Package versions must match the archived analysis. The two stages may also be run separately with `--stage sensitivity` or `--stage calibration` using the same paths.

The entry point executes each stage sequentially with the archived deterministic seeds. It uses a fingerprint of the curated source and a separate output directory; published historical checkpoint contracts are not overwritten. The original publication run used six scheduling shards; sharding did not alter model definitions, priors, data masks, seeds or per-fit sampling settings.

The eight specifications are fixed plex effects; stage subset; stage adjustment within that subset; purity subset; purity adjustment within that subset; mutation adjustment; exclusion of recorded weight loss; and exclusion of adenosquamous histology. Stage and purity adjustment comparisons use their matched subset baselines. Other comparisons use the primary fit.

Calibration uses 50 datasets each for the null, protein-only, conditional-signal, heavy-tail and prior-predictive SBC scenarios. Templates retain the observed paired-cohort covariates and missingness masks. Coverage and interval support are pointwise; they do not establish screen-wide false-discovery control.

Scientific calculations, feature scope and thresholds are retained. One user-directed missing-file message and one instructional comment were neutralized. AST comparison found no computational changes. These checks validate source curation; they do not constitute a new sampling run.
