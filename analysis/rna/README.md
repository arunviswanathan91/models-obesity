# RNA pathway analyses

| Notebook | Analysis |
| --- | --- |
| `continuous_celltype_nu.ipynb` | Immune-coarse all-profile reference model, v2.2; cell-type-specific Student-t tail parameters. |
| `continuous_global_nu.ipynb` | Nonimmune and immune-fine all-profile reference models, v2.1c; shared Student-t tail parameter. |
| `categorical_bmi.ipynb` | Secondary three-category BMI model, v2.2c, including the patient-influence audit. |
| `low_read_sensitivity.ipynb` | Primary low-read-excluded continuous fits and comparison with all-profile traces. |

The two continuous notebooks implement distinct final compartment architectures for the all-profile reference fits. The low-read-excluded fits are the primary continuous analyses after diagnostic review. The global-tail model is required for the nonimmune and immune-fine analyses; it is not an earlier replacement for the immune-coarse model.

## Inputs and execution

These notebooks require the audited BayesPrism normalized-state scores, frozen PCA/QC inclusion tiers, cohort covariates and the corresponding manifests. Configure input paths in the configuration cells. Original Colab paths identify the expected directory structure. Use separate output directories for new runs; do not overwrite the archived manuscript fits.

Install the pinned sampling environment before importing PyMC. Select a compartment, sampling profile and workflow stage, then execute the relevant cells. The continuous workflows retain their preflight, sampling, convergence, posterior-predictive and posterior-summary stages. Paired PPC and patient-influence audits require existing traces and saved PPC outputs. The low-read notebook preserves its immune-fine configuration; set `FAMILIES` to select other compartments.

The scientific model code, priors, sampling settings and statistical calculations are preserved in retained cells. Standalone figure-export cells, provisional effect summaries and redundant or obsolete review cells have been omitted. Saved outputs and execution metadata are cleared. Embedded release hashes describe the curated notebook source and therefore differ from the historical notebook hashes; archived trace identities remain unchanged.

Set `NOTEBOOK_PATH_OVERRIDE` to the current curated notebook when using release verification, rather than its original Colab filename. `BASE` and `PCA_DIR_OVERRIDE` locate the audited inputs; `MODEL_ROOT` selects a separate output directory. The low-read notebook contains its own engine and requires the archived original-run paths in `FAMILY_SPECS` as well as the frozen input manifests.
