# Cellular RNA: Student-t patient × cell-type intercept model

Primary cellular RNA analysis in the paper (Equations 2–5). Inputs, the ≥1,000-read rule, priors, Student-t likelihood and convergence gates are the same as in `analysis/rna`; the shared patient effect is replaced by a Student-t patient × cell-type intercept.

| File | Use in the paper |
|---|---|
| `cellular_rna_ri_hc_refit.ipynb`, `cellular_rna_ri_hc_refit_sharded.ipynb` | Primary fits for coarse immune, fine immune and non-immune (single runtime or sharded) |
| `patch_student_t_intercept.py` | Makes the intercept Student-t; applied after the engine cell of the refit notebook |
| `primary_prior_sensitivity.ipynb` | Prior sensitivity |
| `primary_readscale_sensitivity.ipynb` | Read-dependent residual scale (Figures S7–S8) |
| `primary_categorical_allprofiles.ipynb` | WHO BMI-group model and all-profile fits (Figure 4, S9) |
| `immune_coarse_allprofiles_long.ipynb` | Longer-chain refit of the coarse immune all-profile fit |
| `primary_ppc_diagnostics_screen.ipynb` | Posterior predictive checks, patient diagnostics and screen summaries |
| `extras_student_t.ipynb` | Random-gene control and TAM gene drivers (Figure S11) |
| `within_celltype_contrasts.py` | Within-cell-type contrasts from saved traces |

The orthogonal patient-effect model in `analysis/rna` is kept as the comparison model.
