# Read-dependent residual-scale sensitivity

`read_dependent_scale.ipynb` runs the six manuscript sensitivity fits: immune coarse, immune fine and nonimmune, each with low-read-excluded and all eligible profiles. The log residual scale includes a cell-type-specific standardized log-read term with Normal(0, 0.25²) prior. The orthogonal shared patient effect, score scaling and primary model structure are retained.

`bhm_extension.py` supplies the model and diagnostics; `orthoprec.py` selects the two profile variants, reads archived primary traces, samples the six sensitivity fits and compares results. Other model options in the shared engine are not selected by this workflow.

Install the pinned notebook environment. Place both modules beside the notebook and set `WORK` to that directory. Set `BASE` to the audited data tree; provide `PCA_OVERRIDE`, `REVIEW_OVERRIDE` and `RUN_OVERRIDES` if the frozen inputs have moved. Use a separate `STUDY_ID` for new outputs. Execute the notebook cells in order. Four chains use the archived publication sampling profile; PPC uses 4,000 joint draws.

The source is the executed `BHM_PDAC_published_plus_precision_6_fits.ipynb` workflow. Its fitted model code and analysis settings are preserved. Conversational comments and notebook-specific prose were removed; the driver's engine checksum was updated to the curated module. The six-fit workflow is distinct from the later exploratory IID/factor alternatives, which are not included here.
