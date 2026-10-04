# Compartment QC and PCA

[compartment_qc_and_pca.ipynb](compartment_qc_and_pca.ipynb) contains the final compartment audit, endpoint comparison and weighted publication PCA analysis. It reads completed BayesPrism deconvolution objects and expression components; it does not perform deconvolution or fit the BMI model.

The analysis validates theta and component-expression identities, reconstructs the 127-patient cohort, checks signature mapping, missingness and zero counts, compares score definitions, evaluates empirical dependence, and constructs tiered and strict-QC weighted PCA summaries.

Run the compartment audit for `immune_course`, `immune_fine` and `non_immune` separately using `TARGET_COMPARTMENT`. The combined endpoint and PCA cells require saved audit outputs for all three compartments. Configure `BASE`, `BAYESPRISM_ROOT` and the RDS override in the initial configuration cell to match the mounted input directory. Input RDS objects, patient matrices, clinical data and the frozen signature catalog are required and are not bundled here.

The publication PCA source is the final appended PCA analysis in the September 6 notebook. Its output report identifies `within_component_logcpm_z_cap5`, equal weighting across cell types and within-cell-type correlation clusters, and tiered versus strict-QC comparisons. The earlier standalone August PCA notebook is superseded.

All code-cell ASTs match the source notebook after removal of conversational comments, markdown and execution outputs. Necessary technical comments remain. Validation and the source fingerprint are recorded in [provenance.json](provenance.json). No analysis was rerun during curation.
