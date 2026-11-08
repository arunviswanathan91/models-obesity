# Reference construction, deconvolution and compartment QC

## Reference construction and BayesPrism

- [immune_reference_deconvolution.ipynb](immune_reference_deconvolution.ipynb): CD45+ reference QC, SingleR annotation, immune subtype/TAM refinement and coarse/fine BayesPrism fits.
- [nonimmune_reference_deconvolution.ipynb](nonimmune_reference_deconvolution.ipynb): whole-tumour reference QC, CopyKAT classification, epithelial/stromal annotation, nonimmune BayesPrism fitting and component-expression extraction.

Open each notebook in a Colab Python runtime, mount Drive in the first cell, configure the source paths, and execute retained cells in order. The immune workflow needs the CD45+ count matrix under `BayesPrism/scRNA_data`, three bulk count files under `BayesPrism/data`, and reference annotation packages. The nonimmune workflow needs the GSE242230 matrix/barcode/feature files under `BayesPrism/GSE242230` and the bulk files defined in its fitting cell. The saved-output directories default to `/content/drive/MyDrive/BayesPrism`.

Run the immune and nonimmune fitting workflows before the consolidated expression-extraction cell at the end of the nonimmune notebook. Configure its three `bp_files` paths to the completed objects:

| Compartment | Saved RDS relative to `BayesPrism` |
| --- | --- |
| Immune coarse | `results_annotation/bayesprism_results.rds` |
| Immune fine | `results_annotation/bayesprism_results_FINE_TAM.rds` |
| Nonimmune | `improved_deconvolution_results/bayesprism_object.rds` |

The source extraction cell retains its original local RDS aliases under `/content`; these must be mapped to the saved files above or populated before extraction. These are the same RDS objects resolved by the downstream compartment-QC notebook.

Run the annotation and refinement cells in order; later labels depend on earlier steps.

Source details are in [upstream_provenance.json](upstream_provenance.json). Historical R package versions are not pinned.

## Downstream compartment QC and PCA

[compartment_qc_and_pca.ipynb](compartment_qc_and_pca.ipynb) contains the final compartment audit, endpoint comparison and weighted publication PCA analysis. It reads completed BayesPrism deconvolution objects and expression components; it does not perform deconvolution or fit the BMI model.

The analysis validates theta and component-expression identities, reconstructs the 127-patient cohort, checks signature mapping, missingness and zero counts, compares score definitions, evaluates empirical dependence, and constructs tiered and strict-QC weighted PCA summaries.

Run the compartment audit for `immune_course`, `immune_fine` and `non_immune` separately using `TARGET_COMPARTMENT`. The combined endpoint and PCA cells require saved audit outputs for all three compartments. Configure `BASE`, `BAYESPRISM_ROOT` and the RDS override in the initial configuration cell to match the mounted input directory. Input RDS objects, patient matrices, clinical data and the frozen signature catalog are required and are not bundled here.

PCA uses `within_component_logcpm_z_cap5`, equal weighting across cell types and within-cell-type correlation clusters, and tiered versus strict-QC comparisons.

Source hashes are in [provenance.json](provenance.json).
