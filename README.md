# models-obesity
## Manuscript supporting data

The October 2026 manuscript release adds 20 aggregate result collections (5,381 rows) with full feature identities, paired ASAH1 quantities, 192 matched PTM sensitivity comparisons and metrics from 250 calibration datasets. `web/manuscript-manifest.json` records source names, hashes and row counts. `scripts/prepare_manuscript_release.py` prepares unpublished import batches from the supplied result archive and ASAH1 tables; counts and values must be verified before publication.

The atlas generates these plots on demand in its existing Python worker from Supabase rows. It does not load archived manuscript SVGs. The original `draw_heatmap` function is retained unchanged; manuscript adapters supply the two PTM effect columns. Calibration and sensitivity adapters reuse the supplied extension plotting code, palette and interval conventions. Short PTM labels retain biological identities; peptide sequences distinguish otherwise identical glycoforms. Stage and purity adjustment compare against their own unadjusted subsets.

Renderer checks: `python -m unittest discover -s tests -p 'test_manuscript_renderer.py'`. Downloads continue to use the corresponding-author request dialog.
