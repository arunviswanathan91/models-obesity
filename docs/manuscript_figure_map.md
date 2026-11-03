# Manuscript figure numbering

Current numbering follows the manuscript and presentation decks supplied on 3 November 2026: Figures 1–7 and supplementary Figures S1–S24. The supplementary file is named `Supplementary_Figures_S1-S19.pptx`, but its contents run through S24.

The [machine-readable figure reference](../web/manuscript-figures.json) records the reference-file SHA-256 hashes, current collection titles and panel mapping. The [data manifest](../web/manuscript-manifest.json) retains the original numerical release and source fingerprints.

## Changes from the archived figure package

| Content | Archived numbering | Current numbering |
| --- | --- | --- |
| Named continuous RNA examples | S23a–d | 3b–e |
| Continuous/categorical agreement | S13a–c | S10d–f |
| RNA–protein sensitivity | S14 | S11 |
| RNA–protein model diagnostics | S15 | S12 |
| RNA–protein eligibility and direction agreement | S16 | S13a,b |
| Observed RNA–protein examples | S17 | S13c,d |
| PTM landscape and diagnostic coverage | S18 | S14 |
| STAT3 and ASAH1 feature estimates | S19 | S15 |
| STAT3/ASAH1 sensitivities and bootstrap | S20 | S16 |
| Cellular RNA workflow | S21 | S17a |
| Paired molecular-layer workflows | S22a,b | S17b,c |
| PTM calibration coverage | S24a | S18 |
| PTM calibration bias | S24b | S19 |
| PTM calibration RMSE | S24c | S20 |
| SBC ranks | S24d | S21a |
| Simulation completion | S24e | S21b |
| PTM sensitivity overview/classifications | S25a,b | S22a,b |
| First four PTM sensitivity comparisons | S25c–f | S23a–d |
| Last four PTM sensitivity comparisons | S25g–j | S24a–d |

## Current panel references

| Figure/panel | Content |
| --- | --- |
| 3a | RNA heatmap |
| 3b–e | Named continuous RNA examples, formerly Figure S23 |
| 6c | All 14 conditional phospho candidates |
| 6d | Top 20 of 67 conditional glyco candidates |
| S14a | Protein-conditioned PTM landscape |
| S14b | Candidate diagnostic and sensitivity coverage |
| S15a | STAT3 Y705/S727 estimates |
| S15b | ASAH1 glycofeature estimates and identities |
| S16a | STAT3 continuous sensitivity |
| S16b | STAT3 categorical contrasts |
| S16c | STAT3 and ASAH1 bootstrap benchmarks |
| S17a | Cellular RNA workflow |
| S17b | RNA–protein workflow |
| S17c | Protein–PTM workflow |
| S18 | PTM interval coverage |
| S19 | PTM bias |
| S20 | PTM RMSE |
| S21a | SBC ranks |
| S21b | Simulation completion |
| S22a | PTM sensitivity overview |
| S22b | Interval classifications |
| S23a | Fixed plex versus primary |
| S23b | Stage subset versus primary |
| S23c | Stage adjustment versus matched subset |
| S23d | Purity subset versus primary |
| S24a | Purity adjustment versus matched subset |
| S24b | Mutation adjustment versus primary |
| S24c | Weight-loss exclusion versus primary |
| S24d | Adenosquamous exclusion versus primary |

## Data verification and deployment

The uploaded October ZIP matches the archive SHA-256 in the data manifest. All 18 source tables shared with that archive (5,241 rows) matched the original Supabase records. The two ASAH1 CSVs are separate archived sources and are not included in this ZIP. Public data values were checked against the stored original rows after the documented type conversion. Figures 2–7 in the updated main deck retain 38 exact SVG plot/legend matches to the archive.

Historical SVG/CSV filenames are source provenance, not current figure numbers. Do not rename or rewrite those inputs to implement a caption change. Figure 3 now contains the named examples in panels b–e as well as the heatmap in panel a.

The Atlas reads collection titles and descriptions from Supabase when loaded. The metadata update in [align_manuscript_figure_numbers_20261103.sql](../database/align_manuscript_figure_numbers_20261103.sql) updates those labels and figure recipe references. It preserves dataset IDs, result rows, effect estimates, source hashes, publication flags and renderer functions. The release preparation script uses the same current labels for future imports.
