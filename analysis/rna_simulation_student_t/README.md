# Cellular RNA simulation (Student-t intercept vs orthogonal patient-effect model)

Simulation on the fine immune design reported in Figure 2 and Figures S1–S2. Truth values are posterior means from the primary fine immune Student-t fit.

| File | Use in the paper |
|---|---|
| `k1_validation_simulation.ipynb` | Study set-up, data generator and shard runner (study `k1_validation_fine4_mount5_cpu_v1`) |
| `simulation_finish_student_t.ipynb` | Student-t intercept fits and verified aggregation for the same study and seeds |
| `simulation_part2_published_group_grid.ipynb` | Orthogonal patient-effect model, WHO BMI-group model and effect-size grid |
| `sim_ols_feature_benchmark.ipynb` | Feature-level OLS and Benjamini–Hochberg benchmark |
| `celltype_composite_benchmark_cell.py` | Cell-type composite benchmark |

The study configuration also contains a one-factor (K = 1) variant; it was not used in the paper.
