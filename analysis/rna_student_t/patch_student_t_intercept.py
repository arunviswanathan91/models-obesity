from pathlib import Path

p = Path("/content/ri_hc_engine.py")
s = p.read_text()

edits = [
    ('NOTEBOOK_VERSION = "ri-hc-celltype-intercept-v1"',
     'NOTEBOOK_VERSION = "ri-hc-t-celltype-intercept-v1"'),
    ('RUN_NAME = f"{PRIMARY_RUN_NAME}__ri_hc"',
     'RUN_NAME = f"{PRIMARY_RUN_NAME}__ri_hc_t"'),
    ('''        celltype_intercept_w = pm.Normal("celltype_intercept_w", w_mean, celltype_intercept_sigma[None, :],
                                         dims=("patient", "celltype"))''',
     '''        intercept_nu_minus_two = pm.Exponential("intercept_nu_minus_two", lam=NU_MINUS_TWO_RATE)
        intercept_nu = pm.Deterministic("intercept_nu", intercept_nu_minus_two + 2.0)
        celltype_intercept_w = pm.StudentT("celltype_intercept_w", nu=intercept_nu, mu=w_mean,
                                           sigma=celltype_intercept_sigma[None, :],
                                           dims=("patient", "celltype"))'''),
    ('"log_sigma_feature_scale", "log_sigma_raw", "nu_minus_two"}',
     '"log_sigma_feature_scale", "log_sigma_raw", "nu_minus_two", "intercept_nu_minus_two"}'),
]
for old, new in edits:
    assert s.count(old) == 1, f"pattern not found exactly once:\n{old}"
    s = s.replace(old, new)
p.write_text(s)
print("engine patched: Student-t patient x cell-type intercept, outputs under __ri_hc_t")
