from google.colab import drive, files
drive.mount("/content/drive")

import numpy as np, pandas as pd, arviz as az
from pathlib import Path

BASE = Path("/content/drive/MyDrive/new_dysfunction")
NEW_ROOT = BASE / "normalized_state_bhm_ri_hc_v1"
PRIM_ROOT = BASE / "normalized_state_bhm_lowread_sensitivity_v1"
RUNS = {
    "immune_coarse": "v2p2_g1_fixed_ctnu__immune_coarse__lowread1000",
    "immune_fine": "v2p1c_g1_fixed_gnu__immune_fine__lowread1000",
    "nonimmune": "v2p1c_g1_fixed_gnu__nonimmune__lowread1000",
}
OUT = Path("/content/within_celltype_contrasts"); OUT.mkdir(exist_ok=True)

def contrasts(trace_path):
    post = az.from_netcdf(trace_path).posterior[["feature_bmi_slope", "celltype_bmi_slope"]].load()
    feats = np.asarray(post.feature.values, str)
    cts = np.asarray(post.celltype.values, str)
    parent = np.array([f.rsplit("||", 1)[0] for f in feats])
    missing = set(parent) - set(cts)
    assert not missing, f"feature prefixes not in celltype coord: {sorted(missing)[:3]}"
    m = post["celltype_bmi_slope"].sel(celltype=parent).rename({"celltype": "feature"})
    m = m.assign_coords(feature=feats)
    d = (post["feature_bmi_slope"] - m).rename("contrast")
    h = az.hdi(d.to_dataset(name="contrast"), hdi_prob=0.95)["contrast"].sel(feature=feats).transpose("feature", "hdi").values
    s = np.asarray(d.transpose("feature", "chain", "draw")).reshape(len(feats), -1)
    p_pos = (s > 0).mean(axis=1)
    return pd.DataFrame({"global_feature": feats, "celltype": parent,
                         "mean": s.mean(axis=1), "hdi_low": h[:, 0], "hdi_high": h[:, 1],
                         "p_positive": p_pos, "lfsr": np.minimum(p_pos, 1 - p_pos)})

summary = []
for fam, run in RUNS.items():
    new = contrasts(NEW_ROOT / f"{run}__ri_hc_t" / "publication" / "posterior_trace.nc")
    prim = contrasts(PRIM_ROOT / run / "publication" / "posterior_trace.nc")
    tab = new.merge(prim, on=["global_feature", "celltype"], suffixes=("_new", "_prim"))
    for k in ("new", "prim"):
        tab[f"excl_{k}"] = (tab[f"hdi_low_{k}"] > 0) | (tab[f"hdi_high_{k}"] < 0)
    tab["width_ratio"] = (tab.hdi_high_new - tab.hdi_low_new) / (tab.hdi_high_prim - tab.hdi_low_prim)
    tab.to_csv(OUT / f"{fam}_contrasts.csv", index=False)
    summary.append({"family": fam, "n_features": len(tab),
                    "excl_new": int(tab.excl_new.sum()), "excl_prim": int(tab.excl_prim.sum()),
                    "lfsr05_new": int((tab.lfsr_new <= 0.05).sum()),
                    "lfsr05_prim": int((tab.lfsr_prim <= 0.05).sum()),
                    "expected_if_null": round(0.05 * len(tab), 1),
                    "median_width_ratio": round(float(tab.width_ratio.median()), 3),
                    "corr_mean": round(float(tab.mean_new.corr(tab.mean_prim)), 3)})
    print(f"\n{fam}: contrasts excluding zero (ri_hc_t)")
    cols = ["global_feature", "mean_new", "hdi_low_new", "hdi_high_new", "lfsr_new", "excl_prim"]
    print(tab[tab.excl_new].sort_values("mean_new")[cols].to_string(index=False))

summary = pd.DataFrame(summary)
summary.to_csv(OUT / "summary.csv", index=False)
print("\n", summary.to_string(index=False))

import shutil
shutil.make_archive("/content/within_celltype_contrasts", "zip", OUT)
files.download("/content/within_celltype_contrasts.zip")
