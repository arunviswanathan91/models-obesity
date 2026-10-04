"""Orthogonal-patient precision model: one change from the published (legacy) model.

Variants
- filtered:     same 1,000-read-excluded rows as the legacy fit; only the residual scale changes.
- all_profiles: every profile kept; residual scale depends on allocated reads instead of a cutoff.

The legacy trace is read only. New fits are post hoc sensitivity analyses.
"""
from pathlib import Path
import gc, json, time, traceback
import numpy as np
import pandas as pd
import arviz as az
import bhm_extension as e

EXPECTED_ENGINE_SHA = "e0fb433b32c94c80784924c95538af3e036ee3319969780bed07f1cf5bb7fedd"
VARIANTS = ("filtered", "all_profiles")

def spec(precision_sd):
    return e.ModelSpec("orthogonal_precision", precision=True, patient_mode="orthogonal",
                       precision_sd=precision_sd)

def prepare(cfg, family):
    full, meta, dropped = e.load_frozen(cfg["base"], family, cfg.get("pca_override"),
                                        cfg.get("review_override"), cfg.get("run_overrides", {}).get(family))
    filtered, depth = e.precision_transform(full.loc[~full.low_read].copy())
    everything, depth_all = e.precision_transform(full.copy())
    e.require(depth == depth_all, "Read-depth scaling must come from the same >=1000-read profiles")
    arrays = {"filtered": e.make_arrays(filtered, full), "all_profiles": e.make_arrays(everything, full)}
    e.require(len(arrays["all_profiles"]["y"]) - len(arrays["filtered"]["y"]) == int(full.low_read.sum()),
              "All-profiles variant must add exactly the excluded rows")
    return arrays, meta, dropped, depth

def reports_current(out, trace_sha):
    receipt = out / "report_receipt.json"
    if not receipt.exists():
        return False
    r = json.loads(receipt.read_text())
    return r.get("trace_sha256") == trace_sha and all((out / n).exists() for n in r.get("files", []))

def report(idata, a, model_spec, scope, out, diag, ppc_draws, trace_sha):
    if reports_current(out, trace_sha):
        return pd.read_csv(out / "feature_effects.csv"), pd.read_csv(out / "ppc_summary.csv")
    tab = e.effect_table(idata, a, out, diag["publication_gate"])
    ppc = e.ppc(idata, a, model_spec, scope, out, ppc_draws)
    e.write_json(dict(trace_sha256=trace_sha, files=["feature_effects.csv", "ppc_summary.csv"]),
                 out / "report_receipt.json")
    return tab, ppc

def legacy(cfg, family, a, root):
    out = root / "legacy_readonly"
    audit = out / "legacy_audit.json"
    if audit.exists() and (out / "feature_effects.csv").exists() and (out / "ppc_summary.csv").exists():
        return json.loads(audit.read_text())["diagnostics"]
    idata, diag = e.inspect_legacy(cfg["base"], family, a, out)
    e.ppc(idata, a, e.model_registry()["legacy"], e.SPECS[family]["nu"], out, cfg["ppc_draws"])
    idata.close() if hasattr(idata, "close") else None
    return diag

def eta_summary(idata, a):
    post = idata.posterior["read_log_scale_slope"]
    hdi = az.hdi(idata, var_names=["read_log_scale_slope"], hdi_prob=.95)["read_log_scale_slope"].values
    draws = post.values.reshape(-1, post.shape[-1])
    return pd.DataFrame(dict(global_celltype=a["celltypes"], eta_mean=draws.mean(0),
                             eta_hdi_low=hdi[:, 0], eta_hdi_high=hdi[:, 1], prob_negative=(draws < 0).mean(0),
                             sd_ratio_plus1_read_sd=np.exp(draws.mean(0))))

def compare(family, a, root, variants_done):
    ref_path = root / "legacy_readonly" / "feature_effects.csv"
    if not ref_path.exists():
        return [], []
    ref = pd.read_csv(ref_path).set_index("global_feature").loc[a["features"]]
    ct = pd.Series([a["celltypes"][c] for c in a["f2c"]], index=a["features"])
    rows, ct_rows = [], []
    for variant in variants_done:
        t = pd.read_csv(root / f"orthoprec_{variant}" / "publication" / "feature_effects.csv")
        t = t.set_index("global_feature").loc[a["features"]]
        d = t["mean"] - ref["mean"]
        rows.append(dict(family=family, variant=variant,
                         slope_correlation=float(np.corrcoef(t["mean"], ref["mean"])[0, 1]),
                         mean_absolute_slope_change=float(d.abs().mean()),
                         max_absolute_slope_change=float(d.abs().max()),
                         legacy_mean_slope=float(ref["mean"].mean()), variant_mean_slope=float(t["mean"].mean()),
                         median_interval_width_legacy=float((ref.hdi_high - ref.hdi_low).median()),
                         median_interval_width_variant=float((t.hdi_high - t.hdi_low).median()),
                         hdi_excludes_zero_legacy=int(ref.hdi_excludes_zero.sum()),
                         hdi_excludes_zero_variant=int(t.hdi_excludes_zero.sum()),
                         hdi_status_changes=int((t.hdi_excludes_zero != ref.hdi_excludes_zero).sum()),
                         rope010_selected_legacy=int(ref["rope_0.10_selected"].sum()),
                         rope010_selected_variant=int(t["rope_0.10_selected"].sum())))
        for c in a["celltypes"]:
            m = ct == c
            ct_rows.append(dict(family=family, variant=variant, global_celltype=c,
                                legacy_mean_slope=float(ref.loc[m, "mean"].mean()),
                                variant_mean_slope=float(t.loc[m, "mean"].mean())))
    return rows, ct_rows

def plots(family, a, root, variants_done, figdir):
    import matplotlib.pyplot as plt
    ref_path = root / "legacy_readonly" / "feature_effects.csv"
    if not ref_path.exists() or not variants_done:
        return
    ref = pd.read_csv(ref_path).set_index("global_feature").loc[a["features"]]
    fig, axes = plt.subplots(1, len(variants_done), figsize=(4.4 * len(variants_done), 4.2), squeeze=False)
    for ax, variant in zip(axes[0], variants_done):
        t = pd.read_csv(root / f"orthoprec_{variant}" / "publication" / "feature_effects.csv")
        t = t.set_index("global_feature").loc[a["features"]]
        lim = max(ref["mean"].abs().max(), t["mean"].abs().max()) * 1.15
        ax.scatter(ref["mean"], t["mean"], s=9, alpha=.6)
        ax.plot([-lim, lim], [-lim, lim], color="grey", lw=.8)
        ax.axhspan(-.10, .10, color="grey", alpha=.08)
        ax.set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel="Published model: BMI slope",
               ylabel=f"Precision model ({variant}): BMI slope", title=family)
    fig.tight_layout()
    fig.savefig(figdir / f"{family}_slopes_vs_published.pdf")
    fig.savefig(figdir / f"{family}_slopes_vs_published.svg")
    plt.close(fig)

def run(cfg):
    e.check_environment()
    e.require(e.sha(e.__file__) == EXPECTED_ENGINE_SHA, "Engine differs from the audited extension engine")
    study = Path(cfg["base"]) / "normalized_state_bhm_extensions_v1" / cfg["study_id"]
    study.mkdir(parents=True, exist_ok=True)
    scientific = {k: v for k, v in cfg.items() if k not in ("cores", "stop_on_failed_geometry")}
    lock = dict(config=scientific, engine_sha=e.sha(e.__file__), driver_sha=e.sha(__file__),
                model=spec(cfg["precision_slope_sd"]).__dict__, priors=e.PRIORS,
                design="Legacy orthogonal patient effect plus read-dependent residual scale; one change from the published model")
    lp = study / "orthoprec_lock.json"
    if lp.exists():
        e.require(json.loads(lp.read_text()) == lock, "Study settings or code changed; use a new STUDY_ID")
    else:
        e.write_json(lock, lp)

    status, comparisons, celltype_rows, etas = [], [], [], []
    figdir = study / "figures"
    figdir.mkdir(exist_ok=True)
    model_spec = spec(cfg["precision_slope_sd"])
    stop = False
    for fi, family in enumerate(cfg["families"]):
        if stop:
            break
        root = study / family
        scope = e.SPECS[family]["nu"]
        try:
            arrays, meta, dropped, depth = prepare(cfg, family)
        except Exception as exc:
            status.append(dict(family=family, model="all", status="blocked", error=str(exc)))
            e.csv(pd.DataFrame(status), study / "orthoprec_status.csv")
            print(f"{family}: blocked — {exc}", flush=True)
            continue
        root.mkdir(parents=True, exist_ok=True)
        e.csv(dropped, root / "excluded_profiles.csv")
        e.write_json(depth, root / "read_depth_scaling.json")
        e.write_json(e.environment(), root / "environment.json")

        row = dict(family=family, model="legacy_readonly", status="running", error="")
        try:
            print(f"\n{family} / published model: reading saved trace", flush=True)
            diag = legacy(cfg, family, arrays["filtered"], root)
            row.update(status="passed" if diag["publication_gate"] else "diagnostics_failed", **diag)
        except Exception as exc:
            row.update(status="failed", error=str(exc))
            (root / "legacy_error.txt").write_text(traceback.format_exc())
        status.append(row)
        e.csv(pd.DataFrame(status), study / "orthoprec_status.csv")
        print(f"{family} / published model: {row['status']} {row['error']}", flush=True)

        done = []
        for vi, variant in enumerate(VARIANTS):
            a = arrays[variant]
            out = root / f"orthoprec_{variant}" / "publication"
            row = dict(family=family, model=f"orthoprec_{variant}", status="running", error="",
                       observations=len(a["y"]), reused_trace=(out / "posterior_trace.nc").exists())
            idata = None
            start = time.time()
            try:
                print(f"\n{family} / {variant}: " + ("checking saved trace" if row["reused_trace"] else "new full fit"),
                      flush=True)
                idata, diag = e.fit(a, model_spec, scope, out, "publication",
                                    seed=20261001 + 10 * fi + vi, cores=cfg["cores"])
                trace_sha = json.loads((out / "trace_receipt.json").read_text())["sha256"]
                report(idata, a, model_spec, scope, out, diag, cfg["ppc_draws"], trace_sha)
                eta = eta_summary(idata, a)
                eta.insert(0, "variant", variant)
                eta.insert(0, "family", family)
                e.csv(eta, out / "read_scale_slopes.csv")
                etas.append(eta)
                row.update(status="passed" if diag["publication_gate"] else "diagnostics_failed", **diag)
                done.append(variant)
                if cfg.get("stop_on_failed_geometry", True) and not diag["geometry_ok"]:
                    stop = True
            except Exception as exc:
                row.update(status="failed", error=str(exc))
                out.mkdir(parents=True, exist_ok=True)
                (out / "orthoprec_error.txt").write_text(traceback.format_exc())
            finally:
                if idata is not None and hasattr(idata, "close"):
                    idata.close()
                row["minutes"] = round((time.time() - start) / 60, 1)
                status.append(row)
                e.csv(pd.DataFrame(status), study / "orthoprec_status.csv")
                print(f"{family} / {variant}: {row['status']} {row['error']}", flush=True)
                gc.collect()
            if stop:
                print("Sampling geometry failed; queue stopped so compute is not wasted. Send orthoprec_status.csv.",
                      flush=True)
                break
        rows, ct_rows = compare(family, arrays["filtered"], root, done)
        comparisons += rows
        celltype_rows += ct_rows
        plots(family, arrays["filtered"], root, done, figdir)

    if comparisons:
        e.csv(pd.DataFrame(comparisons), study / "comparison_vs_published.csv")
        e.csv(pd.DataFrame(celltype_rows), study / "celltype_slopes_vs_published.csv")
    if etas:
        e.csv(pd.concat(etas, ignore_index=True), study / "read_scale_slopes.csv")
    ppc_rows = []
    for family in cfg["families"]:
        for model in ["legacy_readonly", *[f"orthoprec_{v}/publication" for v in VARIANTS]]:
            p = study / family / model / "ppc_summary.csv"
            if p.exists():
                t = pd.read_csv(p)
                t.insert(0, "model", model.split("/")[0])
                t.insert(0, "family", family)
                ppc_rows.append(t)
    if ppc_rows:
        e.csv(pd.concat(ppc_rows, ignore_index=True), study / "ppc_comparison.csv")
    result = pd.DataFrame(status)
    e.write_json(dict(counts=result.status.value_counts().to_dict(), stopped_early=stop,
                      interpretation="Sensitivity analysis against the published model; no automatic winner."),
                 study / "orthoprec_summary.json")
    print(f"\nFinished. Review {study}/orthoprec_status.csv and comparison_vs_published.csv", flush=True)
    return result
