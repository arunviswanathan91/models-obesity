"""Build aggregate BayesPrism fraction plots from the September compartment audit.

Usage: python theta_bmi.py cohort.csv coarse/patient_qc.csv fine/patient_qc.csv nonimmune/patient_qc.csv
Writes JSON to stdout; input files contain study IDs, output contains aggregate statistics only.
"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, kruskal, spearmanr

GROUPS = ["normal_weight", "overweight", "obese"]


def summarize(data):
    grouped = [data.loc[data.group == g, "theta"] for g in GROUPS]
    rho, continuous_p = spearmanr(data.bmi, data.theta)
    h, group_p = kruskal(*grouped)
    bins = pd.qcut(data.bmi, 8, duplicates="drop")
    trend = []
    for _, sample in data.groupby(bins, observed=True):
        low, median, high = sample.theta.quantile([0.25, 0.5, 0.75])
        trend.append(dict(bmi=sample.bmi.median(), n=len(sample), low=low, median=median, high=high))
    groups = []
    for name, values in zip(GROUPS, grouped):
        low, median, high = values.quantile([0.25, 0.5, 0.75])
        groups.append(dict(name=name, n=len(values), low=low, median=median, high=high))
    return dict(n=len(data), rho=rho, continuous_p=continuous_p, h=h, group_p=group_p,
                trend=trend, groups=groups)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["cohort", "coarse", "fine", "nonimmune"]:
        parser.add_argument(name)
    args = parser.parse_args()
    cohort = pd.read_csv(args.cohort, usecols=["base_sample_id", "bmi", "bmi_group"])
    cohort["group"] = pd.cut(cohort.bmi, [18.5, 25, 30, np.inf], labels=GROUPS, right=False)
    assert cohort.base_sample_id.is_unique and cohort.group.eq(cohort.bmi_group).all()
    assert cohort.group.value_counts().reindex(GROUPS).tolist() == [52, 57, 18]
    output = []
    for family, path in [("Coarse immune", args.coarse), ("Fine immune", args.fine), ("Nonimmune", args.nonimmune)]:
        qc = pd.read_csv(path, usecols=["base_sample_id", "CellType", "theta"])
        assert not qc.duplicated(["base_sample_id", "CellType"]).any()
        assert qc.theta.between(0, 1).all()
        for component, rows in qc.groupby("CellType"):
            assert set(rows.base_sample_id) == set(cohort.base_sample_id)
            joined = rows.merge(cohort, on="base_sample_id", validate="one_to_one")
            output.append(dict(family=family, component=component, **summarize(joined)))
    for test in ["continuous", "group"]:
        adjusted = false_discovery_control([r[test + "_p"] for r in output], method="bh")
        for row, q in zip(output, adjusted):
            row[test + "_q"] = q
    print(json.dumps(output, allow_nan=False))


if __name__ == "__main__":
    main()
