#!/usr/bin/env python3
"""d2_cross_cohort.py -- Supplementary Table S5: observed cross-cohort transmissibility
transfer, computed from data_download/from_source (str_omic) in integrated/.

For each CPTAC cohort we compute per-gene observed transmissibility (fraction of
amplified cases with elevated tissue-referenced protein) using the same R4 amplification
filter and prot.rel.all>0.83 protein-positive call as 00a. We then correlate the
per-gene observed transmissibility between every cohort pair that shares >=20 genes
(Spearman), quantifying how well the observed cascade transfers across lineages.

Emits tables/m6_observed_cross_cohort.csv (pair, n_genes, rho, p), replacing the static
committed CSV.
"""
import sys
from pathlib import Path
from itertools import combinations
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg

MIN_CASES = 10       # amplified cases per gene per cohort
MIN_SHARED = 20      # shared genes for a cohort-pair correlation
PROT_POS = 0.83      # prot.rel.all ECDF threshold = protein-positive (matches 00a)


def per_cohort_transmissibility():
    so = pd.read_csv(cfg.PATHS["str_omic"],
                     usecols=["gene", "caseid", "tumor_code", "cn", "prot.rel.all"])
    ploidy = pd.read_csv(cfg.DATA_ROOT / "out" / "tables" / "ploidy_table.csv")
    so = so.merge(ploidy[["caseid", "ploidy_continuous"]], on="caseid", how="left")
    so["pr"] = so.ploidy_continuous.round().fillna(2)
    r4 = so[(so.cn >= so.pr + 1) & (so.cn >= 3)].dropna(subset=["prot.rel.all"]).copy()
    r4["target"] = (r4["prot.rel.all"] > PROT_POS).astype(float)
    per = {}
    for ct, gg in r4.groupby("tumor_code"):
        lab = gg.drop_duplicates(["caseid", "gene"]).groupby("gene")["target"].agg(["mean", "size"])
        per[ct] = lab[lab["size"] >= MIN_CASES]["mean"]
    return pd.DataFrame(per)


def main():
    mat = per_cohort_transmissibility()
    rows = []
    for a, b in combinations(sorted(mat.columns), 2):
        m = mat[[a, b]].dropna()
        if len(m) < MIN_SHARED:
            continue
        rho, p = spearmanr(m[a], m[b])
        rows.append({"pair": f"{a}-{b}", "n_genes": len(m), "rho": round(rho, 3), "p": p})
    out = pd.DataFrame(rows).sort_values("rho", ascending=False)
    out.to_csv(cfg.DIR_TAB / "m6_observed_cross_cohort.csv", index=False)
    print(f"m6_observed_cross_cohort.csv: {len(out)} cohort pairs")
    print(out.to_string(index=False))
    return out


if __name__ == "__main__":
    main()
