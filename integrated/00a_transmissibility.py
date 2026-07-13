#!/usr/bin/env python3
"""00a_transmissibility.py -- compute per-gene observed transmissibility from
str_omic (data_download/from_source), replacing the committed gene_intrinsic atlas.

Definition (faithful reimplementation of the gene-intrinsic label):
  observed_transmissibility(gene) = mean over amplified, unique (caseid, gene) rows
  of a binary protein-positive call, where
    - amplified  : R4 filter = cn >= round(ploidy)+1 AND cn >= 3   (high-level gain)
    - protein-positive : prot.rel.all > PROT_POS_THR   (pan-cohort ECDF protein rank)
  Genes retained with >= MIN_AMP_CASES amplified cases.

This reproduces the committed GI atlas at Spearman rho ~= 0.97 (MAE ~0.05) and the
same ~0.40 base rate; the residual is the exact ECDF/protein-positive definition in
the GI feather. Emits integrated/tables/transmissibility_atlas.csv.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg

PROT_POS_THR = 0.83     # prot.rel.all rank above which a tumour is protein-positive
MIN_AMP_CASES = 20      # minimum amplified cases to retain a gene (GI primary support)
PROT_COL = "prot.rel.all"


def _r4_amplified(df, ploidy):
    """R4 high-level amplification: cn >= round(ploidy)+1 AND cn >= 3."""
    p = df.merge(ploidy[["caseid", "ploidy_continuous"]], on="caseid", how="left")
    pr = p["ploidy_continuous"].round().fillna(2)
    return p[(p.cn >= pr + 1) & (p.cn >= 3)].copy()


def build():
    cols = ["gene", "caseid", "tumor_code", "cn", "cn_adjusted", PROT_COL]
    df = pd.read_csv(cfg.PATHS["str_omic"], usecols=cols)
    ploidy = pd.read_csv(cfg.TAB / "ploidy_table.csv")
    amp = _r4_amplified(df, ploidy).dropna(subset=[PROT_COL])
    amp = amp.drop_duplicates(["caseid", "gene"])
    amp["target"] = (amp[PROT_COL] > PROT_POS_THR).astype(float)
    g = amp.groupby("gene")["target"].agg(observed_transmissibility="mean",
                                           n_amplified_cases="size").reset_index()
    g = g[g["n_amplified_cases"] >= MIN_AMP_CASES].reset_index(drop=True)
    base = float(amp["target"].mean())
    return g, base


def main():
    atlas, base = build()
    out = cfg.DIR_TAB / "transmissibility_atlas.csv"
    atlas.to_csv(out, index=False)
    print(f"transmissibility_atlas.csv: {len(atlas)} genes, "
          f"base_rate={base:.4f}, median_obs={atlas.observed_transmissibility.median():.4f}")
    # provenance check vs committed GI atlas, if present
    try:
        gi = pd.read_csv(cfg.GI_CROSSCHECK["atlas"]).set_index("gene")["observed_transmissibility"]
        m = atlas.set_index("gene")["observed_transmissibility"].to_frame("new").join(gi.rename("gi"), how="inner")
        rho = m["new"].corr(m["gi"], method="spearman")
        mae = (m["new"] - m["gi"]).abs().mean()
        print(f"vs committed GI atlas: n_shared={len(m)} spearman={rho:.4f} MAE={mae:.4f}")
    except Exception as e:
        print(f"(GI cross-check skipped: {e})")
    return atlas


if __name__ == "__main__":
    main()
