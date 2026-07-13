#!/usr/bin/env python3
"""d1_threshold_sensitivity.py -- Supplementary Table S3: amplification-threshold
sensitivity of the nomination funnel, computed in integrated/ from data_download/from_source.

Re-runs the recurrent-amplicon + co-elevation pipeline (00c logic) and the transmitted /
surviving-antigen counts at two amplification thresholds (cn_adjusted >= 1.4 the paper
basis, and >= 2.0 a stringent alternative), so a reader can see the funnel is not an
artefact of the 1.4 cutoff. Every count is recomputed from str_omic + the integrated
transmissibility atlas + the 22 nominated antigens; nothing is a static literal.

Emits tables/d1_threshold_sensitivity.csv (metric, thresh_1.4, thresh_2.0).
"""
import sys, importlib
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg
import cnt_shared as ish
import importlib.util

def _load_00c():
    spec = importlib.util.spec_from_file_location("m00c", Path(__file__).resolve().parent / "00c_amplicons.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

TRANSMIT_MIN = 0.40


def funnel_at(threshold, m00c, atlas, nominated):
    """Return the stage counts at a given amplification threshold."""
    orig = cfg.AMP_THRESHOLD
    cfg.AMP_THRESHOLD = threshold
    try:
        d = ish.load_str_omic().dropna(subset=["cn_adjusted","cytogenetic_location","prot.rel.tissue"]).copy()
        d["amplified"] = (d.cn_adjusted >= threshold).astype(int)
        d["prot_high"] = (d["prot.rel.tissue"] > cfg.REL_TISSUE_HI).astype(int)
        br, band_samp = m00c.recurrent_amplicons(d)
        rec = br[br.recurrent].copy()
        rec_keys = set(zip(rec.tumor_code, rec.cytogenetic_location))
        co = m00c.coelevation(d, band_samp, rec_keys)
        co_sig = co[co.fdr < cfg.FDR_ALPHA]
        n_recurrent = len(rec)
        n_coelev = co_sig.gene.nunique()
        # transmitted genes = observed transmissibility >= floor in the integrated atlas
        n_transmit = int((atlas.observed_transmissibility >= TRANSMIT_MIN).sum())
        # co-elevated AND transmitted
        transmit_genes = set(atlas.loc[atlas.observed_transmissibility >= TRANSMIT_MIN, "gene"])
        n_coelev_transmit = len(set(co_sig.gene) & transmit_genes)
        # of the 22 nominated antigens, how many are co-elevated at this threshold
        n_surviving = len(set(nominated) & set(co_sig.gene) & transmit_genes)
        return dict(recurrent=n_recurrent, coelev=n_coelev, transmit=n_transmit,
                    coelev_transmit=n_coelev_transmit, surviving=n_surviving)
    finally:
        cfg.AMP_THRESHOLD = orig


def main():
    m00c = _load_00c()
    atlas = ish.load_atlas()
    nominated = pd.read_csv(cfg.DIR_TAB / "adc_target_antigens.csv").antigen.unique().tolist()
    n_nom = len(nominated)
    r14 = funnel_at(1.4, m00c, atlas, nominated)
    r20 = funnel_at(2.0, m00c, atlas, nominated)
    rows = [
        ("recurrent amplicons", r14["recurrent"], r20["recurrent"]),
        ("co-elevated genes (FDR<0.1)", r14["coelev"], r20["coelev"]),
        ("transmitted genes (obs>=0.40)", r14["transmit"], r20["transmit"]),
        ("co-elevated & transmitted", r14["coelev_transmit"], r20["coelev_transmit"]),
        (f"of {n_nom} nominated antigens surviving", r14["surviving"], r20["surviving"]),
    ]
    out = pd.DataFrame(rows, columns=["metric", "thresh_1.4", "thresh_2.0"])
    out.to_csv(cfg.DIR_TAB / "d1_threshold_sensitivity.csv", index=False)
    print(f"d1_threshold_sensitivity.csv ({n_nom} nominated antigens):")
    print(out.to_string(index=False))
    return out


if __name__ == "__main__":
    main()
