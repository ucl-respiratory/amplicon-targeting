# =============================================================================
# _mechanism.py  --  shared dosage-transmission / attenuation computation
# -----------------------------------------------------------------------------
# Foundation for the mechanism figures (7-9,13). Computes, from the wide feature
# table, per-(gene,tumor_type) Pearson correlations of mRNA and protein with
# copy number, then a pan-cancer per-gene estimate by Fisher-z averaging across
# types (weighted by n-3). Cached so the mechanism stages share one computation.
#
#   transmission   = corr(cn, rna)      [CN -> mRNA]      (raw copy number)
#   responsiveness = corr(cn, prot)     [CN -> protein]
#   attenuation    = transmission - responsiveness
#
# Matches the manuscript computation: raw `cn` (not ploidy-adjusted); one joint
# non-NA mask over cn/rna/prot; per-type correlations carry no min-n filter, but
# the pan-cancer per-gene aggregation keeps only per-type estimates with n>=20
# and averages (Fisher-z, weighted by n-3) genes seen in >=2 types.
# =============================================================================
import functools
import numpy as np, pandas as pd
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cnt_io

MIN_N_PER_TYPE = 20   # per-type estimates with fewer tumours dropped pre-aggregation
MIN_TYPES = 2         # a per-gene estimate needs >=2 tumor types (paper genes have all 6)

def _pearson_by_gene(sub, xcol, ycol):
    """Per-gene Pearson r of ycol on xcol across tumours (one type). No min-n."""
    g = sub.groupby("gene")
    n = g.size()
    sx = g[xcol].sum(); sy = g[ycol].sum()
    sxx = g[xcol].apply(lambda t: (t**2).sum())
    syy = g[ycol].apply(lambda t: (t**2).sum())
    sxy = g.apply(lambda t: (t[xcol]*t[ycol]).sum(), include_groups=False)
    cov = sxy - sx*sy/n
    vx = sxx - sx**2/n; vy = syy - sy**2/n
    return cov/np.sqrt(vx*vy), n

@functools.lru_cache(maxsize=1)
def per_type_attenuation():
    """Per (gene,tumor_type): r_cn_rna, r_cn_prot, n. Uses RAW cn, joint dropna."""
    df = cnt_io.load_str_omic()
    d = df.dropna(subset=["cn", "rna", "prot"]).copy()
    rows = []
    for tc, sub in d.groupby("tumor_code"):
        s = sub[["gene", "cn", "rna", "prot"]]
        r_rna, n  = _pearson_by_gene(s, "cn", "rna")
        r_prot, _ = _pearson_by_gene(s, "cn", "prot")
        rows.append(pd.DataFrame({"gene": r_rna.index, "r_cn_rna": r_rna.values,
                                  "r_cn_prot": r_prot.values, "n": n.values,
                                  "tumor_code": tc}))
    return pd.concat(rows, ignore_index=True)

def _fisher_avg(r, n):
    r = np.clip(np.asarray(r, float), -0.999, 0.999)
    z = np.arctanh(r); w = np.asarray(n, float) - 3
    return np.tanh((z*w).sum()/w.sum())

@functools.lru_cache(maxsize=1)
def per_gene_attenuation():
    """Pan-cancer per-gene cn_rna_corr / cn_prot_corr / attenuation."""
    pt = per_type_attenuation().dropna(subset=["r_cn_rna", "r_cn_prot"])
    pt = pt[pt.n >= MIN_N_PER_TYPE]
    agg = []
    for gene, sub in pt.groupby("gene"):
        if len(sub) < MIN_TYPES: continue
        r_rna  = _fisher_avg(sub.r_cn_rna.values,  sub.n.values)
        r_prot = _fisher_avg(sub.r_cn_prot.values, sub.n.values)
        agg.append((gene, r_rna, r_prot, r_rna - r_prot, len(sub), int(sub.n.sum())))
    return pd.DataFrame(agg, columns=["gene","cn_rna_corr","cn_prot_corr",
                                      "attenuation","n_types","n_samples_total"])
