# =============================================================================
# _features.py  --  gene-level CORUM + regulatory features for the mechanism
#                   variance decomposition (Figs 9, 11).
# -----------------------------------------------------------------------------
# CORUM complex features (from data_download stage 05 humanComplexes.txt):
#   is_complex_subunit, corum_n_complexes, corum_max_complex_size
# Regulatory features:
#   mean_promoter_meth, meth_rna_corr  -- from the `meth` column of str_omic
#   log_gene_length, gene_density_1mb  -- from Ensembl gene coordinates (cached)
#
# The Ensembl lookup is cached to OUT/tables/_ensembl_coords.parquet so the
# pipeline is deterministic and offline-friendly after the first run.
# =============================================================================
import functools, json, time
import numpy as np, pandas as pd
from scipy.stats import pearsonr
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cnt_io
from cnt_io import cfg

@functools.lru_cache(maxsize=1)
def corum_gene_features():
    """Per-gene CORUM complex membership features."""
    complexes = cnt_io.load_corum()   # dict complex_id -> set(genes)
    from collections import defaultdict
    n_cx = defaultdict(int); max_sz = defaultdict(int)
    for cid, genes in complexes.items():
        sz = len(genes)
        for g in genes:
            n_cx[g] += 1
            if sz > max_sz[g]: max_sz[g] = sz
    genes = sorted(n_cx)
    return pd.DataFrame({"gene": genes,
                         "is_complex_subunit": [1]*len(genes),
                         "corum_n_complexes": [n_cx[g] for g in genes],
                         "corum_max_complex_size": [max_sz[g] for g in genes]})

@functools.lru_cache(maxsize=1)
def ensembl_coords():
    """Gene -> (chrom,start,end,strand). Read from the data_download stage-17
    cache (ensembl_gene_coords.parquet); the analysis pipeline never queries the
    network. Run `data_download` stage 17 first if this file is missing."""
    cache = cfg.PATHS["ensembl_coords"]
    if not cache.exists():
        raise FileNotFoundError(
            f"{cache} missing. Run data_download stage 17 (17_ensembl_coords.py) "
            "to fetch and cache Ensembl gene coordinates before running analysis.")
    cc = pd.read_parquet(cache)
    return cc[cc.chrom.isin([str(i) for i in range(1,23)]+["X","Y"])].copy()

@functools.lru_cache(maxsize=1)
def regulatory_gene_features():
    """Methylation + genomic regulatory features per gene (paper-exact)."""
    df = cnt_io.load_str_omic()
    d = df.dropna(subset=["meth","rna","cn_adjusted"]).copy()
    recs = []
    for gene, g in d.groupby("gene"):
        if len(g) < 40: continue
        zs, ws = [], []
        for tc, gt in g.groupby("tumor_code"):
            if len(gt) < 20 or gt.meth.std() == 0 or gt.rna.std() == 0: continue
            r = pearsonr(gt.meth, gt.rna)[0]
            if np.isfinite(r): zs.append(np.arctanh(np.clip(r,-0.999,0.999))); ws.append(len(gt)-3)
        meth_rna = np.tanh(np.average(zs, weights=ws)) if zs else np.nan
        recs.append({"gene":gene,"mean_promoter_meth":g.meth.mean(),
                     "std_promoter_meth":g.meth.std(),"meth_rna_corr":meth_rna})
    methfeat = pd.DataFrame(recs)
    cc = ensembl_coords().copy()
    cc["gene_length"] = (cc.end - cc.start).abs()
    cc["tss"] = np.where(cc.strand == 1, cc.start, cc.end)
    cc["log_gene_length"] = np.log10(cc.gene_length + 1)
    cc = cc.sort_values(["chrom","tss"]).reset_index(drop=True)
    dens = np.zeros(len(cc), dtype=int)
    for chrom, idx in cc.groupby("chrom").groups.items():
        tss = cc.loc[idx, "tss"].values
        for j, t in enumerate(tss):
            dens[idx[j]] = np.sum(np.abs(tss - t) <= 500_000) - 1
    cc["gene_density_1mb"] = dens
    return methfeat.merge(cc[["gene","gene_length","log_gene_length","gene_density_1mb","chrom"]],
                          on="gene", how="outer")
