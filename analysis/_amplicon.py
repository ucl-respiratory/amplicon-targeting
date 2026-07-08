# =============================================================================
# _amplicon.py  --  recurrent-amplicon co-elevation nomination (shared)
# -----------------------------------------------------------------------------
# Per (tumor_type, cytoband, gene) one-sided Fisher test: is high tissue-relative
# protein (prot.rel.tissue > 0.8) enriched in band-amplified vs non-amplified
# tumours? Amplified band = >=50% of its genes at cn_adjusted>=1.4; recurrent =
# band amplified in >=15% AND >=8 samples of a type. BH-FDR across all tests.
# Used by stage 13 (mechanism) and stage 21 (nomination). Paper-exact.
# =============================================================================
import functools
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cnt_io
from cnt_io import cfg

def _bh_fdr(p):
    p = np.asarray(p); n = len(p); order = np.argsort(p)
    ranked = np.empty(n); cummin = 1.0
    for i in range(n-1, -1, -1):
        idx = order[i]; val = p[idx]*n/(i+1); cummin = min(cummin, val); ranked[idx] = cummin
    return np.clip(ranked, 0, 1)

@functools.lru_cache(maxsize=1)
def amplicon_coelevation():
    df = cnt_io.load_str_omic()
    d = df.dropna(subset=["cn_adjusted","cytogenetic_location","prot.rel.tissue"]).copy()
    d["amplified"] = (d.cn_adjusted >= cfg.AMP_THRESHOLD).astype(int)
    d["prot_high"] = (d["prot.rel.tissue"] > cfg.REL_TISSUE_HI).astype(int)

    band_samp = d.groupby(["tumor_code","cytogenetic_location","caseid"]).agg(
        frac_genes_amp=("amplified","mean"), n_genes=("amplified","size")).reset_index()
    band_samp["band_amplified"] = (band_samp.frac_genes_amp >= cfg.COELEV_MIN_FRAC).astype(int)

    band_recur = band_samp.groupby(["tumor_code","cytogenetic_location"]).agg(
        n_amp_samples=("band_amplified","sum"), n_samples=("band_amplified","size")).reset_index()
    band_recur["recurrence"] = band_recur.n_amp_samples / band_recur.n_samples
    recurrent = band_recur[(band_recur.recurrence >= 0.15) & (band_recur.n_amp_samples >= 8)]

    d2 = d.merge(band_samp[["tumor_code","cytogenetic_location","caseid","band_amplified"]],
                 on=["tumor_code","cytogenetic_location","caseid"], how="left")
    rec_keys = set(zip(recurrent.tumor_code, recurrent.cytogenetic_location))
    d2["is_rec"] = [(t,c) in rec_keys for t,c in zip(d2.tumor_code, d2.cytogenetic_location)]
    dd = d2[d2.is_rec].copy()

    recs = []
    for (tc, band, gene), g in dd.groupby(["tumor_code","cytogenetic_location","gene"]):
        amp = g[g.band_amplified == 1]; noamp = g[g.band_amplified == 0]
        if len(amp) < 8 or len(noamp) < 5: continue
        a = amp.prot_high.sum(); b = len(amp)-a
        c = noamp.prot_high.sum(); e = len(noamp)-c
        try: _, pval = fisher_exact([[a,b],[c,e]], alternative="greater")
        except Exception: pval = np.nan
        recs.append({"tumor_code":tc,"band":band,"gene":gene,
                     "p_high_amp":amp.prot_high.mean(),"p_high_noamp":noamp.prot_high.mean(),
                     "lift":amp.prot_high.mean()-noamp.prot_high.mean(),
                     "n_amp":len(amp),"n_noamp":len(noamp),"fisher_p":pval,
                     "is_surface":int(g.is_surface.iloc[0]),"is_secreted":int(g.is_secreted.iloc[0])})
    co = pd.DataFrame(recs).dropna(subset=["fisher_p"]).reset_index(drop=True)
    co["fdr"] = _bh_fdr(co.fisher_p.values)
    return co

@functools.lru_cache(maxsize=1)
def recurrent_amplicons():
    """(tumor_type, cytoband) recurrence table — for the amplicon landscape (Fig 2)."""
    df = cnt_io.load_str_omic()
    d = df.dropna(subset=["cn_adjusted","cytogenetic_location"]).copy()
    d["amplified"] = (d.cn_adjusted >= cfg.AMP_THRESHOLD).astype(int)
    band_samp = d.groupby(["tumor_code","cytogenetic_location","caseid"]).agg(
        frac_genes_amp=("amplified","mean")).reset_index()
    band_samp["band_amplified"] = (band_samp.frac_genes_amp >= cfg.COELEV_MIN_FRAC).astype(int)
    br = band_samp.groupby(["tumor_code","cytogenetic_location"]).agg(
        n_amp_samples=("band_amplified","sum"), n_samples=("band_amplified","size")).reset_index()
    br["recurrence"] = br.n_amp_samples / br.n_samples
    br["recurrent"] = (br.recurrence >= 0.15) & (br.n_amp_samples >= 8)
    return br
