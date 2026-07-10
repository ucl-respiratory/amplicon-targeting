# =============================================================================
# integrated/_io.py  --  shared loaders + transmission/co-elevation computations
# -----------------------------------------------------------------------------
# Self-contained (does not import analysis/ or gene_intrinsic/). Reads the
# data_download/from_source tables via config.PATHS and reproduces the two
# source pipelines' load-bearing computations in one integrated module:
#   - per-gene dosage transmission / responsiveness / attenuation
#   - recurrent-amplicon co-elevation (one-sided Fisher, BH-FDR)
# Parameters all come from config.py (CN>=1.4 basis).
# =============================================================================
import functools, importlib.util, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import fisher_exact

_cfg_path = Path(__file__).resolve().parent / "config.py"
_spec = importlib.util.spec_from_file_location("integrated_config", _cfg_path)
cfg = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(cfg)

# ---- loaders ----------------------------------------------------------------
@functools.lru_cache(maxsize=1)
def load_str_omic():
    """Wide (gene x caseid) feature matrix from data_download stage 13."""
    p = cfg.PATHS["str_omic"]
    if not p.exists():
        raise FileNotFoundError(
            f"{p} missing — run: cd data_download && Rscript ../integrated/data_prep.R")
    return pd.read_csv(p, low_memory=False)

@functools.lru_cache(maxsize=1)
def load_atac():
    return pd.read_parquet(cfg.PATHS["atac"])

@functools.lru_cache(maxsize=1)
def load_topology():
    return pd.read_parquet(cfg.PATHS["topology"])

@functools.lru_cache(maxsize=1)
def load_dep_prob():
    df = pd.read_csv(cfg.PATHS["dep_prob"], index_col=0)
    df.columns = [c.split(" (")[0] for c in df.columns]
    return df

# ---- GI committed exports (predictor outputs; reused until fresh re-fit) -----
@functools.lru_cache(maxsize=1)
def load_atlas():
    """Transmissibility atlas: measured observed + predicted OOF prior, per gene."""
    return pd.read_csv(cfg.GI_PATHS["atlas"])

@functools.lru_cache(maxsize=1)
def load_feature_table():
    return pd.read_csv(cfg.GI_PATHS["feature_table"])

@functools.lru_cache(maxsize=1)
def load_transfer():
    return pd.read_csv(cfg.GI_PATHS["transfer"])

@functools.lru_cache(maxsize=1)
def load_leads_validation():
    return pd.read_csv(cfg.GI_PATHS["leads_valid"])

# ---- mechanism: transmission / responsiveness / attenuation -----------------
_MIN_N_PER_TYPE = 20
_MIN_TYPES = 2

def _pearson_by_gene(sub, xcol, ycol):
    g = sub.groupby("gene"); n = g.size()
    sx = g[xcol].sum(); sy = g[ycol].sum()
    sxx = g[xcol].apply(lambda t: (t**2).sum())
    syy = g[ycol].apply(lambda t: (t**2).sum())
    sxy = g.apply(lambda t: (t[xcol]*t[ycol]).sum(), include_groups=False)
    cov = sxy - sx*sy/n
    vx = sxx - sx**2/n; vy = syy - sy**2/n
    return cov/np.sqrt(vx*vy), n

def _fisher_avg(r, n):
    r = np.clip(np.asarray(r, float), -0.999, 0.999)
    z = np.arctanh(r); w = np.asarray(n, float) - 3
    return np.tanh((z*w).sum()/w.sum())

@functools.lru_cache(maxsize=1)
def per_type_attenuation():
    """Per (gene, tumor_code): r_cn_rna, r_cn_prot, n. RAW cn, joint dropna."""
    d = load_str_omic().dropna(subset=["cn", "rna", "prot"]).copy()
    rows = []
    for tc, sub in d.groupby("tumor_code"):
        s = sub[["gene", "cn", "rna", "prot"]]
        r_rna, n  = _pearson_by_gene(s, "cn", "rna")
        r_prot, _ = _pearson_by_gene(s, "cn", "prot")
        rows.append(pd.DataFrame({"gene": r_rna.index, "r_cn_rna": r_rna.values,
                                  "r_cn_prot": r_prot.values, "n": n.values,
                                  "tumor_code": tc}))
    return pd.concat(rows, ignore_index=True)

@functools.lru_cache(maxsize=1)
def per_gene_attenuation():
    """Pan-cancer per-gene cn_rna_corr / cn_prot_corr / attenuation."""
    pt = per_type_attenuation().dropna(subset=["r_cn_rna", "r_cn_prot"])
    pt = pt[pt.n >= _MIN_N_PER_TYPE]
    agg = []
    for gene, sub in pt.groupby("gene"):
        if len(sub) < _MIN_TYPES: continue
        r_rna  = _fisher_avg(sub.r_cn_rna.values,  sub.n.values)
        r_prot = _fisher_avg(sub.r_cn_prot.values, sub.n.values)
        agg.append((gene, r_rna, r_prot, r_rna - r_prot, len(sub), int(sub.n.sum())))
    return pd.DataFrame(agg, columns=["gene","cn_rna_corr","cn_prot_corr",
                                      "attenuation","n_types","n_samples_total"])

# ---- recurrent-amplicon co-elevation ----------------------------------------
def _bh_fdr(p):
    p = np.asarray(p); n = len(p); order = np.argsort(p)
    ranked = np.empty(n); cummin = 1.0
    for i in range(n-1, -1, -1):
        idx = order[i]; val = p[idx]*n/(i+1); cummin = min(cummin, val); ranked[idx] = cummin
    return np.clip(ranked, 0, 1)

@functools.lru_cache(maxsize=1)
def amplicon_coelevation():
    """Per (tumor_code, band, gene) one-sided Fisher: high tissue-relative protein
    enriched in band-amplified vs non-amplified tumours. CN>=1.4 amplified;
    recurrent band = amplified in >=15% and >=8 samples. BH-FDR across tests."""
    d = load_str_omic().dropna(
        subset=["cn_adjusted","cytogenetic_location","prot.rel.tissue"]).copy()
    d["amplified"] = (d.cn_adjusted >= cfg.AMP_THRESHOLD).astype(int)
    d["prot_high"] = (d["prot.rel.tissue"] > cfg.REL_TISSUE_HI).astype(int)

    band_samp = d.groupby(["tumor_code","cytogenetic_location","caseid"]).agg(
        frac_genes_amp=("amplified","mean")).reset_index()
    band_samp["band_amplified"] = (band_samp.frac_genes_amp >= cfg.COELEV_MIN_FRAC).astype(int)

    br = band_samp.groupby(["tumor_code","cytogenetic_location"]).agg(
        n_amp=("band_amplified","sum"), n=("band_amplified","size")).reset_index()
    br["recurrence"] = br.n_amp / br.n
    recurrent = br[(br.recurrence >= 0.15) & (br.n_amp >= cfg.RECURRENCE_MIN_N)]
    rec_keys = set(zip(recurrent.tumor_code, recurrent.cytogenetic_location))

    d2 = d.merge(band_samp[["tumor_code","cytogenetic_location","caseid","band_amplified"]],
                 on=["tumor_code","cytogenetic_location","caseid"], how="left")
    d2["is_rec"] = [(t,c) in rec_keys for t,c in
                    zip(d2.tumor_code, d2.cytogenetic_location)]
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
                     "is_surface":int(g.is_surface.iloc[0]) if "is_surface" in g else 0,
                     "is_secreted":int(g.is_secreted.iloc[0]) if "is_secreted" in g else 0})
    co = pd.DataFrame(recs).dropna(subset=["fisher_p"]).reset_index(drop=True)
    co["fdr"] = _bh_fdr(co.fisher_p.values)
    return co

def record(stage, d):
    """Write a per-stage values JSON into reports/values/ for the manifest."""
    import json
    vd = cfg.DIR_REP / "values"; vd.mkdir(parents=True, exist_ok=True)
    (vd / f"{stage}.json").write_text(json.dumps(d, indent=2, default=float))

# ---- self-contained publication style (so modules run standalone) -----------
def apply_style(sizes=(9, 8, 7)):
    """Role-mapped font ladder + clean spines. Mirrors the figure-style skill so
    modules produce publication-grade figures when run as scripts."""
    import matplotlib as mpl
    base, mid, small = sizes
    mpl.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 200,
        "font.size": base, "axes.titlesize": base, "axes.labelsize": base,
        "xtick.labelsize": small, "ytick.labelsize": small,
        "legend.fontsize": mid, "axes.spines.top": False, "axes.spines.right": False,
        "axes.titleweight": "regular", "axes.grid": False,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "svg.fonttype": "none", "pdf.fonttype": 42,
    })

def panel_letter(ax, letter, case="lower", dx=-0.08, dy=1.06):
    ax.text(dx, dy, letter if case == "lower" else letter.upper(),
            transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", ha="right")
