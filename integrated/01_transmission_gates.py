#!/usr/bin/env python3
# =============================================================================
# 01_transmission_gates.py  --  The transmission cascade and where it is gated.
# -----------------------------------------------------------------------------
# Quantifies each step of amplicon -> mRNA -> protein -> surface from the real
# data_download tables, and produces Figure 2: the "where is the signal"
# diagnostic. Every arrow is a measured transmitted fraction; every gate is a
# measured determinant.
#
#   step 1  CN -> mRNA      per-gene corr(cn, rna), pan-cancer Fisher-z mean
#   step 2  mRNA -> protein rank preserved: corr of CN->mRNA vs CN->protein
#           transmission across genes; attenuation = transmission - responsiveness
#   gate    chromatin       nested-CV R^2 of transmission from gene features
#           with vs without tumour promoter ATAC accessibility (5 ATAC types)
#
# Methylation features are part of the gene-feature block; in this build they
# resolve to NA (ChAMP unavailable), so the CPTAC-feature R^2 is reported as
# "gene features (ATAC excluded)" and the ATAC delta is the load-bearing result.
#
# Outputs: figures/fig2_transmission_gates.png, tables/transmission_gates.csv
# =============================================================================
import sys, warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cnt_shared as ish
from cnt_shared import cfg
ish.apply_style(sizes=(9, 8, 7))

ATAC_TYPES = ["LUAD", "LSCC", "CCRCC", "GBM", "UCEC"]

# ---- step 1+2: transmission / responsiveness / attenuation ------------------
def cascade_summary():
    pg = ish.per_gene_attenuation()
    return {
        "n_genes": int(len(pg)),
        "cn_rna_median":  float(pg.cn_rna_corr.median()),
        "cn_prot_median": float(pg.cn_prot_corr.median()),
        "attenuation_median": float(pg.attenuation.median()),
        "frac_attenuated": float((pg.attenuation > 0).mean()),
        # rank preservation: does the CN->mRNA ranking survive to CN->protein?
        "rank_preserved_rho": float(spearmanr(pg.cn_rna_corr, pg.cn_prot_corr).correlation),
    }, pg

# ---- chromatin gate: nested-CV R^2 with vs without ATAC ---------------------
def _cv_r2(X, y, seed):
    kf = KFold(5, shuffle=True, random_state=seed)
    ps = np.zeros(len(y))
    for tr, te in kf.split(X):
        ps[te] = LinearRegression().fit(X[tr], y[tr]).predict(X[te])
    return r2_score(y, ps)

def _transmission(df, ct):
    g = df[df.tumor_code == ct].dropna(subset=["cn_adjusted", "rna"])
    tr = {}
    for gene, gg in g.groupby("gene"):
        if len(gg) < 20 or gg.cn_adjusted.std() == 0 or gg.rna.std() == 0: continue
        tr[gene] = pearsonr(gg.cn_adjusted, gg.rna)[0]
    return pd.DataFrame({"gene": list(tr), "transmission": list(tr.values())})

# gene-property features available WITHOUT methylation (local gene density is
# positional/structural, derived from genomic coordinates regardless of ChAMP).
GENE_FEATS = ["gene_density_1mb"]

def _gene_density():
    """Genes within +/-500 kb of each gene's locus, from the genemap coordinates
    (structural baseline for the chromatin gate; independent of any omic layer)."""
    import pyreadr
    gm = pyreadr.read_r(str(cfg.PARSED / "genemap.RData"))["genemap"]
    gm = gm[["symbol", "Chromosome", "POS"]].dropna().copy()
    gm["POS"] = gm["POS"].astype(float)
    out = {}
    for chrom, cc in gm.groupby("Chromosome"):
        pos = cc["POS"].values
        for sym, p in zip(cc["symbol"], pos):
            out[sym] = int(((pos >= p - 5e5) & (pos <= p + 5e5)).sum()) - 1
    return pd.DataFrame({"gene": list(out), "gene_density_1mb": list(out.values())})

def chromatin_gate():
    df = ish.load_str_omic()
    # positional features from the feature table if present, else derive here
    reg_cols = [c for c in ["mean_promoter_meth", "meth_rna_corr"] if c in df.columns]
    reg = (df[["gene"] + reg_cols].drop_duplicates("gene")
           if reg_cols else pd.DataFrame({"gene": df.gene.unique()}))
    # methylation features assembled from GDC EPIC betas via the manifest probe->
    # gene map (build_meth_features.py, no ChAMP). Merge if the standalone table
    # exists and the omic layer didn't already carry them.
    mfp = cfg.DIR_TAB / "meth_gene_features.csv"
    if mfp.exists() and "mean_promoter_meth" not in reg.columns:
        mf = pd.read_csv(mfp)[["gene", "mean_promoter_meth", "meth_rna_corr"]]
        reg = reg.merge(mf, on="gene", how="left")
    # structural gene-density baseline (genome coordinates), always available
    try:
        reg = reg.merge(_gene_density(), on="gene", how="left")
    except Exception as e:
        print("gene-density derivation skipped:", str(e)[:80])
    # structural baseline stays clean (gene density only); methylation is its own
    # gate. mean_promoter_meth has full gene coverage; meth_rna_corr is sparser, so
    # it is used only where present and never forced into the row-wise dropna.
    struct_feats = [f for f in GENE_FEATS if f in reg.columns]
    have_meth = "mean_promoter_meth" in reg.columns and reg["mean_promoter_meth"].notna().any()
    meth_feats = [f for f in ["mean_promoter_meth"] if f in reg.columns] if have_meth else []

    accdf = ish.load_atac()
    # ATAC parquet carries an explicit `gene` column; use it as the index so the
    # per-type accessibility Series carries gene labels (not the RangeIndex).
    if "gene" in accdf.columns:
        accdf = accdf.set_index("gene")
    elif accdf.index.name is None:
        accdf.index.name = "gene"
    rows = []
    for ct in ATAC_TYPES:
        if ct not in accdf.columns: continue
        t = _transmission(df, ct).merge(reg, on="gene", how="left")
        acc_ct = accdf[ct].rename("atac").reset_index()
        acc_ct.columns = ["gene", "atac"]
        t = t.merge(acc_ct, on="gene", how="left")
        t = t.dropna(subset=["transmission", "atac"] + struct_feats + meth_feats)
        if len(t) < 50: continue
        y = t.transmission.values
        r2_struct = _cv_r2(t[struct_feats].values, y, cfg.SEED) if struct_feats else 0.0
        r2_meth  = _cv_r2(t[struct_feats + meth_feats].values, y, cfg.SEED) if meth_feats else r2_struct
        r2_atac_only = _cv_r2(t[["atac"]].values, y, cfg.SEED)
        r2_struct_atac = _cv_r2(t[struct_feats + ["atac"]].values, y, cfg.SEED)
        r2_all = _cv_r2(t[struct_feats + meth_feats + ["atac"]].values, y, cfg.SEED)
        rows.append({"tumor_code": ct, "n_genes": len(t),
                     "R2_gene_feats": r2_struct,
                     "R2_plus_meth": r2_meth,
                     "R2_plus_ATAC": r2_struct_atac,
                     "R2_ATAC_alone": r2_atac_only,
                     "R2_all_gates": r2_all,
                     "meth_delta": r2_meth - r2_struct,
                     "atac_delta": r2_struct_atac - r2_struct,
                     "raw_corr": pearsonr(t.transmission, t.atac)[0]})
    return pd.DataFrame(rows), have_meth

def main():
    summ, pg = cascade_summary()
    gate, have_meth = chromatin_gate()
    gate.to_csv(cfg.DIR_TAB / "transmission_gates.csv", index=False)

    # ---- Figure 2: three-panel cascade diagnostic --------------------------
    fig = plt.figure(figsize=(13, 4.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1], wspace=0.34)
    C_RNA, C_PROT, C_ATAC = "#2c6fbb", "#c0392b", "#e67e22"

    # panel a: per-gene CN->mRNA vs CN->protein transmission (rank preserved)
    axa = fig.add_subplot(gs[0, 0])
    axa.scatter(pg.cn_rna_corr, pg.cn_prot_corr, s=5, alpha=0.25,
                color="#555", rasterized=True, edgecolors="none")
    lim = [-0.4, 1.0]
    axa.plot(lim, lim, "--", color="#999", lw=1, zorder=1)
    axa.set_xlim(lim); axa.set_ylim(lim)
    axa.set_xlabel("CN → mRNA transmission  (per-gene r)")
    axa.set_ylabel("CN → protein responsiveness  (per-gene r)")
    axa.set_title(f"Rank preserved to protein (ρ = {summ['rank_preserved_rho']:.2f}),\n"
                  f"attenuated below the diagonal in {summ['frac_attenuated']*100:.0f}% of genes",
                  fontsize=9, loc="left")
    axa.annotate(f"median transmission {summ['cn_rna_median']:.2f}",
                 xy=(summ["cn_rna_median"], -0.34), color=C_RNA, fontsize=7.5, ha="center")
    axa.annotate(f"median responsiveness {summ['cn_prot_median']:.2f}",
                 xy=(-0.36, summ["cn_prot_median"]), color=C_PROT, fontsize=7.5,
                 rotation=90, va="center")

    # panel b: attenuation distribution
    axb = fig.add_subplot(gs[0, 1])
    axb.hist(pg.attenuation, bins=50, color="#888", alpha=0.85)
    axb.axvline(0, color="#333", lw=1)
    axb.axvline(summ["attenuation_median"], color=C_PROT, lw=1.6)
    axb.set_xlabel("Attenuation  (transmission − responsiveness)")
    axb.set_ylabel(f"genes  (n = {summ['n_genes']:,})")
    axb.set_title(f"Post-transcriptional loss:\nmedian attenuation {summ['attenuation_median']:.2f}",
                  fontsize=9, loc="left")

    # panel c: regulatory gates — R^2 from structure, +methylation, +ATAC, all,
    # per cancer type (methylation and accessibility are complementary gates)
    axc = fig.add_subplot(gs[0, 2])
    C_METH = "#159090"
    if len(gate):
        x = np.arange(len(gate)); w = 0.20
        axc.bar(x - 1.5*w, gate.R2_gene_feats, w, color="#aab", label="gene structure")
        axc.bar(x - 0.5*w, gate.R2_plus_meth, w, color=C_METH, label="+ promoter methylation")
        axc.bar(x + 0.5*w, gate.R2_plus_ATAC, w, color=C_ATAC, label="+ tumour ATAC")
        axc.bar(x + 1.5*w, gate.R2_all_gates, w, color="#6c5b9c", label="all gates")
        axc.set_xticks(x); axc.set_xticklabels(gate.tumor_code, rotation=0, fontsize=7.5)
        m_mult = gate.R2_plus_meth.mean() / max(gate.R2_gene_feats.mean(), 1e-6)
        a_mult = gate.R2_plus_ATAC.mean() / max(gate.R2_gene_feats.mean(), 1e-6)
        axc.set_ylabel("cross-validated R²  of transmission")
        axc.set_title(f"Transmission is gated at chromatin:\nmethylation ×{m_mult:.0f}, "
                      f"accessibility ×{a_mult:.0f} over structure",
                      fontsize=9, loc="left")
        axc.legend(frameon=False, fontsize=6.5, loc="upper right")
    for ax, L in zip([axa, axb, axc], "abc"):
        ish.panel_letter(ax, L)

    fig.savefig(cfg.DIR_FIG / "fig2_transmission_gates.png", dpi=200, bbox_inches="tight")

    ish.record("01_transmission_gates", {
        **summ,
        "atac_types": list(gate.tumor_code) if len(gate) else [],
        "R2_gene_feats_mean": float(gate.R2_gene_feats.mean()) if len(gate) else None,
        "R2_plus_meth_mean": float(gate.R2_plus_meth.mean()) if len(gate) else None,
        "R2_plus_ATAC_mean": float(gate.R2_plus_ATAC.mean()) if len(gate) else None,
        "R2_all_gates_mean": float(gate.R2_all_gates.mean()) if len(gate) else None,
        "meth_delta_mean": float(gate.meth_delta.mean()) if len(gate) else None,
        "atac_delta_mean": float(gate.atac_delta.mean()) if len(gate) else None,
        "meth_n_genes": int(gate.n_genes.mean()) if len(gate) else None,
        "methylation_available": bool(have_meth),
    })
    print("cascade:", summ)
    print(gate.to_string(index=False) if len(gate) else "ATAC gate: no data yet")
    print("methylation features available:", have_meth)

if __name__ == "__main__":
    main()
