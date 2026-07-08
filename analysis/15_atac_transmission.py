# =============================================================================
# 15_atac_transmission.py  --  Fig 12: tumour chromatin accessibility explains
# dosage transmission that CPTAC features cannot. Nested cross-validated R^2 of
# CN->mRNA transmission from CPTAC features vs +tumour promoter ATAC, across the
# 5 ATAC-matched cancer types. Headline: accessibility ~2.7x the explained R^2.
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.stats import pearsonr, gaussian_kde
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cnt_io, _features
from cnt_io import DIR_TAB, savefig, record_values, cfg

ATAC_TYPES = ["LUAD","LSCC","CCRCC","GBM","UCEC"]
CPTAC_FEATS = ['mean_promoter_meth','meth_rna_corr','log_gene_length','gene_density_1mb',
               'is_complex_subunit','corum_n_complexes','corum_max_complex_size']

def _cv_r2(X, y):
    kf = KFold(5, shuffle=True, random_state=cfg.SEED_BOOTSTRAP)
    ps = np.zeros(len(y))
    for tr, te in kf.split(X):
        ps[te] = LinearRegression().fit(X[tr], y[tr]).predict(X[te])
    return r2_score(y, ps)

def _transmission(df, ct):
    g = df[(df.tumor_code == ct)].dropna(subset=["cn_adjusted","rna"])
    trans = {}
    for gene, gg in g.groupby("gene"):
        if len(gg) < 20 or gg.cn_adjusted.std() == 0 or gg.rna.std() == 0: continue
        trans[gene] = pearsonr(gg.cn_adjusted, gg.rna)[0]
    return pd.DataFrame({"gene": list(trans), "transmission": list(trans.values())})

def _build(df, reg, corum, accdf, ct):
    t = _transmission(df, ct)
    t = t.merge(reg[["gene","mean_promoter_meth","meth_rna_corr","log_gene_length","gene_density_1mb"]],
                on="gene", how="left")
    t = t.merge(corum, on="gene", how="left")
    for c in ["is_complex_subunit","corum_n_complexes","corum_max_complex_size"]:
        t[c] = t[c].fillna(0)
    t = t.merge(accdf[ct].rename("atac").reset_index().rename(columns={"index":"gene"}),
                on="gene", how="left")
    return t

def main():
    df = cnt_io.load_str_omic()
    reg = _features.regulatory_gene_features()
    corum = _features.corum_gene_features()
    accdf = cnt_io.load_atac()
    if accdf.index.name is None: accdf.index.name = "gene"

    rows = []
    for ct in ATAC_TYPES:
        t = _build(df, reg, corum, accdf, ct).dropna(subset=["transmission"]+CPTAC_FEATS+["atac"])
        y = t.transmission.values
        r2_cptac = _cv_r2(t[CPTAC_FEATS].values, y)
        r2_both  = _cv_r2(t[CPTAC_FEATS+["atac"]].values, y)
        r2_atac  = _cv_r2(t[["atac"]].values, y)
        rr = pearsonr(t.transmission, t.atac)[0]
        rows.append({"tumor_code":ct,"n_genes":len(t),"raw_corr":rr,"R2_CPTAC_feats":r2_cptac,
                     "R2_plus_ATAC":r2_both,"R2_ATAC_alone":r2_atac,"delta":r2_both-r2_cptac})
    comp = pd.DataFrame(rows)
    comp.to_csv(DIR_TAB / "atac_transmission_model.csv", index=False)

    plt.rcParams.update({"font.size":8,"axes.titlesize":8.5,"axes.labelsize":8,"legend.fontsize":7,
        "xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,
        "axes.titlelocation":"left"})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    ax = axes[0]
    t = _build(df, reg, corum, accdf, "LSCC").dropna(subset=["transmission","atac"])
    x = t.atac.values; y = t.transmission.values
    xy = np.vstack([x,y]); z = gaussian_kde(xy)(xy); idx = z.argsort()
    ax.scatter(x[idx],y[idx],c=z[idx],s=5,cmap="viridis",rasterized=True)
    b1,b0 = np.polyfit(x,y,1); xs = np.array([x.min(),x.max()]); ax.plot(xs,b1*xs+b0,color="#c0504d",lw=1.5)
    ax.set_xlabel("Promoter chromatin accessibility (LUSC tumour ATAC, log2)")
    ax.set_ylabel("CN\u2192mRNA transmission (per gene)")
    ax.set_title(f"A  Accessible promoters transmit dosage better\nLSCC: r={pearsonr(x,y)[0]:.2f}, n={len(t)} genes")
    ax = axes[1]
    w = 0.38; yb = np.arange(len(comp))
    ax.barh(yb+w/2, comp.R2_CPTAC_feats, w, color="#8aa1b4", label="CPTAC features")
    ax.barh(yb-w/2, comp.R2_plus_ATAC, w, color="#c07a2d", label="CPTAC + tumour ATAC")
    ax.set_yticks(yb); ax.set_yticklabels(comp.tumor_code)
    for i, r in enumerate(comp.itertuples()):
        ax.text(r.R2_CPTAC_feats+0.002, i+w/2, f"{r.R2_CPTAC_feats:.3f}", va="center", fontsize=6.3)
        ax.text(r.R2_plus_ATAC+0.002, i-w/2, f"{r.R2_plus_ATAC:.3f}", va="center", fontsize=6.3)
    ax.set_xlabel("Cross-validated R\u00b2 (predicting CN\u2192mRNA transmission)")
    ax.set_title("B  Accessibility ~2.7\u00d7 the explained transmission")
    ax.legend(frameon=False, fontsize=6.8, loc="lower right")
    ax.set_xlim(0, 0.14); ax.invert_yaxis()
    fig.suptitle("Tumour chromatin accessibility explains dosage transmission that CPTAC features "
                 "cannot (5 matched cancer types)", fontsize=8.6, y=1.02)
    fig.tight_layout(); savefig(fig, "atac_transmission.png"); plt.close(fig)

    record_values("15_atac_transmission", {
        "mean_R2_CPTAC": round(float(comp.R2_CPTAC_feats.mean()),4),
        "mean_R2_plus_ATAC": round(float(comp.R2_plus_ATAC.mean()),4),
        "fold_gain": round(float(comp.R2_plus_ATAC.mean()/comp.R2_CPTAC_feats.mean()),3),
        "per_type": {r.tumor_code: {"R2_CPTAC":round(r.R2_CPTAC_feats,4),
                     "R2_plus_ATAC":round(r.R2_plus_ATAC,4),"raw_corr":round(r.raw_corr,4),
                     "n_genes":int(r.n_genes)} for r in comp.itertuples()},
    })
    print(comp[["tumor_code","n_genes","R2_CPTAC_feats","R2_plus_ATAC","delta"]].to_string(index=False))
    print(f"mean CPTAC={comp.R2_CPTAC_feats.mean():.4f} +ATAC={comp.R2_plus_ATAC.mean():.4f} "
          f"fold={comp.R2_plus_ATAC.mean()/comp.R2_CPTAC_feats.mean():.2f}")

if __name__ == "__main__":
    main()
