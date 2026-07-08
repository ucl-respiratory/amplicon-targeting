# =============================================================================
# 10_transmission_reconciliation.py  --  Fig 7: reconciling attenuation with
# inheritance. (A) protein tracks CN ~half as strongly as mRNA; (B) between-gene,
# protein responsiveness is inherited from mRNA responsiveness (r=0.73);
# (C) transmission explains most between-gene variance, buffering adds little.
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.stats import pearsonr, gaussian_kde
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

import _mechanism
from cnt_io import DIR_TAB, savefig, record_values

def main():
    att = _mechanism.per_gene_attenuation()
    a = att.dropna(subset=["cn_rna_corr","cn_prot_corr"]).copy()

    # total_r2 from the Shapley decomposition (stage 12). Recompute if absent.
    dec_path = DIR_TAB / "mechanism_variance_decomposition.csv"
    if dec_path.exists():
        total_r2 = float(pd.read_csv(dec_path).total_r2.iloc[0])
    else:
        import importlib.util
        spec = importlib.util.spec_from_file_location("s12", Path(__file__).parent / "12_mechanism_decomposition.py")
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
        gl = m.build_gene_table(); _, total_r2 = m.shapley(gl)

    r = pearsonr(a.cn_rna_corr, a.cn_prot_corr)[0]

    plt.rcParams.update({"font.size":8,"axes.titlesize":8.5,"axes.labelsize":8,"legend.fontsize":7,
        "xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,
        "axes.titlelocation":"left"})
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3))

    ax = axes[0]
    ax.hist(a.cn_rna_corr, bins=50, alpha=0.6, color="#4a6fa5", density=True,
            label=f"CN\u2192mRNA (mean {a.cn_rna_corr.mean():.2f})")
    ax.hist(a.cn_prot_corr, bins=50, alpha=0.6, color="#5a9367", density=True,
            label=f"CN\u2192protein (mean {a.cn_prot_corr.mean():.2f})")
    ax.axvline(a.cn_rna_corr.mean(), color="#4a6fa5", lw=1.2, ls="--")
    ax.axvline(a.cn_prot_corr.mean(), color="#5a9367", lw=1.2, ls="--")
    ax.axvline(0, color="k", lw=0.6, ls=":")
    ax.set_xlabel("Per-gene correlation with copy number"); ax.set_ylabel("Density")
    ax.set_title("A  Protein tracks CN ~half as strongly\nas mRNA (attenuation)")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")

    ax = axes[1]
    x = a.cn_rna_corr.values; y = a.cn_prot_corr.values
    xy = np.vstack([x, y]); z = gaussian_kde(xy)(xy); idx = z.argsort()
    ax.scatter(x[idx], y[idx], c=z[idx], s=5, cmap="viridis", rasterized=True)
    lo, hi = -0.3, 0.8
    ax.plot([lo, hi], [lo, hi], color="grey", lw=0.8, ls=":", label="protein = mRNA (no buffering)")
    b1, b0 = np.polyfit(x, y, 1); xs = np.array([lo, hi])
    ax.plot(xs, b1*xs+b0, color="#c0504d", lw=1.5, label=f"fit: r={r:.2f}")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("CN\u2192mRNA correlation (transmission)")
    ax.set_ylabel("CN\u2192protein correlation (responsiveness)")
    ax.set_title("B  Which genes have responsive protein\nfollows which have responsive mRNA")
    ax.legend(frameon=False, fontsize=6.3, loc="upper left")

    ax = axes[2]
    r2_rna_alone = r**2
    vals = [r2_rna_alone, total_r2 - r2_rna_alone]
    ax.bar([0], [vals[0]], color="#4a6fa5", label="transmission (CN\u2192mRNA) alone")
    ax.bar([0], [vals[1]], bottom=[vals[0]], color="#c07a2d", label="all other features add")
    ax.set_xticks([0]); ax.set_xticklabels(["between-gene\nvariance in protein\nresponsiveness"], fontsize=7)
    ax.set_ylabel("R\u00b2 explained"); ax.set_ylim(0, 0.65)
    ax.text(0, r2_rna_alone/2, f"{r2_rna_alone:.2f}", ha="center", va="center", color="white", fontsize=8, weight="bold")
    ax.text(0, r2_rna_alone+vals[1]/2, f"+{vals[1]:.2f}", ha="center", va="center", color="white", fontsize=7)
    ax.set_title("C  Transmission explains most\nbetween-gene variance; buffering adds little")
    ax.legend(frameon=False, fontsize=6.3, loc="upper right")

    fig.suptitle("Reconciling the two facts: protein is attenuated relative to mRNA (A), "
                 "yet protein responsiveness is inherited from mRNA responsiveness (B,C)",
                 fontsize=8.5, y=1.02)
    fig.tight_layout(); savefig(fig, "transmission_reconciliation.png"); plt.close(fig)

    record_values("10_transmission_reconciliation", {
        "mean_cn_rna_corr": round(float(a.cn_rna_corr.mean()), 4),
        "mean_cn_prot_corr": round(float(a.cn_prot_corr.mean()), 4),
        "between_gene_r": round(float(r), 4),
        "between_gene_r2": round(float(r**2), 4),
        "mean_attenuation": round(float(a.attenuation.mean()), 4),
        "n_genes": int(len(a)),
    })
    print(f"mean_cn_rna={a.cn_rna_corr.mean():.4f} mean_cn_prot={a.cn_prot_corr.mean():.4f} r={r:.4f} n={len(a)}")

if __name__ == "__main__":
    main()
