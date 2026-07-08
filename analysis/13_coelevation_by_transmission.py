# =============================================================================
# 13_coelevation_by_transmission.py  --  Fig 10: on real amplicons, whether a
# co-amplified gene's protein co-elevates is driven by CN->mRNA transmission
# (r~+0.42), not post-transcriptional buffering (r~+0.07).
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.stats import pearsonr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import _amplicon, _mechanism
from cnt_io import DIR_TAB, savefig, record_values

def main():
    co = _amplicon.amplicon_coelevation()
    att = _mechanism.per_gene_attenuation()
    m = co.merge(att[["gene","cn_rna_corr","cn_prot_corr","attenuation"]], on="gene",
                 how="left").dropna(subset=["cn_rna_corr"])
    r_trans = pearsonr(m.cn_rna_corr, m.p_high_amp)[0]
    r_atten = pearsonr(m.attenuation, m.p_high_amp)[0]
    m["trans_tercile"] = pd.qcut(m.cn_rna_corr, 3, labels=["low","mid","high"])
    tab = m.groupby("trans_tercile", observed=True).agg(
        mean_coelev_P=("p_high_amp","mean"),
        frac_reliable=("fdr", lambda x:(x<0.1).mean())).reset_index()

    plt.rcParams.update({"font.size":8,"axes.titlesize":8.5,"axes.labelsize":8,"legend.fontsize":7,
        "xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,
        "axes.titlelocation":"left"})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    ax.bar(range(3), tab.mean_coelev_P, color=["#c7ccd1","#7a97b5","#4a6fa5"])
    ax.set_xticks(range(3)); ax.set_xticklabels(["low","mid","high"])
    ax.set_xlabel("CN\u2192mRNA transmission tercile"); ax.set_ylabel("P(protein high | amplicon present)")
    for i,v in enumerate(tab.mean_coelev_P): ax.text(i, v+0.01, f"{v:.2f}", ha="center", fontsize=7)
    ax.set_title("A  Protein co-elevation rises with transmission"); ax.set_ylim(0, 0.75)
    ax = axes[1]
    ax.bar(range(3), tab.frac_reliable, color=["#c7ccd1","#c9a97a","#c07a2d"])
    ax.set_xticks(range(3)); ax.set_xticklabels(["low","mid","high"])
    ax.set_xlabel("CN\u2192mRNA transmission tercile"); ax.set_ylabel("Fraction reliably co-elevated (FDR<0.1)")
    for i,v in enumerate(tab.frac_reliable): ax.text(i, v+0.005, f"{v:.0%}", ha="center", fontsize=7)
    ax.set_title("B  Reliable co-targets concentrate in high-transmission genes"); ax.set_ylim(0, 0.45)
    fig.suptitle(f"On amplicons: protein co-elevation is driven by CN\u2192mRNA transmission "
                 f"(r={r_trans:+.2f}), not buffering (r={r_atten:+.2f})", fontsize=8, y=1.03)
    fig.tight_layout(); savefig(fig, "coelevation_by_transmission.png"); plt.close(fig)

    record_values("13_coelevation_by_transmission", {
        "r_transmission_vs_coelev": round(float(r_trans),4),
        "r_attenuation_vs_coelev": round(float(r_atten),4),
        "n_amplicon_gene_events": int(len(m)),
        "coelev_P_by_tercile": {r.trans_tercile: round(r.mean_coelev_P,4) for r in tab.itertuples()},
        "frac_reliable_by_tercile": {r.trans_tercile: round(r.frac_reliable,4) for r in tab.itertuples()},
    })
    print(f"r_trans={r_trans:+.4f} r_atten={r_atten:+.4f} n={len(m)}")

if __name__ == "__main__":
    main()
