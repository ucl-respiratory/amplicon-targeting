#!/usr/bin/env python3
# =============================================================================
# 00b_target_funnel.py -- Figure 1: target-nomination funnel. Every count is
# COMPUTED from the transmissibility atlas (all six CPTAC cohorts), the
# gene-property feature table (surfaceome proxy: TM domain / signal peptide,
# DepMap annotation), and the ectodomain-gated nominated antigen table. No
# hardcoded counts. Universe spans CCRCC, GBM, LSCC, LUAD, PDA, UCEC.
# Output: figures/fig2_target_funnel.png, reports/values/00b_target_funnel.json
# =============================================================================
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg
import cnt_shared as ish

TRANSMIT_MIN = 0.40

def main():
    # The funnel now traces the ACTUAL pan-cancer nomination path:
    #   genome universe -> transmitted -> co-elevated on a recurrent amplicon
    #   (Fisher FDR<0.1, from 00c) -> transmitted co-elevated -> surface-accessible
    #   (UniProt ectodomain gate) -> nominated antigens. Every count is computed.
    atlas = pd.read_csv(cfg.GI_PATHS["atlas"])
    n1 = len(atlas)
    n2 = int((atlas["observed_transmissibility"] >= TRANSMIT_MIN).sum())

    co = pd.read_csv(cfg.DIR_TAB / "amplicon_coelevation.csv")
    co_sig = co[co.fdr < cfg.FDR_ALPHA]
    n3 = int(co_sig.gene.nunique())                         # co-elevated on a recurrent amplicon
    trans_genes = set(atlas.loc[atlas.observed_transmissibility >= TRANSMIT_MIN, "gene"])
    n4 = int(len(set(co_sig.gene) & trans_genes))           # co-elevated AND transmitted

    nom = pd.read_csv(cfg.DIR_TAB / "adc_target_antigens.csv")
    n5 = int(nom.antigen.nunique())                          # surface-accessible (UniProt gate)
    n6 = int(nom.amplicon.nunique())                         # distinct amplicons carrying >=1

    _figure([n1, n2, n3, n4, n5, n6])
    ish.record("00b_target_funnel", dict(
        universe=n1, transmitted=n2, coelevated=n3, coelev_transmitted=n4,
        nominated_antigens=n5, nominated_amplicons=n6,
        n_cohorts_nominated=int(nom.cohort.nunique()),
        transmit_floor=TRANSMIT_MIN, fdr_alpha=cfg.FDR_ALPHA))
    print("funnel:", n1, n2, n3, n4, n5, n6, "| cohorts:", sorted(nom.cohort.unique()))

def _figure(counts):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    ish.apply_style(sizes=(9,8,7))
    n1,n2,n3,n4,n5,n6 = counts
    steps = [
     ("Genes with measured + predicted\ntransmissibility (CPTAC pan-cancer)", f"{n1:,} genes"),
     (f"Transmitted to protein\n(observed transmissibility \u2265 {TRANSMIT_MIN:.2f})", f"{n2:,} genes"),
     ("Co-elevated on a recurrent amplicon\n(Fisher FDR < 0.1, any cohort)", f"{n3:,} genes"),
     ("Co-elevated AND transmitted", f"{n4:,} genes"),
     ("Surface-accessible ectodomain\n(UniProt \u2265 50 aa, \u22651 TM)", f"{n5} antigens"),
     ("On distinct recurrent amplicons\n(pan-cancer nomination set)", f"{n6} amplicons"),
    ]
    fig, ax = plt.subplots(figsize=(8.6,6.8)); ax.axis("off")
    n=len(steps); ymax=n
    widths=np.linspace(1.0,0.42,n); cols=plt.cm.YlGnBu(np.linspace(0.35,0.85,n)); cx=0.5
    for i,((label,gc),w,c) in enumerate(zip(steps,widths,cols)):
        y=ymax-i-1; x0=cx-w/2
        ax.add_patch(FancyBboxPatch((x0,y+0.12),w,0.76,boxstyle="round,pad=0.01,rounding_size=0.03",
                     linewidth=1,edgecolor="white",facecolor=c))
        tc="white" if i>=3 else "#1a1a1a"
        ax.text(cx,y+0.60,label,ha="center",va="center",fontsize=8.8,color=tc,fontweight="bold")
        ax.text(cx,y+0.28,gc,ha="center",va="center",fontsize=11,color=tc)
        if i<n-1:
            ax.annotate("",xy=(cx,y+0.02),xytext=(cx,y+0.12),
                        arrowprops=dict(arrowstyle="-|>",color="#888",lw=1.4))
    ax.set_xlim(-0.05,1.15); ax.set_ylim(-0.2,ymax)
    ax.set_title("Target-nomination funnel: from pan-cancer transmissibility\n"
                 "to accessible multi-antigen ADC targets",
                 fontsize=10.5,fontweight="bold",pad=12,loc="center")
    ax.text(cx,-0.12,"Co-elevation tested per recurrent amplicon across the four proteome-supported "
            "cohorts (CCRCC, LSCC, LUAD, UCEC).\nNominated antigens span all four; single-cell "
            "co-detection is demonstrated where a malignant census slice exists (LUAD, LSCC, GBM).",
            ha="center",va="center",fontsize=6.4,color="#555",style="italic")
    fig.savefig(cfg.DIR_FIG / "fig2_target_funnel.png", dpi=200, bbox_inches="tight")

if __name__ == "__main__":
    main()