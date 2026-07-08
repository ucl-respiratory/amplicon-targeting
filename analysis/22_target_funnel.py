# =============================================================================
# 22_target_funnel.py  --  Fig 1: target nomination funnel, from genes on
# recurrent amplicons through to surface-accessible multi-antigen sets. Every
# count is COMPUTED from the co-elevation + dependency + topology pipeline
# (no hardcoded numbers), so the funnel tracks the data exactly.
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import _amplicon, _depmap, cnt_io
from cnt_io import savefig, record_values, cfg

def main():
    co = _amplicon.amplicon_coelevation()
    topo = cnt_io.load_topology()

    # ---- funnel counts (events = gene x type x amplicon; genes = unique) ----
    s1_g, s1_e = co.gene.nunique(), len(co)                                   # on recurrent amplicons
    sig = co[co.fdr < 0.1]
    s2_g, s2_e = sig.gene.nunique(), len(sig)                                 # reliably co-elevated
    surf = sig[(sig.is_surface == 1) | (sig.is_secreted == 1)]
    s3_g, s3_e = surf.gene.nunique(), len(surf)                              # surface / secreted
    strong = surf[surf.p_high_amp >= cfg.COELEV_MIN_FRAC]
    s4_g, s4_e = strong.gene.nunique(), len(strong)                         # elevated in >=50% amplified
    gs = _depmap.cotarget_dependency_annotated().drop_duplicates(["gene","tumor_code","arm"])
    s5_g, s5_e = gs.gene.nunique(), len(gs)                                  # DepMap dependency annotated
    s6_g = int(topo.surface_accessible.sum())                               # surface-accessible ectodomain

    steps = [
        ("Genes on recurrent amplicons\n(6 CPTAC cancer types)", f"{s1_g:,} genes", f"{s1_e:,} gene\u00d7amplicon tests"),
        ("Reliably co-elevated\n(Fisher FDR < 0.1)", "\u2014", f"{s2_e:,} events"),
        ("Surface / secreted\n(in-silico surfaceome)", f"{s3_g} genes", f"{s3_e} events"),
        ("Elevated in \u226550% of\namplified tumours", f"{s4_g} genes", f"{s4_e} events"),
        ("DepMap dependency\nannotation available", f"{s5_g} genes", f"{s5_e} events"),
        ("Surface-accessible ectodomain\n(UniProt-curated topology)", f"{s6_g} genes", "\u2014"),
    ]
    plt.rcParams.update({"font.size":9,"figure.dpi":150})
    fig, ax = plt.subplots(figsize=(8.4, 6.6)); ax.axis("off")
    n = len(steps); ymax = n
    widths = np.linspace(1.0, 0.42, n); cols = plt.cm.YlGnBu(np.linspace(0.35, 0.85, n)); cx = 0.5
    for i, ((label, gcount, ecount), w, c) in enumerate(zip(steps, widths, cols)):
        y = ymax - i - 1; x0 = cx - w/2
        ax.add_patch(FancyBboxPatch((x0, y+0.12), w, 0.76, boxstyle="round,pad=0.01,rounding_size=0.03",
                                    linewidth=1, edgecolor="white", facecolor=c))
        txtcol = "white" if i >= 3 else "#1a1a1a"
        ax.text(cx, y+0.62, label, ha="center", va="center", fontsize=8.8, color=txtcol, fontweight="bold")
        disp = gcount if gcount != "\u2014" else ecount
        ax.text(cx, y+0.30, disp, ha="center", va="center", fontsize=10.5, color=txtcol)
        if gcount != "\u2014" and ecount != "\u2014":
            ax.text(cx+w/2+0.02, y+0.5, ecount, ha="left", va="center", fontsize=7.2, color="#555", style="italic")
        if i < n-1:
            ax.annotate("", xy=(cx, y+0.02), xytext=(cx, y+0.12), arrowprops=dict(arrowstyle="-|>", color="#888", lw=1.4))
    ax.set_xlim(-0.05, 1.15); ax.set_ylim(0, ymax)
    ax.set_title("Target nomination funnel: from amplicon genes to accessible multi-antigen sets",
                 fontsize=10.5, fontweight="bold", pad=12, loc="center")
    ax.text(cx, -0.15, "Grey italics = gene\u00d7cancer-type events; a gene co-elevated on the same arm in two cancer types counts once here.",
            ha="center", va="top", fontsize=6.8, color="#666", transform=ax.transData)
    fig.tight_layout(); savefig(fig, "target_funnel.png"); plt.close(fig)

    record_values("22_target_funnel", {
        "s1_recurrent_amplicon": {"genes": int(s1_g), "tests": int(s1_e)},
        "s2_reliably_coelevated_events": int(s2_e), "s2_genes": int(s2_g),
        "s3_surface_secreted": {"genes": int(s3_g), "events": int(s3_e)},
        "s4_elevated_50pct": {"genes": int(s4_g), "events": int(s4_e)},
        "s5_depmap_annotated": {"genes": int(s5_g), "events": int(s5_e)},
        "s6_surface_accessible_genes": int(s6_g),
    })
    print(f"S1 {s1_g}g/{s1_e}t  S2 {s2_g}g/{s2_e}e  S3 {s3_g}g/{s3_e}e  "
          f"S4 {s4_g}g/{s4_e}e  S5 {s5_g}g/{s5_e}e  S6 {s6_g}g")

if __name__ == "__main__":
    main()
