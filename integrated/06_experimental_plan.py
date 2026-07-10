#!/usr/bin/env python3
# =============================================================================
# 06_experimental_plan.py  --  Figure 7: from in-silico nomination to the bench.
# -----------------------------------------------------------------------------
# A schematic that maps each in-silico output to the concrete wet-lab step that
# tests it, with explicit go/no-go criteria. The paper's forward plan, drawn so a
# clinical reader sees exactly what would be measured next and what would kill or
# advance a lead. No data dependency.
# Output: figures/fig7_experimental_plan.png
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import cnt_shared as ish
ish.apply_style(sizes=(9, 8, 7))

STEPS = [
    ("In-silico nomination",
     "Amplicon co-target set\n(measured + predicted,\ntagged by confidence tier)",
     "#34495e"),
    ("1 · Surface-protein\nco-expression",
     "Flow / CyTOF on amplicon+ lines\n(named panel below) vs amplicon−.\nGO: all antigens surface-positive\nin >50% of amplicon+ cells",
     "#2c6fbb"),
    ("2 · Ectodomain binding\n(marginal antigens)",
     "Cell-surface ELISA / BLI on\nmarginal ectodomains (e.g. PLSCR1).\nGO: measurable antibody binding\nto native surface epitope",
     "#e67e22"),
    ("3 · AND-gate construct",
     "Bispecific / avidity ADC;\nco-engagement required for stable\nbinding. GO: tumour killing needs\nboth antigens; single-antigen spares",
     "#6c5b9c"),
    ("4 · Selectivity window",
     "Amplicon+ tumour vs matched\nnormal (low-antigen tissue).\nGO: therapeutic index > single-\nantigen ADC at equal payload",
     "#27ae60"),
]

# Cell-line entry points are amplicon-positive lines identified from DepMap 23Q4
# copy number (all construct anchor genes at relative CN >= 1.4 in the matching
# lineage); amplicon-negative LUAD lines serve as the single-antigen control arm.
PANELS = [
    ("LUAD 7p (EGFR+ITGB8+\nTSPAN13+TTYH3)",
     "amplicon+ : NCI-H3255, NCI-H1838, HCC4006\ncontrol (amplicon\u2212): NCI-H441, NCI-H358\nsame-cell co-detection confirmed"),
    ("LUAD 1q (ADAM15+NCSTN)\n& LSCC 1q (HSD17B7+\nMPZL1+NCSTN)",
     "amplicon+ : NCI-H23, HCC1833 (LUAD);\nLUDLU-1, GT3TKB (LSCC)\nco-detection confirmed in both"),
    ("LSCC 3q (ATP11B+TFRC)\n& predicted-only extension",
     "amplicon+ : HCC95, LC-1/sq (squamous 3q,\nSOX2 locus). Predicted-only sets in a\nproteomics-poor type extend arm 1 first"),
]

def main():
    fig = plt.figure(figsize=(13.5, 7.2))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 13.5); ax.set_ylim(0, 7.2); ax.axis("off")
    ax.text(6.75, 6.95, "From in-silico nomination to the bench",
            ha="center", fontsize=14, fontweight="bold")
    ax.text(6.75, 6.6, "each in-silico output maps to one wet-lab test with an explicit go/no-go criterion",
            ha="center", fontsize=9.5, style="italic", color="#555")

    # flow row
    n = len(STEPS); x0 = 0.35; w = 2.35; gap = 0.28; y = 3.9; h = 1.85
    boxes = []
    for i, (title, body, color) in enumerate(STEPS):
        x = x0 + i*(w+gap)
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.10",
                                    lw=1.6, edgecolor=color, facecolor=color, alpha=0.12, zorder=2))
        ax.text(x+w/2, y+h-0.28, title, ha="center", va="top", fontsize=9,
                fontweight="bold", color=color)
        ax.text(x+w/2, y+h-0.72, body, ha="center", va="top", fontsize=7.3, color="#222")
        boxes.append((x, x+w))
    for i in range(n-1):
        ax.add_patch(FancyArrowPatch((boxes[i][1]+0.02, y+h/2), (boxes[i+1][0]-0.02, y+h/2),
                     arrowstyle="-|>", mutation_scale=16, lw=1.8, color="#888", zorder=1))

    # go/no-go strip
    ax.text(6.75, 3.5, "A no-go at any step returns the set to nomination with the failing "
            "antigen dropped or the set re-formed — the tiering makes that cheap.",
            ha="center", fontsize=8.2, style="italic", color="#666")

    # named-panel row
    ax.text(0.35, 2.9, "Named cell-line / cohort panels (arm 1 entry points):",
            fontsize=9.5, fontweight="bold", color="#34495e")
    pw = 4.15; py = 0.6; ph = 1.9
    for j, (name, desc) in enumerate(PANELS):
        px = 0.35 + j*(pw+0.3)
        ax.add_patch(FancyBboxPatch((px, py), pw, ph, boxstyle="round,pad=0.04,rounding_size=0.08",
                                    lw=1.3, edgecolor="#34495e", facecolor="#f4f6f8", zorder=2))
        ax.text(px+pw/2, py+ph-0.3, name, ha="center", va="top", fontsize=9,
                fontweight="bold", color="#2c6fbb")
        ax.text(px+pw/2, py+ph-0.75, desc, ha="center", va="top", fontsize=7.6, color="#333")

    fig.savefig(ish.cfg.DIR_FIG / "fig7_experimental_plan.png", dpi=200, bbox_inches="tight")
    print("saved fig7_experimental_plan.png")

if __name__ == "__main__":
    main()
