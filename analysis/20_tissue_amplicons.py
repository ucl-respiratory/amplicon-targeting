# =============================================================================
# 20_tissue_amplicons.py  --  Fig 2: tissue-specific amplicon landscape
# (fraction of samples with each chromosome arm amplified, per cancer type).
# =============================================================================
import sys, re
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cnt_io
from cnt_io import DIR_TAB, savefig, record_values, cfg

def main():
    df = cnt_io.load_str_omic()
    d = df.dropna(subset=["cn_adjusted","cytogenetic_location"]).copy()
    d["amplified"] = (d.cn_adjusted >= cfg.AMP_THRESHOLD).astype(int)

    band_samp = d.groupby(["tumor_code","cytogenetic_location","caseid"]).amplified.mean().reset_index()
    band_samp["band_amp"] = (band_samp.amplified >= cfg.COELEV_MIN_FRAC).astype(int)
    peaks = band_samp.groupby(["tumor_code","cytogenetic_location"]).agg(
        freq=("band_amp","mean"), n_amp=("band_amp","sum"), n=("band_amp","size")).reset_index()
    peaks["arm"] = peaks.cytogenetic_location.str.extract(r"^(\d+[pq])")
    peaks.to_csv(DIR_TAB / "tissue_amplicon_peaks.csv", index=False)

    arm_samp = d.groupby(["tumor_code","arm","caseid"]).amplified.mean().reset_index()
    arm_freq = arm_samp.groupby(["tumor_code","arm"]).apply(
        lambda x: (x.amplified >= cfg.COELEV_MIN_FRAC).mean(), include_groups=False).reset_index(name="freq")
    arm_mat = arm_freq.pivot(index="arm", columns="tumor_code", values="freq").fillna(0)
    def armkey(a):
        m = re.match(r"(\d+)([pq])", a); return (int(m.group(1)), 0 if m.group(2)=="p" else 1) if m else (99,0)
    arm_mat = arm_mat.loc[sorted(arm_mat.index, key=armkey)]

    plt.rcParams.update({"font.size":8,"axes.titlesize":9,"axes.labelsize":8,"xtick.labelsize":7,
        "ytick.labelsize":6.5,"axes.spines.top":False,"axes.spines.right":False})
    fig, ax = plt.subplots(figsize=(6.2, 8))
    im = ax.imshow(arm_mat.values, aspect="auto", cmap="OrRd", vmin=0, vmax=0.7)
    ax.set_xticks(range(len(arm_mat.columns))); ax.set_xticklabels(arm_mat.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(arm_mat.index))); ax.set_yticklabels(arm_mat.index)
    ax.set_xlabel("Tumour type"); ax.set_ylabel("Chromosome arm")
    ax.set_title("Tissue-specific amplicon landscape\n(fraction of samples with arm amplified)", fontsize=9, loc="left")
    known = {"7p":"GBM","7q":"GBM","3q":"LSCC","5q":"CCRCC","5p":"LSCC/LUAD"}
    for arm_, tis in known.items():
        if arm_ in arm_mat.index:
            yi = list(arm_mat.index).index(arm_); ax.text(len(arm_mat.columns)-0.3, yi, f"\u2190 {tis}", va="center", fontsize=6, color="#333")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.16); cbar.set_label("amplification frequency", fontsize=7)
    fig.tight_layout(); savefig(fig, "tissue_amplicon_heatmap.png"); plt.close(fig)

    n_recurrent = int(((peaks.freq >= 0.15) & (peaks.n_amp >= 8)).sum())
    record_values("20_tissue_amplicons", {
        "n_recurrent_amplicon_events": n_recurrent,
        "n_tumor_types": int(d.tumor_code.nunique()),
        "top_arm_by_type": {tc: arm_mat[tc].idxmax() for tc in arm_mat.columns},
    })
    print(f"recurrent amplicon (type,band) events: {n_recurrent}")
    print("top arm by type:", {tc: arm_mat[tc].idxmax() for tc in arm_mat.columns})

if __name__ == "__main__":
    main()
