# =============================================================================
# 30_normal_singlecell_threshold.py  --  Fig 5: magnitude-aware normal-tissue
# co-expression. Normal-cell binding of an AND-gate is set by the WEAKEST
# (limiting) antigen; raising the per-antigen threshold from detection (10 nTPM)
# to binding-relevant (25/50 nTPM) collapses the normal-cell burden. Also shows
# the amplified-tumour vs normal-tissue RNA fold change (median ~64x).
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cnt_io
from cnt_io import DIR_TAB, savefig, record_values, cfg

SETS = {"LUAD_1q":["ADAM15","CD46","EFNA1","MUC1","NCSTN","XPR1"],
        "LUAD_7p":["DAGLB","EGFR","ITGB8","TSPAN13"],
        "LSCC_1q":["F11R","HSD17B7","NCSTN"],
        "LSCC_20q":["GGT7","SDC4","TM9SF4"],
        "PDA_1q":["F11R","MUC1","NCSTN"]}
THRESH = [10, 25, 50, 100]

def main():
    sc = cnt_io.load_hpa_singlecell()
    scw = sc.pivot_table(index="Gene name", columns="Cell type", values="nTPM", aggfunc="max")
    main = cnt_io.load_str_omic()[["gene","caseid","tumor_code","cn_adjusted","rna","gtex","prot.rel.tissue"]]

    def limiting_count(genes, T):
        p = [g for g in genes if g in scw.index]
        if not p: return np.nan
        return int((scw.loc[p].fillna(0).min(axis=0) >= T).sum())

    def tumor_and(genes, tc):
        d = main[(main.tumor_code==tc)&(main.gene.isin(genes))]
        pp = d.pivot_table(index="caseid",columns="gene",values="prot.rel.tissue",aggfunc="first")
        pc = d.pivot_table(index="caseid",columns="gene",values="cn_adjusted",aggfunc="first")
        gp = [g for g in genes if g in pp.columns and g in pc.columns]
        common = pp[gp].dropna().index.intersection(pc[gp].dropna().index)
        if not len(common): return np.nan, 0
        pp = pp.loc[common,gp]; pc = pc.loc[common,gp]
        amp = pc.mean(axis=1) >= cfg.AMP_THRESHOLD; high = pp > cfg.REL_TISSUE_HI
        return (float(high.all(axis=1)[amp].mean()) if amp.sum() else np.nan, int(amp.sum()))

    rows = []
    for sn, genes in SETS.items():
        tc = sn.split("_")[0]; tand, namp = tumor_and(genes, tc)
        meds = {g: scw.loc[g].median() for g in genes if g in scw.index}
        lim = min(meds, key=meds.get)
        r = {"set":sn,"module":"+".join(genes),"valence":len(genes),"limiting_antigen":lim,
             "limiting_normal_median_nTPM":round(meds[lim],1),
             "tumour_AND_amp":round(tand,3) if not np.isnan(tand) else None,"n_amp":namp}
        for T in THRESH: r[f"normal_celltypes_minGE{T}"] = limiting_count(genes, T)
        rows.append(r)
    mag = pd.DataFrame(rows)
    mag.to_csv(DIR_TAB / "normal_singlecell_threshold.csv", index=False)

    fc = []
    for sn, genes in SETS.items():
        tc = sn.split("_")[0]
        for g in genes:
            d = main[(main.tumor_code==tc)&(main.gene==g)]; amp = d[d.cn_adjusted >= cfg.AMP_THRESHOLD]
            if len(amp) < 5: continue
            tum = amp.rna.median(); nrm = d.gtex.median()
            fc.append({"set":sn,"gene":g,"tumRNA_amp":round(tum,0),"gtexRNA_normal":round(nrm,1),
                       "fold_tum_over_normal":round(tum/nrm,1) if nrm and nrm>0 else None,
                       "normal_sc_median_nTPM":round(scw.loc[g].median(),1) if g in scw.index else None})
    fcdf = pd.DataFrame(fc)
    median_fold = float(fcdf.fold_tum_over_normal.median())

    plt.rcParams.update({"font.size":8,"axes.titlesize":9,"axes.labelsize":8,"legend.fontsize":7,
        "xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,
        "axes.titlelocation":"left"})
    setcols = {"LUAD_1q":"#4c72b0","LUAD_7p":"#55a868","LSCC_1q":"#c44e52","LSCC_20q":"#dd8452","PDA_1q":"#8172b3"}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7))
    ax = axes[0]
    for _, r in mag.iterrows():
        y = [r[f"normal_celltypes_minGE{T}"] for T in THRESH]
        ax.plot(THRESH, y, "-o", ms=4, lw=1.4, color=setcols[r["set"]], label=f"{r['set']} ({r.limiting_antigen})")
    ax.set_xlabel("Binding-relevant threshold on limiting antigen (nTPM)")
    ax.set_ylabel("Normal cell types where AND-gate binds\n(weakest member \u2265 threshold)")
    ax.set_title("A  Normal-cell burden collapses when threshold\nreflects binding, not detection")
    ax.set_xticks(THRESH); ax.legend(frameon=False, fontsize=6.3, loc="upper right")
    ax = axes[1]
    fcs = fcdf.sort_values("fold_tum_over_normal")
    islim = fcs.gene.isin(["MUC1","DAGLB","HSD17B7","GGT7"]); cols = ["#c44e52" if l else "#b0b0b0" for l in islim]
    ax.barh(range(len(fcs)), np.log10(fcs.fold_tum_over_normal), color=cols)
    ax.set_yticks(range(len(fcs))); ax.set_yticklabels(fcs.gene, fontsize=6.2)
    ax.axvline(0, color="k", lw=0.8); ax.axvline(np.log10(10), color="grey", ls=":", lw=0.8)
    ax.set_xlabel("log10(amplified-tumour RNA / normal-tissue RNA)")
    ax.set_title(f"B  Every antigen: amplified tumour \u226b normal\n(median {median_fold:.0f}\u00d7; red = limiting antigens)")
    ax.set_xlim(0, 3.2)
    for i, (_, r) in enumerate(fcs.iterrows()):
        ax.text(np.log10(r.fold_tum_over_normal)+0.03, i, f"{r.fold_tum_over_normal:.0f}\u00d7", va="center", fontsize=5.8)
    ax = axes[2]
    x = np.arange(len(mag)); w = 0.38
    det = [mag.iloc[i]["normal_celltypes_minGE10"] for i in range(len(mag))]
    binr = [mag.iloc[i]["normal_celltypes_minGE25"] for i in range(len(mag))]
    ax.bar(x-w/2, det, w, label="detection (\u226510 nTPM)", color="#bdbdbd")
    ax.bar(x+w/2, binr, w, label="binding-relevant (\u226525 nTPM)", color="#4c72b0")
    ax.set_xticks(x); ax.set_xticklabels(mag.set, rotation=30, ha="right", fontsize=6.5)
    ax.set_ylabel("Normal cell types co-expressing whole set")
    ax.set_title("C  Detection threshold overstated burden:\n\u226510 vs \u226525 nTPM on the limiting member")
    ax.legend(frameon=False, loc="upper right")
    for i in range(len(mag)):
        ax.text(x[i]-w/2, det[i]+0.5, str(det[i]), ha="center", fontsize=6)
        ax.text(x[i]+w/2, binr[i]+0.5, str(binr[i]), ha="center", fontsize=6)
    fig.suptitle("Magnitude-aware single-cell co-expression: normal binding is gated by the weakest antigen",
                 fontsize=9.5, y=1.02)
    fig.tight_layout(); savefig(fig, "singlecell_valence_purity.png"); plt.close(fig)

    record_values("30_normal_singlecell_threshold", {
        "median_fold_tumour_over_normal": round(median_fold,1),
        "sets": {r.set: {"limiting_antigen": r.limiting_antigen,
                 "limiting_normal_median_nTPM": r.limiting_normal_median_nTPM,
                 "normal_celltypes_ge10": int(r.normal_celltypes_minGE10),
                 "normal_celltypes_ge25": int(r.normal_celltypes_minGE25),
                 "normal_celltypes_ge50": int(r.normal_celltypes_minGE50)} for r in mag.itertuples()},
    })
    print(mag[["set","limiting_antigen","limiting_normal_median_nTPM","normal_celltypes_minGE10","normal_celltypes_minGE25","normal_celltypes_minGE50"]].to_string(index=False))
    print(f"median fold tumour/normal = {median_fold:.1f}x")

if __name__ == "__main__":
    main()
