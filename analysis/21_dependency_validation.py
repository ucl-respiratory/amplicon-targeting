# =============================================================================
# 21_dependency_validation.py  --  Fig 3 + core nomination tables.
# Cross-cohort DepMap validation (CCLE proteomics), CRISPR-dependency overlay
# separating dispensable passengers from essential drivers, and the clean
# multi-target set list. Writes depmap_validation.csv, gene_dependency_summary.csv,
# cotarget_dependency_annotated.csv.
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import _depmap
from cnt_io import DIR_TAB, savefig, record_values

def main():
    gene_dep = _depmap.gene_dependency_summary()
    val = _depmap.depmap_validation()
    gs = _depmap.cotarget_dependency_annotated()
    gene_dep.to_csv(DIR_TAB / "gene_dependency_summary.csv", index=False)
    val.to_csv(DIR_TAB / "depmap_validation.csv", index=False)
    gs.to_csv(DIR_TAB / "cotarget_dependency_annotated.csv", index=False)

    # --- headline validation values ---
    n_testable = len(val)
    frac_up = float((val.delta > 0).mean())
    med_delta = float(val.delta.median())
    q1, q3 = float(val.delta.quantile(.25)), float(val.delta.quantile(.75))
    n_fdr10 = int((val.bh_fdr < 0.10).sum()); n_fdr05 = int((val.bh_fdr < 0.05).sum())
    n_big = int((val.delta > 0.5).sum())
    # surface-accessible subset = validated targets whose gene has an accessible
    # ectodomain by UniProt-curated topology (surface_accessible == True).
    import cnt_io as _cio
    _topo = _cio.load_topology()
    _acc_genes = set(_topo.loc[_topo.surface_accessible, "gene"])
    med_delta_surf = float(val[val.gene.isin(_acc_genes)].delta.median())

    # nominated surface targets are counted at the (gene x type x arm) EVENT level
    # (matching the manuscript: "75 nominated surface targets, 70 dispensable, 4 essential")
    ev = gs.drop_duplicates(["gene","tumor_code","arm"])
    n_events = len(ev); n_genes = ev.gene.nunique()
    cls = ev.dependency_class.value_counts()
    n_passenger = int(cls.get("dispensable_passenger", 0))
    n_essential = int(cls.get("essential_driver-like", 0))
    n_unknown   = int(cls.get("unknown", 0))
    # tgt = events with a dependency value, for the scatter panels
    tgt = ev.dropna(subset=["mean_dependency"])
    n_targets = n_events

    setcols = gs.groupby(["tumor_code","arm"]).apply(
        lambda x: pd.Series({
            "n_total": x.gene.nunique(),
            "n_passenger": x[x.dependency_class=="dispensable_passenger"].gene.nunique(),
            "passengers": ";".join(sorted(set(x[x.dependency_class=="dispensable_passenger"].gene))),
            "essential_members": ";".join(sorted(set(x[x.dependency_class=="essential_driver-like"].gene))),
            "mean_dep": x.mean_dependency.mean()}), include_groups=False).reset_index()
    clean_sets = setcols[setcols.n_passenger >= 2].sort_values("n_passenger", ascending=False)
    clean_sets.to_csv(DIR_TAB / "clean_multitarget_sets.csv", index=False)

    plt.rcParams.update({"font.size":8,"axes.titlesize":8.5,"axes.labelsize":8,"legend.fontsize":7,
        "xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,
        "axes.titlelocation":"left"})
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4))
    ax = axes[0]
    ax.hist(gene_dep[gene_dep.common_essential==1].mean_dependency, bins=30, alpha=0.5, density=True, color="#c0504d", label="common-essential")
    ax.hist(gene_dep[gene_dep.nonessential==1].mean_dependency, bins=30, alpha=0.5, density=True, color="#8aa1b4", label="nonessential")
    ax.scatter(tgt.mean_dependency, np.full(len(tgt), -0.4), s=14, color="#5a9367", zorder=5, label="nominated targets", clip_on=False)
    ax.axvline(0.5, color="k", lw=0.8, ls="--")
    ax.set_xlabel("Mean CRISPR dependency"); ax.set_ylabel("Density")
    ax.set_title(f"A  {n_passenger}/{n_events} targets are dispensable"); ax.legend(frameon=False, fontsize=6.3); ax.set_ylim(bottom=-0.8)
    ax = axes[1]
    cc = np.where(tgt.dependency_class=="essential_driver-like", "#c0504d", "#5a9367")
    ax.scatter(tgt.mean_dependency, tgt.p_high_amp, s=22, c=cc, alpha=0.75, edgecolors="none")
    ax.axvline(0.5, color="k", lw=0.8, ls="--")
    ax.set_xlabel("Mean CRISPR dependency"); ax.set_ylabel("Co-elevation reliability P(high|amp)")
    ax.set_title("B  ADC-preferred = low dep, high reliability")
    for _, r in tgt[tgt.dependency_class=="essential_driver-like"].iterrows():
        ax.annotate(r.gene, (r.mean_dependency, r.p_high_amp), fontsize=6, color="#8a2820", xytext=(2,2), textcoords="offset points")
    ax = axes[2]
    cs = clean_sets.head(10).iloc[::-1]
    ax.barh(range(len(cs)), cs.n_passenger, color="#5a9367")
    ax.set_yticks(range(len(cs))); ax.set_yticklabels([f"{r.tumor_code} {r.arm}" for r in cs.itertuples()], fontsize=6.5)
    ax.set_xlabel("# passenger surface targets (ADC-preferred)"); ax.set_title("C  Clean multi-target ADC sets")
    for i, r in enumerate(cs.itertuples()): ax.text(r.n_passenger+0.05, i, str(r.n_passenger), va="center", fontsize=6.5)
    fig.suptitle("CRISPR-dependency overlay: separating dispensable passengers (preferred ADC targets) from essential driver-like genes", fontsize=9, y=1.02)
    fig.tight_layout(); savefig(fig, "dependency_overlay.png"); plt.close(fig)

    record_values("21_dependency_validation", {
        "depmap_n_testable": n_testable, "depmap_frac_up": round(frac_up,4),
        "depmap_median_delta_log2": round(med_delta,4),
        "depmap_delta_iqr": [round(q1,4), round(q3,4)],
        "depmap_median_delta_surface": round(med_delta_surf,4),
        "depmap_n_fdr10": n_fdr10, "depmap_n_fdr05": n_fdr05, "depmap_n_delta_gt0.5": n_big,
        "n_nominated_surface_events": n_events, "n_nominated_unique_genes": int(n_genes),
        "n_dispensable_passenger": n_passenger, "n_essential_driver_like": n_essential,
        "n_unknown_dependency": n_unknown,
        "n_clean_multitarget_sets": int(len(clean_sets)),
    })
    print(f"testable={n_testable} frac_up={frac_up:.2%} med_delta={med_delta:.3f} IQR=[{q1:.2f},{q3:.2f}]")
    print(f"FDR<0.1={n_fdr10} FDR<0.05={n_fdr05} delta>0.5={n_big} med_delta_surface={med_delta_surf:.3f}")
    print(f"events={n_events} genes={n_genes} passenger={n_passenger} essential={n_essential} unknown={n_unknown}")
    print(f"essential genes: {sorted(tgt[tgt.dependency_class=='essential_driver-like'].gene.unique())}")

if __name__ == "__main__":
    main()
