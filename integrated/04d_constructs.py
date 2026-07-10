# =============================================================================
# 04d_constructs.py -- multivalent ADC constructs from pan-cancer nominations
# -----------------------------------------------------------------------------
# Assembles one multivalent construct per arm-amplicon that carries >=2
# nominated surface antigens (from 04c). Same-cell co-detection enrichment is
# computed where a malignant single-cell census slice exists (LUAD, LSCC);
# constructs on cohorts without a census slice (CCRCC, UCEC) are reported as
# NOMINATED but flagged not-yet-single-cell-tested -- the single-cell demo is
# currently lung-only because only LUAD/LSCC malignant slices were extracted.
#
# Emits: tables/adc_constructs.csv  (with 'single_cell_tested' flag)
# =============================================================================
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cnt_shared as ish
import config as cfg

def enrich(df, genes, seed=cfg.SEED, n_boot=1000, n_perm=1000):
    genes = [g for g in genes if g in df.columns]
    if len(genes) < 2: return None
    rng = np.random.default_rng(seed)
    D = (df[genes].values > 0).astype(int)
    codet = (D.sum(1) == len(genes)).mean()
    indep = np.prod(D.mean(0))
    enr = codet / indep if indep > 0 else np.nan
    donors = df["donor_id"].values; ud = np.unique(donors)
    idx_by_donor = {d: np.where(donors == d)[0] for d in ud}
    bs = []
    for _ in range(n_boot):
        samp = np.concatenate([idx_by_donor[d] for d in rng.choice(ud, len(ud), True)])
        Db = D[samp]; i = np.prod(Db.mean(0))
        bs.append((Db.sum(1) == len(genes)).mean() / i if i > 0 else np.nan)
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    ge = 0
    for _ in range(n_perm):
        Pm = np.column_stack([rng.permutation(D[:, j]) for j in range(len(genes))])
        ip = np.prod(Pm.mean(0))
        if (((Pm.sum(1) == len(genes)).mean() / ip) if ip > 0 else 0) >= enr: ge += 1
    return dict(enrich=round(enr, 2), ci_lo=round(lo, 2), ci_hi=round(hi, 2),
                perm_p=round((ge + 1) / (n_perm + 1), 4),
                n_cells=int(len(df)), n_donors=int(len(ud)))

COHORT_COLORS = {"LUAD":"#4C72B0","LSCC":"#55A868","CCRCC":"#C44E52","UCEC":"#8172B3"}

def _figure(nom, con):
    import matplotlib.pyplot as plt
    ish.apply_style(sizes=(9,8,7))
    atlas = pd.read_csv(cfg.GI_PATHS["atlas"])
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.1))

    # (a) nominated antigens on the transmissibility plane vs genome background
    ax = axes[0]
    ax.scatter(atlas.predicted_transmissibility_oof, atlas.observed_transmissibility,
               s=4, c="#d0d0d0", alpha=0.4, edgecolors="none", rasterized=True)
    for co, g in nom.groupby("cohort"):
        ax.scatter(g.pred_transmit, g.obs_transmit, s=34, c=COHORT_COLORS.get(co,"#333"),
                   edgecolors="k", linewidths=0.4, label=co)
    ax.axhline(0.40, ls="--", c="grey", lw=0.8)
    ax.set_xlabel("predicted transmissibility"); ax.set_ylabel("observed transmissibility")
    ax.set_title("Nominated antigens vs genome"); ax.legend(frameon=False, fontsize=7, title="cohort")
    ish.panel_letter(ax, "a")

    # (b) antigens per amplicon, colored by cohort
    ax = axes[1]
    order = (nom.groupby("amplicon").obs_transmit.max().sort_values().index.tolist())
    yl = []
    for i, amp in enumerate(order):
        g = nom[nom.amplicon == amp]; co = g.cohort.iloc[0]
        ax.barh(i, g.obs_transmit.max(), color=COHORT_COLORS.get(co,"#333"), alpha=0.85)
        ax.text(g.obs_transmit.max()+0.01, i, "+".join(g.antigen), va="center", fontsize=6)
        yl.append(amp)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(yl, fontsize=6)
    ax.set_xlabel("max observed transmissibility on amplicon")
    ax.set_title("Nominated antigens per amplicon (pan-cancer)")
    ish.panel_letter(ax, "b")

    # (c) construct co-detection: tested (filled) vs nominated-only (hatched)
    ax = axes[2]
    ct = con.dropna(subset=["enrich"]).sort_values("enrich")
    y = np.arange(len(ct))
    ax.barh(y, ct.enrich, xerr=[ct.enrich-ct.ci_lo, ct.ci_hi-ct.enrich],
            color=[COHORT_COLORS.get(c,"#333") for c in ct.cohort], alpha=0.85,
            error_kw=dict(lw=0.8, capsize=2))
    ax.axvline(1.0, ls="--", c="grey", lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(ct.construct, fontsize=6)
    ax.set_xlabel("same-cell co-detection enrichment (x)")
    n_nom = int((~con.single_cell_tested).sum())
    ax.set_title(f"Construct co-detection (LUAD/LSCC tested;\n{n_nom} nominated non-lung not shown)")
    ish.panel_letter(ax, "c")

    fig.tight_layout()
    fig.savefig(cfg.DIR_FIG / "fig5_surface_targets.png", dpi=200, bbox_inches="tight")

def main():
    nom = pd.read_csv(cfg.DIR_TAB / "adc_target_antigens.csv")
    census = {}
    for code in ("LUAD", "LSCC"):
        p = cfg.PATHS[f"cxg_{code.lower()}"]
        if p.exists(): census[code] = pd.read_parquet(p)

    rows = []
    for amp, sub in nom.groupby("amplicon"):
        genes = sorted(sub.antigen.unique()); cohort = sub.cohort.iloc[0]
        if len(genes) < 2:  # single-antigen amplicon: not a multivalent construct
            continue
        valence = {2:"bivalent",3:"trivalent"}.get(len(genes), f"{len(genes)}-valent")
        rec = dict(construct=f"{amp} ({'+'.join(genes)})", cohort=cohort,
                   amplicon=amp, antigens="+".join(genes), valence=valence, k=len(genes))
        df = census.get(cohort)
        if df is not None and all(g in df.columns for g in genes):
            e = enrich(df, genes)
            if e:
                rec.update(e); rec["single_cell_tested"] = True
        if "single_cell_tested" not in rec:
            rec.update(dict(enrich=np.nan, ci_lo=np.nan, ci_hi=np.nan, perm_p=np.nan,
                            n_cells=np.nan, n_donors=np.nan, single_cell_tested=False))
        rows.append(rec)

    con = pd.DataFrame(rows)
    con = con.sort_values(["single_cell_tested","enrich"], ascending=[False, False])
    con.to_csv(cfg.DIR_TAB / "adc_constructs.csv", index=False)

    _figure(nom, con)

    tested = con[con.single_cell_tested]
    print(f"constructs: {len(con)} total across {con.cohort.nunique()} cohorts")
    print(f"  single-cell tested (LUAD/LSCC): {len(tested)}")
    print(f"  nominated, not yet single-cell tested: {len(con)-len(tested)}")
    if len(tested):
        print("  tested enrichments:", dict(zip(tested.construct, tested.enrich)))

    ish.record("04d_constructs", {
        "n_constructs": int(len(con)),
        "n_tested": int(len(tested)),
        "n_nominated_only": int(len(con) - len(tested)),
        "tested_enrich_min": float(tested.enrich.min()) if len(tested) else None,
        "tested_enrich_max": float(tested.enrich.max()) if len(tested) else None,
        "cohorts_with_constructs": sorted(con.cohort.unique().tolist()),
        "cohorts_single_cell": ["LUAD","LSCC"],
    })
    return con

if __name__ == "__main__":
    main()
