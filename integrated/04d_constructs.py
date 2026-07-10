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

def _depth_bins(nnz, n_bins=cfg.SC_DEPTH_BINS):
    """Quantile (decile) bins on per-cell detected-gene count (sequencing depth)."""
    qs = np.quantile(nnz, np.linspace(0, 1, n_bins + 1)); qs[0] -= 1; qs[-1] += 1
    return np.digitize(nnz, qs[1:-1])

def _perm_codet_depth(D, binid, rng):
    """One depth-stratified permutation: within each depth bin, shuffle each
    antigen's detection column independently. Preserves BOTH per-gene marginal
    detection rate AND per-cell sequencing depth, so the null is co-detection
    expected from depth structure alone."""
    Dp = np.empty_like(D); k = D.shape[1]
    for b in np.unique(binid):
        ix = np.where(binid == b)[0]
        if len(ix) == 0: continue
        for j in range(k):
            Dp[ix, j] = D[ix[rng.permutation(len(ix))], j]
    return (Dp.sum(1) == k).mean()

def enrich(df, genes, seed=cfg.SEED, n_boot=1000, n_perm=1000):
    """Same-cell co-detection enrichment against a DEPTH-STRATIFIED null.

    Co-detection of any gene pair is inflated by per-cell sequencing depth
    (deep cells detect nearly everything), so a marginal-only column shuffle
    tests against an independence null that depth alone would violate. We
    instead permute within per-cell depth (nnz) deciles, so the reported
    enrichment is observed co-detection divided by the co-detection expected
    from depth structure, and the permutation p asks whether observed exceeds
    that depth-driven expectation."""
    genes = [g for g in genes if g in df.columns]
    if len(genes) < 2: return None
    rng = np.random.default_rng(seed)
    D = (df[genes].values > 0).astype(int)
    k = len(genes)
    codet = (D.sum(1) == k).mean()
    nnz = df["nnz"].values.astype(float)
    binid = _depth_bins(nnz)
    perm = np.array([_perm_codet_depth(D, binid, rng) for _ in range(n_perm)])
    exp_depth = perm.mean()
    enr = codet / exp_depth if exp_depth > 0 else np.nan
    pval = (np.sum(perm >= codet) + 1) / (n_perm + 1)
    # donor-block bootstrap of the depth-corrected enrichment
    donors = df["donor_id"].values; ud = np.unique(donors)
    idx_by_donor = {d: np.where(donors == d)[0] for d in ud}
    bs = []
    for _ in range(n_boot):
        samp = np.concatenate([idx_by_donor[d] for d in rng.choice(ud, len(ud), True)])
        Db = D[samp]; nb = nnz[samp]
        ob = (Db.sum(1) == k).mean()
        eb = _perm_codet_depth(Db, _depth_bins(nb), rng)
        if eb > 0: bs.append(ob / eb)
    lo, hi = np.nanpercentile(bs, [2.5, 97.5]) if bs else (np.nan, np.nan)
    return dict(enrich=round(enr, 2), ci_lo=round(lo, 2), ci_hi=round(hi, 2),
                perm_p=round(pval, 4),
                enrich_marginal=round(codet / np.prod(D.mean(0)), 2) if np.prod(D.mean(0)) > 0 else np.nan,
                n_cells=int(len(df)), n_donors=int(len(ud)))

COHORT_COLORS = {"LUAD":"#4C72B0","LSCC":"#55A868","CCRCC":"#C44E52","UCEC":"#8172B3"}

def _figure(nom, con):
    import matplotlib.pyplot as plt
    ish.apply_style(sizes=(9,8,7))
    atlas = pd.read_csv(cfg.GI_PATHS["atlas"])
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.1))

    # Panels (a) and (b) show the MEASURED nomination plane (observed
    # transmissibility from proteome); GBM is prediction-only (no proteome) and
    # therefore has no observed value to plot here. GBM enters only in panel (c),
    # the single-cell co-detection test, and in Supplementary Table S4.
    nom_m = nom[nom.get("evidence", "measured") == "measured"] if "evidence" in nom else nom

    # (a) nominated antigens on the transmissibility plane vs genome background
    ax = axes[0]
    ax.scatter(atlas.predicted_transmissibility_oof, atlas.observed_transmissibility,
               s=4, c="#d0d0d0", alpha=0.4, edgecolors="none", rasterized=True)
    for co, g in nom_m.groupby("cohort"):
        ax.scatter(g.pred_transmit, g.obs_transmit, s=34, c=COHORT_COLORS.get(co,"#333"),
                   edgecolors="k", linewidths=0.4, label=co)
    ax.axhline(0.40, ls="--", c="grey", lw=0.8)
    ax.set_xlabel("predicted transmissibility"); ax.set_ylabel("observed transmissibility")
    ax.set_title("Nominated antigens vs genome"); ax.legend(frameon=False, fontsize=7, title="cohort")
    ish.panel_letter(ax, "a")

    # (b) antigens per amplicon, colored by cohort
    ax = axes[1]
    order = (nom_m.groupby("amplicon").obs_transmit.max().sort_values().index.tolist())
    yl = []
    for i, amp in enumerate(order):
        g = nom_m[nom_m.amplicon == amp]; co = g.cohort.iloc[0]
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
    _sc = "/".join(sorted(con.loc[con.single_cell_tested, "cohort"].unique()))
    ax.set_title(f"Construct co-detection ({_sc} tested;\n{n_nom} nominated, no single-cell slice)")
    ish.panel_letter(ax, "c")

    fig.tight_layout()
    fig.savefig(cfg.DIR_FIG / "fig5_surface_targets.png", dpi=200, bbox_inches="tight")

# Cohorts with a malignant single-cell census slice. LUAD/LSCC carry measured
# nominations; GBM is prediction-only (no CPTAC GBM proteome) with a large
# malignant slice, and demonstrates the full predict->gate->single-cell-verify
# path with no protein data. CCRCC/UCEC/PDA have no malignant primary cells in
# the census and cannot be tested.
SC_COHORTS = ("LUAD", "LSCC", "GBM")

def _load_census():
    census = {}
    for code in SC_COHORTS:
        p = cfg.PATHS.get(f"cxg_{code.lower()}")
        if p is not None and p.exists():
            census[code] = pd.read_parquet(p)
    return census

def _gbm_prediction_constructs():
    """One prediction-only construct per GBM arm carrying >=2 nominated antigens
    (from the M2 table). Returns a DataFrame shaped like the measured nominations
    (cohort/amplicon/antigen) so it flows through the same construct assembly."""
    f = cfg.DIR_TAB / "m2_gbm_prediction_only.csv"
    if not f.exists():
        return pd.DataFrame()
    g = pd.read_csv(f)
    g = g.rename(columns={"pred_transmit": "obs_transmit"})  # no measured value; carry predicted in same slot
    g["cohort"] = "GBM"; g["amplicon"] = "GBM_" + g["arm"].astype(str)
    return g[["cohort", "amplicon", "antigen", "obs_transmit"]]

def main():
    nom = pd.read_csv(cfg.DIR_TAB / "adc_target_antigens.csv")
    nom["evidence"] = "measured"
    gbm = _gbm_prediction_constructs()
    if len(gbm):
        gbm["evidence"] = "prediction-only"
        nom = pd.concat([nom, gbm], ignore_index=True)
    census = _load_census()

    rows = []
    for amp, sub in nom.groupby("amplicon"):
        genes = sorted(sub.antigen.unique()); cohort = sub.cohort.iloc[0]
        if len(genes) < 2:  # single-antigen amplicon: not a multivalent construct
            continue
        valence = {2:"bivalent",3:"trivalent"}.get(len(genes), f"{len(genes)}-valent")
        rec = dict(construct=f"{amp} ({'+'.join(genes)})", cohort=cohort,
                   amplicon=amp, antigens="+".join(genes), valence=valence, k=len(genes),
                   evidence=sub.evidence.iloc[0])
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
    sc_present = sorted(census.keys())
    print(f"constructs: {len(con)} total across {con.cohort.nunique()} cohorts")
    print(f"  single-cell tested ({'/'.join(sc_present)}): {len(tested)}")
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
        "cohorts_single_cell": sc_present,
    })
    return con

if __name__ == "__main__":
    main()
