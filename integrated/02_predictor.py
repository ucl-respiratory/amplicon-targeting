#!/usr/bin/env python3
# =============================================================================
# 02_predictor.py  --  Predicting transmissibility WITHOUT protein data.
# -----------------------------------------------------------------------------
# Trains a gradient-boosted predictor of per-gene copy-number-to-protein
# transmissibility from gene properties alone (dosage sensitivity, protein
# biophysics, mRNA, evolution, complex membership, expression breadth, network
# centrality). NO protein-derived feature is a predictor; the outcome is the
# protein-derived transmissibility. Leave-gene-out out-of-fold prediction gives
# the honest generalisation the extrapolation claim rests on.
#
# It also reproduces the two gene-intrinsic controls that license using the
# prediction as a genome-wide prior:
#   - positional control: leave-chromosome-arm-out vs leave-gene-out (no same-arm
#     neighbour leaks into a gene's prediction)  [needs arm annotation]
#   - transfer: cross-lineage rank concordance of predicted transmissibility
#     (from the committed transfer table)
#
# Inputs (committed GI exports; a fresh data_download+refit would regenerate them):
#   feature table (gene properties + transmissibility outcome), transfer table.
# Outputs: figures/fig3_predictor.png, tables/predictor_oof.csv
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.stats import spearmanr, kendalltau
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
import xgboost as xgb
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cnt_shared as ish
from cnt_shared import cfg
ish.apply_style(sizes=(9, 8, 7))

# Gene-property feature columns (NO protein-derived feature). Flatten the config
# groups and keep those present in the feature table.
FEATURE_COLS = [c for grp in cfg.PREDICTOR_GROUPS.values() for c in grp]

def _model():
    return xgb.XGBRegressor(n_estimators=600, max_depth=4, learning_rate=0.03,
                            subsample=0.8, colsample_bytree=0.8, n_jobs=1,
                            random_state=cfg.SEED, tree_method="hist")

def leave_gene_out(ft, feats):
    """5-fold out-of-fold prediction (genes grouped only by fold; no gene in its
    own training set). Returns oof predictions aligned to ft.gene."""
    X = ft[feats].values; y = ft["transmissibility"].values
    oof = np.zeros(len(y))
    kf = KFold(5, shuffle=True, random_state=cfg.SEED)
    for tr, te in kf.split(X):
        m = _model().fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
    return oof

def leave_arm_out(ft, feats, arm_col):
    """Hold out whole chromosome arms in turn: no gene is predicted using a
    same-arm neighbour. Positional control."""
    X = ft[feats].values; y = ft["transmissibility"].values
    oof = np.full(len(y), np.nan)
    arms = ft[arm_col].fillna("NA").values
    for arm in np.unique(arms):
        te = arms == arm; tr = ~te
        if te.sum() == 0 or tr.sum() < 100: continue
        m = _model().fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
    return oof

def main():
    ft = ish.load_feature_table()
    feats = [c for c in FEATURE_COLS if c in ft.columns]
    ft = ft.dropna(subset=["transmissibility"]).copy()
    # impute feature NaNs with column medians (predictor is gene-property only)
    ft[feats] = ft[feats].apply(lambda s: s.fillna(s.median()))

    oof = leave_gene_out(ft, feats)
    ft["predicted_oof"] = oof
    rho = spearmanr(ft.transmissibility, oof).correlation
    r2 = r2_score(ft.transmissibility, oof)

    # positional control: merge gene->arm map (from genemap, exported by data_prep)
    arm_map_path = cfg.DIR_TAB / "gene_arm_map.csv"
    if arm_map_path.exists() and "arm" not in ft.columns:
        am = pd.read_csv(arm_map_path)[["gene", "arm"]].drop_duplicates("gene")
        ft = ft.merge(am, on="gene", how="left")
    arm_col = next((c for c in ["arm", "chr_arm", "cytoband_arm"] if c in ft.columns), None)
    rho_arm = None
    if arm_col:
        oof_arm = leave_arm_out(ft, feats, arm_col)
        m = ~np.isnan(oof_arm)
        rho_arm = spearmanr(ft.transmissibility[m], oof_arm[m]).correlation

    # transfer: cross-lineage concordance from the committed transfer table
    transfer = ish.load_transfer()
    pct_cols = [c for c in transfer.columns if c.startswith("pct_")]
    W = None
    if len(pct_cols) >= 2:
        R = transfer[pct_cols].rank()
        # Kendall's W across lineages
        n, k = R.shape
        S = ((R.sum(axis=1) - R.sum(axis=1).mean())**2).sum()
        W = 12*S / (k**2 * (n**3 - n))

    out = ft[["gene", "transmissibility", "predicted_oof"]].copy()
    out.to_csv(cfg.DIR_TAB / "predictor_oof.csv", index=False)

    # ---- Figure 3: predictor + gene-intrinsic evidence ---------------------
    fig = plt.figure(figsize=(12, 4.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1, 1], wspace=0.34)
    C = "#2c6fbb"

    # panel a: predicted vs observed transmissibility (leave-gene-out)
    axa = fig.add_subplot(gs[0, 0])
    axa.scatter(oof, ft.transmissibility, s=5, alpha=0.25, color="#555",
                rasterized=True, edgecolors="none")
    lim = [min(oof.min(), ft.transmissibility.min()) - 0.02,
           max(oof.max(), ft.transmissibility.max()) + 0.02]
    axa.plot(lim, lim, "--", color="#999", lw=1)
    axa.set_xlim(lim); axa.set_ylim(lim)
    axa.set_xlabel("predicted transmissibility  (leave-gene-out)")
    axa.set_ylabel("observed transmissibility")
    axa.set_title(f"Predictable from gene properties alone\nρ = {rho:.2f},  R² = {r2:.2f}  "
                  f"(n = {len(ft):,} genes)", fontsize=9, loc="left")

    # panel b: positional control (leave-gene-out vs leave-arm-out)
    axb = fig.add_subplot(gs[0, 1])
    if rho_arm is not None:
        bars = axb.bar(["leave-\ngene-out", "leave-\narm-out"], [rho, rho_arm],
                       color=[C, "#7aa8d8"], width=0.6)
        axb.set_ylim(0, max(rho, rho_arm)*1.25)
        axb.set_ylabel("rank accuracy (ρ)")
        axb.set_title(f"Not positional:\nΔ = {abs(rho-rho_arm):.3f} holding out whole arms",
                      fontsize=9, loc="left")
        for b, v in zip(bars, [rho, rho_arm]):
            axb.text(b.get_x()+b.get_width()/2, v+0.01, f"{v:.2f}", ha="center", fontsize=8)
    else:
        axb.text(0.5, 0.5, "arm annotation\nnot in feature table\n(positional control run in\ndata_download refit)",
                 ha="center", va="center", fontsize=8, color="#777", transform=axb.transAxes)
        axb.axis("off")

    # panel c: transfer across lineages
    axc = fig.add_subplot(gs[0, 2])
    if W is not None:
        # pairwise rank concordance heat-ish summary: show per-lineage rank corr to mean
        meanrank = transfer[pct_cols].rank().mean(axis=1)
        corrs = [spearmanr(transfer[c], meanrank).correlation for c in pct_cols]
        labels = [c.replace("pct_", "") for c in pct_cols]
        axc.bar(range(len(corrs)), corrs, color="#6c5b9c", width=0.7)
        axc.set_xticks(range(len(corrs))); axc.set_xticklabels(labels, rotation=45, ha="right")
        axc.set_ylim(0, 1)
        axc.set_ylabel("ρ to consensus ranking")
        axc.set_title(f"Transfers across lineages\nKendall W = {W:.2f}", fontsize=9, loc="left")
    else:
        axc.axis("off")
    for ax, L in zip([axa, axb, axc], "abc"):
        ish.panel_letter(ax, L)

    fig.savefig(cfg.DIR_FIG / "fig3_predictor.png", dpi=200, bbox_inches="tight")

    ish.record("02_predictor", {
        "n_genes": int(len(ft)), "n_features": len(feats),
        "leave_gene_out_rho": float(rho), "leave_gene_out_R2": float(r2),
        "leave_arm_out_rho": None if rho_arm is None else float(rho_arm),
        "positional_delta": None if rho_arm is None else float(abs(rho - rho_arm)),
        "kendall_W": None if W is None else float(W),
    })
    print(f"leave-gene-out: rho={rho:.4f} R2={r2:.4f} n={len(ft)} feats={len(feats)}")
    print(f"positional: leave-arm-out rho={rho_arm} (delta from LGO)")
    print(f"transfer: Kendall W={W}")

if __name__ == "__main__":
    main()
