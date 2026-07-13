#!/usr/bin/env python3
# =============================================================================
# 05_single_cell_andgate.py  --  Why AND-gated multi-antigen targeting is
# selective: same-cell co-detection in tumour, threshold collapse in normal.
# -----------------------------------------------------------------------------
# Two single-cell arguments, one figure:
#
#  (A) On-tumour co-detection. In malignant cells, are all antigens of a co-target
#      set co-detected in the SAME cell more often than if detections were
#      independent?  enrichment = P(all detected) / prod(P(detected)_gene.
#      Donor-block bootstrap CI; marginal-preserving permutation p.  An AND-gate
#      only fires where antigens co-occur, so >1 enrichment is what makes a
#      multi-antigen construct efficient on tumour.
#
#  (B) Off-tumour selectivity. A normal cell binds an AND-gate only if its
#      LIMITING (weakest) antigen clears the binding threshold. Raising the
#      per-antigen threshold from detection (10 nTPM) to binding-relevant
#      (25/50 nTPM) collapses the normal-cell co-expression burden, provided each
#      set contains at least one antigen that is low in normal tissue.
#
# Inputs: CELLxGENE malignant-cell slices (stage 16) for (A); HPA single-cell
# normal for (B). Runs once data_download provides them.
# Outputs: figures/fig6_andgate.png, tables/andgate_enrichment.csv
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cnt_shared as ish
from cnt_shared import cfg
ish.apply_style(sizes=(9, 8, 7))

# Co-target sets testable in single-cell: sourced from the pan-cancer construct
# table (04d) so this figure always matches Table 2, restricted to the cohorts
# with a malignant single-cell slice (LUAD, LSCC measured; GBM prediction-only).
# Each set is (cohort, [antigens]).
SC_COHORTS = ("LUAD", "LSCC", "GBM")

def _load_sets():
    con = pd.read_csv(cfg.DIR_TAB / "adc_constructs.csv")
    con = con[con.single_cell_tested & con.cohort.isin(SC_COHORTS)]
    sets, order = {}, []
    # order: cohort group, then most-enriched first within cohort
    for _, r in con.sort_values(["cohort", "enrich"], ascending=[True, False]).iterrows():
        amp = r["amplicon"]
        sets[amp] = (r["cohort"], [g.strip() for g in str(r["antigens"]).split("+")])
        order.append(amp)
    return sets, order

SETS, SET_ORDER = _load_sets()

def _depth_bins(nnz, n_bins=cfg.SC_DEPTH_BINS):
    qs = np.quantile(nnz, np.linspace(0, 1, n_bins + 1)); qs[0] -= 1; qs[-1] += 1
    return np.digitize(nnz, qs[1:-1])

def _perm_codet_depth(D, binid, rng):
    """One depth-stratified permutation: within each per-cell depth bin, shuffle
    each antigen column independently — preserves per-gene marginals AND per-cell
    depth, so co-detection expected from depth alone is the null."""
    Dp = np.empty_like(D); k = D.shape[1]
    for b in np.unique(binid):
        ix = np.where(binid == b)[0]
        if len(ix) == 0: continue
        for j in range(k):
            Dp[ix, j] = D[ix[rng.permutation(len(ix))], j]
    return Dp.all(1).mean()

def enrich_ci(D, donors, nnz, rng, n_boot=cfg.N_BOOTSTRAP, n_perm=cfg.N_PERMUTATION):
    """Same-cell co-detection enrichment against a DEPTH-STRATIFIED null.
    Co-detection is inflated by per-cell sequencing depth, so we permute each
    antigen column within per-cell depth (nnz) deciles rather than globally: the
    enrichment is observed / depth-expected co-detection, and the permutation p
    asks whether observed exceeds what depth structure alone produces. Also
    reports the naive marginal-null enrichment for comparison."""
    n = D.shape[0]
    obs = D.all(1).mean()
    exp_marg = float(np.prod(D.mean(0)))
    binid = _depth_bins(nnz)
    perm = np.array([_perm_codet_depth(D, binid, rng) for _ in range(n_perm)])
    exp_depth = float(perm.mean())
    e0 = obs / exp_depth if exp_depth > 0 else np.nan
    pval = (np.sum(perm >= obs) + 1) / (n_perm + 1)
    uniq = np.unique(donors); d_to_rows = {d: np.where(donors == d)[0] for d in uniq}
    boots = []
    for _ in range(n_boot):
        samp = rng.choice(uniq, size=len(uniq), replace=True)
        rows = np.concatenate([d_to_rows[d] for d in samp])
        Db = D[rows]; nb = nnz[rows]
        eb = _perm_codet_depth(Db, _depth_bins(nb), rng)
        if eb > 0: boots.append(Db.all(1).mean() / eb)
    lo, hi = np.percentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan)
    return dict(n=int(n), n_donors=int(len(uniq)), obs_pct=obs*100,
                exp_pct=exp_depth*100, exp_marg_pct=exp_marg*100,
                enrich=e0, enrich_marginal=obs/exp_marg if exp_marg>0 else np.nan,
                ci_lo=lo, ci_hi=hi, perm_p=pval)

def load_cxg(code):
    p = cfg.PATHS[f"cxg_{code.lower()}"]
    return pd.read_parquet(p) if p.exists() else None

def normal_burden():
    """HPA single-cell normal: fraction of normal cell types where the LIMITING
    antigen of each set exceeds threshold T (detection 10 vs binding 25/50)."""
    p = cfg.PATHS["hpa_singlecell"]
    if not p.exists(): return None
    sc = pd.read_csv(p, sep="\t")
    gcol = "Gene name" if "Gene name" in sc.columns else sc.columns[1]
    ccol = "Cell type" if "Cell type" in sc.columns else sc.columns[2]
    scw = sc.pivot_table(index=gcol, columns=ccol, values="nTPM", aggfunc="max")
    rows = []
    for sn in SET_ORDER:
        genes = [g for g in SETS[sn][1] if g in scw.index]
        if not genes: continue
        for T in (cfg.SC_DETECT_NTPM, cfg.SC_BINDING_NTPM, 50):
            limiting = (scw.loc[genes] >= T).all(0)   # all antigens clear T in a cell type
            rows.append({"set": sn, "threshold": T,
                         "frac_normal_celltypes": float(limiting.mean())})
    return pd.DataFrame(rows)

def main():
    have_cxg = all(load_cxg(c) is not None for c in SC_COHORTS)
    rng = np.random.default_rng(cfg.SEED)
    ci = None
    if have_cxg:
        slices = {c: load_cxg(c) for c in SC_COHORTS}
        # Enrichment/CI/perm_p are the SINGLE authoritative construct values computed
        # by 04d_constructs.py (adc_constructs.csv); this figure reads them so Table 2
        # and Figure 6a are identical by construction rather than two independent
        # bootstraps of the same quantity. We recompute only the ABSOLUTE co-detection
        # fraction (obs_pct) here, which 04d does not store and the text quotes.
        con = pd.read_csv(cfg.DIR_TAB / "adc_constructs.csv")
        con = con[con.single_cell_tested & con.amplicon.isin(SET_ORDER)]
        rows = []
        for sn in SET_ORDER:
            r = con[con.amplicon == sn]
            if not len(r):
                continue
            r = r.iloc[0]
            code, genes = SETS[sn]; df = slices[code]
            genes = [g for g in genes if g in df.columns]
            D = (df[genes].to_numpy() > 0).astype(int) if len(genes) >= 2 else None
            obs_pct = float(D.all(1).mean() * 100) if D is not None else np.nan
            rows.append(dict(set=sn, enrich=r["enrich"], ci_lo=r["ci_lo"],
                             ci_hi=r["ci_hi"],
                             ci_lo_raw=r.get("ci_lo_raw", r["ci_lo"]),
                             perm_p=r["perm_p"],
                             enrich_marginal=r.get("enrich_marginal", np.nan),
                             n=int(r["n_cells"]), n_donors=int(r["n_donors"]),
                             obs_pct=obs_pct))
        ci = pd.DataFrame(rows)
        ci.to_csv(cfg.DIR_TAB / "andgate_enrichment.csv", index=False)
    else:
        print("CELLxGENE slices absent — run data_download stage 16 "
              "(cnt-census env). AND-gate figure deferred.")

    nb = normal_burden()

    # ---- Figure 6: (a) on-tumour co-detection, (b) normal-tissue collapse ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
    # cohort colours threaded with the manuscript's other figures
    C_COH = {"LUAD": "#2c6fbb", "LSCC": "#c0392b", "GBM": "#3a3a3a"}

    axa = axes[0]
    if ci is not None and len(ci):
        ciI = ci.set_index("set").reindex([s for s in SET_ORDER if s in ci.set.values])
        y = np.arange(len(ciI))
        err = np.vstack([ciI.enrich - ciI.ci_lo, ciI.ci_hi - ciI.enrich])
        coh = [SETS[s][0] for s in ciI.index]
        cols = [C_COH.get(c, "#777") for c in coh]
        # significance: permutation p<0.05 AND bootstrap CI clear of 1.0
        # (use full-precision lower bound so 2-dp rounding doesn't flip a call)
        _lo = ciI["ci_lo_raw"] if "ci_lo_raw" in ciI.columns else ciI["ci_lo"]
        sig = (ciI.perm_p < 0.05) & (_lo > 1.0)
        bars = axa.barh(y, ciI.enrich, xerr=err, color=cols,
                        error_kw=dict(lw=1.0, capsize=2.5))
        # non-significant bars: hollow (face lightened) so they read as "no enrichment"
        for b, s in zip(bars, sig):
            if not s:
                b.set_alpha(0.35); b.set_hatch("///")
        axa.axvline(1, color="k", ls="--", lw=0.9)
        axa.set_yticks(y); axa.set_yticklabels(ciI.index); axa.invert_yaxis()
        axa.set_xlabel("same-cell co-detection enrichment\n(observed / depth-expected)")
        axa.set_title("On tumour: antigens co-detected in the same malignant\n"
                      "cell above a depth-matched null (hatched = not significant)",
                      fontsize=9, loc="left")
        axa.annotate("independent", xy=(1, len(ciI)-0.5), xytext=(4, 0),
                     textcoords="offset points", fontsize=7, color="#555")
        # cohort legend (threaded colours)
        from matplotlib.patches import Patch
        seen = [c for c in ("LUAD", "LSCC", "GBM") if c in coh]
        lab = {"LUAD": "LUAD (measured)", "LSCC": "LSCC (measured)",
               "GBM": "GBM (prediction-only)"}
        axa.legend(handles=[Patch(facecolor=C_COH[c], label=lab[c]) for c in seen],
                   frameon=False, fontsize=7, loc="lower right")
    else:
        axa.text(0.5, 0.5, "CELLxGENE malignant-cell slices\nnot yet built\n"
                 "(data_download stage 16)", ha="center", va="center",
                 fontsize=9, color="#777", transform=axa.transAxes)
        axa.axis("off")

    axb = axes[1]
    if nb is not None and len(nb):
        piv = nb.pivot(index="set", columns="threshold", values="frac_normal_celltypes")
        piv = piv.reindex([s for s in SET_ORDER if s in piv.index])
        x = np.arange(len(piv)); w = 0.25
        for i, T in enumerate(sorted(piv.columns)):
            axb.bar(x + (i-1)*w, piv[T]*100, w,
                    label=f"{'detection' if T==10 else 'binding'} ≥{T} nTPM")
        axb.set_xticks(x); axb.set_xticklabels(piv.index, rotation=45, ha="right")
        axb.set_ylabel("normal cell types with all antigens ≥ T  (%)")
        axb.set_title("Off tumour: raising to binding threshold\ncollapses the normal burden",
                      fontsize=9, loc="left")
        axb.legend(frameon=False, fontsize=7)
    else:
        axb.text(0.5, 0.5, "HPA single-cell normal\nnot yet built\n(data_download stage 04)",
                 ha="center", va="center", fontsize=9, color="#777", transform=axb.transAxes)
        axb.axis("off")
    for ax, L in zip(axes, "ab"):
        if ax.axison: ish.panel_letter(ax, L)

    fig.savefig(cfg.DIR_FIG / "fig6_andgate.png", dpi=200, bbox_inches="tight")

    ish.record("05_single_cell_andgate", {
        "have_cellxgene": bool(have_cxg),
        "enrichment": ci.to_dict("records") if ci is not None else None,
        "normal_burden": nb.to_dict("records") if nb is not None else None,
    })
    if ci is not None:
        print(ci[["set", "n", "n_donors", "enrich", "ci_lo", "ci_hi", "perm_p"]].to_string(index=False))

if __name__ == "__main__":
    main()
