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
# with a malignant single-cell slice (LUAD, LSCC). Each set is (cohort, [antigens]).
def _load_sets():
    con = pd.read_csv(cfg.DIR_TAB / "adc_constructs.csv")
    con = con[con.single_cell_tested & con.cohort.isin(["LUAD", "LSCC"])]
    sets, order = {}, []
    for _, r in con.sort_values(["cohort", "enrich"], ascending=[True, False]).iterrows():
        amp = r["amplicon"]
        sets[amp] = (r["cohort"], [g.strip() for g in str(r["antigens"]).split("+")])
        order.append(amp)
    return sets, order

SETS, SET_ORDER = _load_sets()

def enrich_ci(D, donors, rng, n_boot=cfg.N_BOOTSTRAP, n_perm=cfg.N_PERMUTATION):
    """Same-cell co-detection enrichment with donor-block bootstrap CI and a
    marginal-preserving permutation p (shuffle each gene's column independently)."""
    n = D.shape[0]
    obs = D.all(1).mean()
    exp = float(np.prod(D.mean(0)))
    e0 = obs / exp if exp > 0 else np.nan
    uniq = np.unique(donors); d_to_rows = {d: np.where(donors == d)[0] for d in uniq}
    boots = []
    for _ in range(n_boot):
        samp = rng.choice(uniq, size=len(uniq), replace=True)
        Db = D[np.concatenate([d_to_rows[d] for d in samp])]
        ee = float(np.prod(Db.mean(0)))
        if ee > 0: boots.append(Db.all(1).mean() / ee)
    lo, hi = np.percentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan)
    perm = np.empty(n_perm)
    for i in range(n_perm):
        Dp = np.empty_like(D)
        for j in range(D.shape[1]): Dp[:, j] = D[rng.permutation(n), j]
        perm[i] = Dp.all(1).mean()
    pval = (np.sum(perm >= obs) + 1) / (n_perm + 1)
    return dict(n=int(n), n_donors=int(len(uniq)), obs_pct=obs*100, exp_pct=exp*100,
                enrich=e0, ci_lo=lo, ci_hi=hi, perm_p=pval)

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
    have_cxg = all(load_cxg(c) is not None for c in ("LUAD", "LSCC"))
    rng = np.random.default_rng(cfg.SEED)
    ci = None
    if have_cxg:
        slices = {c: load_cxg(c) for c in ("LUAD", "LSCC")}
        rows = []
        for sn in SET_ORDER:
            code, genes = SETS[sn]; df = slices[code]
            genes = [g for g in genes if g in df.columns]
            if len(genes) < 2: continue
            D = df[genes].to_numpy().astype(int); don = df["donor_id"].to_numpy()
            r = enrich_ci(D, don, rng); r["set"] = sn; rows.append(r)
        ci = pd.DataFrame(rows)
        ci.to_csv(cfg.DIR_TAB / "andgate_enrichment.csv", index=False)
    else:
        print("CELLxGENE slices absent — run data_download stage 16 "
              "(cnt-census env). AND-gate figure deferred.")

    nb = normal_burden()

    # ---- Figure 6: (a) on-tumour co-detection, (b) normal-tissue collapse ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    C_LUAD, C_LSCC = "#2c6fbb", "#c0392b"

    axa = axes[0]
    if ci is not None and len(ci):
        ciI = ci.set_index("set").reindex([s for s in SET_ORDER if s in ci.set.values])
        y = np.arange(len(ciI))
        err = np.vstack([ciI.enrich - ciI.ci_lo, ciI.ci_hi - ciI.enrich])
        cols = [C_LUAD if s.startswith("LUAD") else C_LSCC for s in ciI.index]
        axa.barh(y, ciI.enrich, xerr=err, color=cols, error_kw=dict(lw=1.1, capsize=3))
        axa.axvline(1, color="k", ls="--", lw=0.9)
        axa.set_yticks(y); axa.set_yticklabels(ciI.index); axa.invert_yaxis()
        axa.set_xlabel("same-cell co-detection enrichment\n(observed / independent)")
        axa.set_title("On tumour: antigens co-occur in the same malignant cell\n"
                      "far above independence", fontsize=9, loc="left")
        axa.annotate("independent", xy=(1, len(ciI)-0.5), xytext=(4, 0),
                     textcoords="offset points", fontsize=7, color="#555")
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
