# =============================================================================
# 31_tumour_samecell.py  --  Fig 6: same-cell co-expression in MALIGNANT cells.
# Reads the data_download stage-16 CELLxGENE detection slices (LUAD, LSCC) and,
# for each co-target set, tests whether all antigens are co-detected in the SAME
# malignant cell more often than expected if detections were independent:
#   enrichment = observed co-detection / product of per-gene detection rates.
# CIs from a donor-block bootstrap (resample donors); significance from a
# marginal-preserving permutation (shuffle each gene's column independently).
# Writes tumour_singlecell_enrichment_CI.csv. Skipped if slices are absent.
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cnt_io
from cnt_io import DIR_TAB, PATHS, savefig, record_values, cfg

# co-target sets testable in single-cell (LUAD + LSCC have malignant cells)
SETS = {"LUAD_1q":("LUAD",["ADAM15","CD46","EFNA1","MUC1","NCSTN","XPR1"]),
        "LUAD_7p":("LUAD",["DAGLB","EGFR","ITGB8","TSPAN13"]),
        "LUAD_5p":("LUAD",["CLPTM1L","SLC12A7"]),
        "LSCC_1q":("LSCC",["F11R","HSD17B7","NCSTN"]),
        "LSCC_20q":("LSCC",["GGT7","SDC4","TM9SF4"])}
SET_ORDER = ["LUAD_1q","LUAD_7p","LUAD_5p","LSCC_1q","LSCC_20q"]

def enrich_ci(D, donors, rng, n_boot=cfg.N_BOOTSTRAP, n_perm=cfg.N_PERMUTATION):
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
    lo, hi = np.percentile(boots, [2.5, 97.5])
    perm = np.empty(n_perm)
    for i in range(n_perm):
        Dp = np.empty_like(D)
        for j in range(D.shape[1]): Dp[:, j] = D[rng.permutation(n), j]
        perm[i] = Dp.all(1).mean()
    pval = (np.sum(perm >= obs) + 1) / (n_perm + 1)
    return dict(n=int(n), n_donors=int(len(uniq)), obs_pct=round(obs*100,3), exp_pct=round(exp*100,3),
                enrich=round(e0,2), ci_lo=round(lo,2), ci_hi=round(hi,2), perm_p=pval)

def main():
    if not (PATHS["cxg_luad"].exists() and PATHS["cxg_lscc"].exists()):
        sys.stderr.write("[31] CELLxGENE slices absent; run data_download stage 16. Skipping Fig 6.\n")
        return 0
    rng = np.random.default_rng(cfg.SEED_BOOTSTRAP)
    slices = {c: cnt_io.load_cellxgene(c) for c in ("LUAD","LSCC")}
    rows = []
    for sn in SET_ORDER:
        code, genes = SETS[sn]; df = slices[code]
        D = df[genes].to_numpy().astype(int); don = df["donor_id"].to_numpy()
        r = enrich_ci(D, don, rng); r["set"] = sn; rows.append(r)
    ci = pd.DataFrame(rows)[["set","n","n_donors","obs_pct","exp_pct","enrich","ci_lo","ci_hi","perm_p"]]
    ci.to_csv(DIR_TAB / "tumour_singlecell_enrichment_CI.csv", index=False)

    plt.rcParams.update({"font.size":8,"axes.titlesize":9,"axes.labelsize":8,"legend.fontsize":7,
        "xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,
        "axes.titlelocation":"left"})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    ciI = ci.set_index("set").loc[SET_ORDER]
    ax = axes[0]; y = np.arange(len(SET_ORDER))
    err = np.vstack([ciI.enrich - ciI.ci_lo, ciI.ci_hi - ciI.enrich])
    cols = ["#4c72b0" if s.startswith("LUAD") else "#c44e52" for s in SET_ORDER]
    ax.barh(y, ciI.enrich, xerr=err, color=cols, error_kw=dict(lw=1.1, capsize=3))
    ax.axvline(1, color="k", ls="--", lw=0.8, label="no enrichment (independent)")
    ax.set_yticks(y); ax.set_yticklabels(SET_ORDER); ax.invert_yaxis()
    ax.set_xlabel("Same-cell co-detection enrichment\n(observed / independent expectation)")
    ax.set_title("A  All antigens co-detected in the same malignant cell\nmore than chance (95% donor-bootstrap CI)")
    ax.legend(frameon=False, loc="lower right", fontsize=6.5)
    for i, s in enumerate(SET_ORDER):
        ax.text(ciI.ci_hi.iloc[i]+0.3, i, f"{ciI.enrich.iloc[i]:.1f}\u00d7", va="center", fontsize=7)
    ax = axes[1]; w = 0.38
    ax.bar(y-w/2, ciI.obs_pct, w, label="observed same-cell", color="#2c6fb0")
    ax.bar(y+w/2, ciI.exp_pct, w, label="independent expectation", color="#b9c9de")
    ax.set_xticks(y); ax.set_xticklabels(SET_ORDER, rotation=30, ha="right"); ax.set_ylabel("Malignant cells co-expressing all (%)")
    ax.set_title("B  Observed vs independent co-detection rate"); ax.legend(frameon=False, loc="upper right")
    for i in range(len(SET_ORDER)):
        ax.text(y[i]-w/2, ciI.obs_pct.iloc[i]+0.1, f"{ciI.obs_pct.iloc[i]:.1f}", ha="center", fontsize=6)
    fig.suptitle("On-tumour same-cell co-expression of co-target sets in malignant cells (CELLxGENE census)",
                 fontsize=9.5, y=1.02)
    fig.tight_layout(); savefig(fig, "tumour_singlecell_coexpression.png"); plt.close(fig)

    record_values("31_tumour_samecell", {
        "sets": {r.set: {"n_cells":int(r.n),"n_donors":int(r.n_donors),"obs_pct":float(r.obs_pct),
                 "exp_pct":float(r.exp_pct),"enrich":float(r.enrich),"ci_lo":float(r.ci_lo),
                 "ci_hi":float(r.ci_hi),"perm_p":float(r.perm_p)} for r in ci.itertuples()},
    })
    print(ci.to_string(index=False))

if __name__ == "__main__":
    main()
