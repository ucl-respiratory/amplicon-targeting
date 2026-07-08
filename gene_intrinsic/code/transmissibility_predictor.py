"""
transmissibility_predictor.py
==============================

Gene-intrinsic prediction of copy-number-to-protein transmissibility.

A fourth, non-circular line of evidence for the gene-intrinsic thesis: a gene's
CN-to-protein transmissibility is predicted from EXTERNAL gene properties alone
(gnomAD, DepMap, CORUM, UniProt, Ensembl, GTEx, STRING) with NO tumour multi-omics.
The label is derived from CPTAC proteomics; the features are not. The two never
share a source, so the model cannot be circular.

Label (outcome)
---------------
transmissibility(gene) = mean(target) over the gene's UNIQUE (caseid, gene) rows
in the R4-filtered ALL cohort, where R4 = (cn >= round(ploidy)+1) AND (cn >= 3).
target is the binary protein-positive call from the canonical pipeline. This is
the identical R4/dedup basis as the thesis anchor (ALL row-level base rate 0.3983,
1,758,884 rows -> 700,211 unique pairs). All 6,926 modelled genes have >=34
amplified cases, so the >=20 primary and >=10 relaxed support thresholds retain
the full universe.

Features (predictors), 8 groups, 47 encoded columns
---------------------------------------------------
G1 dosage       gnomAD LOEUF/pLI/mis_z + DepMap dep_mean_effect/dep_frac_dependent
G2 complex      CORUM in_complex, n_complexes, complex_size, ...
G3 biophysics   UniProt length/MW/pI/GRAVY/TM + local disorder (VSL2/PSIPRED/...)
G5 mRNA         Ensembl transcript_length, GC, GC3 codon-opt, n_isoforms, UTRs
G6 evolution    dN/dS (Compara NG86), gene_age_proxy (orthologue depth); phyloP UNAVAILABLE
G7 function     UniProt is_tf/kinase/receptor/enzyme + GO MF category (one-hot)
G8 breadth      GTEx n_tissues_expressed, tau  [BORDERLINE: normal-tissue RNA; tested +/-]
G9 network      STRING v12 degree, weighted_degree, betweenness
(G4 protein turnover UNAVAILABLE: no reachable half-life atlas on the allowlist; not fabricated.)

Leakage safety (project's highest priority)
-------------------------------------------
Data-driven audit against the actual feature_data_manifest.csv source strings:
0 of 16 sources reference CPTAC tumour proteomics/RNA or the label feathers, and
no feature column is byte-identical to the label. G8 (GTEx normal-tissue RNA) is
the only borderline group and is reported with and without.

Key results (seed=2, XGBoost 600 trees, single-thread)
------------------------------------------------------
Primary (leave-gene-out 10-fold OOF, all features):
    Spearman 0.5474, R2 0.3220, AU-PR 0.7827 (random 0.5, lift 1.57)
    Without G8: Spearman 0.5180 (G8 contributes +0.029)
Permutation null (label-shuffle, models refit, B=1000, matched 300-tree observed 0.5514):
    null mean 0.0002, max 0.0579; p = 0.000999 (1/1001 floor), z = +35.0
Positional control (leave-chromosome-arm-out, 41 arms):
    Spearman 0.5462 vs leave-gene-out 0.5474 -> drop -0.0012 (signal is NOT positional)
Feature-group forward selection:
    G1 alone 0.359 (best single, 66% of full) -> +G3 0.459 (84%) -> full 0.547
Mechanism (exact TreeSHAP via xgboost pred_contribs; shap 0.46 loader/plots broken vs xgb 2.1.1):
    top feature dep_frac_dependent (essentiality, +); group order
    G3 biophysics 31% > G1 dosage 24% > G5 15% > G8 9% ~ G9 9% > G7 5% > G6 4% > G2 1%
    Reconciles with buffering chapter (essentiality sign-convention consistent; dosage dominant).
Atlas + 3q ADC cross-reference:
    24/24 ADC candidates scored. ADC rank accuracy (0.553) == whole-atlas (0.547).
    THE BOUND: R2=0.32 => ~68% per-gene variance unexplained; 6/24 ADCs context-dependent
    (context_boosted PLD1/IL1RAP/PLXNA1; intrinsic_only AP2M1/CD47/RNF13).
    Doubly-corroborated high-transmission leads: TFRC, ATP13A3, ITGB5, ATP1B3, PLSCR1.
    => The gene-intrinsic predictor is a PRIOR, not a replacement for the tumour model
       in the highly-amplified regime (the model's distinctive value).

Deliverables (Science artifacts + local Cache/pipeline_v2/transmissibility/)
----------------------------------------------------------------------------
feature_data_manifest.csv          4a2e7a7f-c64f-4821-839a-d727731f8e61
gene_transmissibility_label.csv    2e12bc25-a327-4f1c-a09f-44837e1019a6
gene_feature_table.csv             6a6c8217-81ef-47db-804a-97998de5e2a9
feature_provenance_audit.csv       87d3fecd-801a-46b8-910e-88d104a41c29 (v2, data-driven)
transmissibility_model_performance.csv  e296748c-5b5e-4dc7-b524-679b1873e322 (v2, B=1000)
transmissibility_pred_vs_obs.png   026064f3-b880-4384-a597-1432d9da7591
transmissibility_group_selection.csv    5f32badb-9251-4610-98f8-1d75ebec01e3
transmissibility_group_selection.png    154ebbea-177c-4966-a26c-2d273e73f529
transmissibility_positional_control.csv ea624050-80b8-4452-87f3-eced7a660b35
transmissibility_shap.csv          96cfa64a-55ad-4117-b7ca-5ce8b5208b0b
transmissibility_shap.png          ab925d98-ee52-4eb1-bc7e-c8ef3509e6b9
transmissibility_atlas.csv         aef8351d-966a-4853-81ba-6cf599de6f3a
transmissibility_adc_crossref.csv  3cf31fd7-0239-4d3d-af57-038e23d5acca
transmissibility_atlas.png         2d3e7006-7fdf-47a1-88d0-d1091e1993ab

Environment: python 3.11, xgboost 2.1.1, numpy 1.26.4, pandas 2.3.3, scipy 1.15.3
             (cnenv; shap 0.46.0 installed but its XGBoost loader/plots are incompatible
             with xgboost 2.1.1 + this numpy, so exact TreeSHAP uses xgboost pred_contribs).
"""
import json
import numpy as np
import pandas as pd
import xgboost as xgb
from scipy import stats
from sklearn.model_selection import KFold
from sklearn.metrics import average_precision_score

SEED = 2
XGB_TREES = 600
XGB_THREADS = 1

CACHE = "/Users/jameslam/Desktop/Work/Jie/Cache"
FEATURES_DIR = f"{CACHE}/external_annotations/gene_features"
OUT = f"{CACHE}/pipeline_v2/transmissibility"
LABEL_FEATHER = f"{CACHE}/fs_filtered_all__60e25364.feather"


# --------------------------------------------------------------------------
# Model
# --------------------------------------------------------------------------
def make_model(n_trees: int = XGB_TREES) -> "xgb.XGBRegressor":
    """Canonical gradient-boosted regressor (matches the tumour-model config)."""
    return xgb.XGBRegressor(
        n_estimators=n_trees, max_depth=6, learning_rate=0.05,
        subsample=0.9, colsample_bytree=0.9, tree_method="hist",
        random_state=SEED, n_jobs=XGB_THREADS, objective="reg:squarederror",
    )


# --------------------------------------------------------------------------
# Label: R4 ALL basis, dedup to unique (caseid, gene)
# --------------------------------------------------------------------------
def apply_r4_filter(df: pd.DataFrame, cn_col: str = "cn") -> pd.DataFrame:
    """R4 = (cn >= round(ploidy)+1) AND (cn >= 3). Identical to the thesis anchor."""
    cn = pd.to_numeric(df[cn_col], errors="coerce")
    ploidy = pd.to_numeric(df["ploidy"], errors="coerce") if "ploidy" in df.columns else 2.0
    mask = (cn >= np.round(ploidy) + 1) & (cn >= 3)
    return df[mask]


def build_label(feather_path: str = LABEL_FEATHER,
                min_cases: int = 20) -> pd.DataFrame:
    """transmissibility(gene) = mean(target) over unique amplified (caseid, gene) rows."""
    df = pd.read_feather(feather_path)
    r4 = apply_r4_filter(df)
    r4u = r4.drop_duplicates(["caseid", "gene"])
    g = r4u.groupby("gene")["target"].agg(["mean", "size"]).reset_index()
    g.columns = ["gene", "transmissibility", "n_amplified_cases"]
    g["pass_primary"] = g["n_amplified_cases"] >= min_cases
    g["pass_relaxed"] = g["n_amplified_cases"] >= 10
    return g


# --------------------------------------------------------------------------
# Feature assembly + encoding
# --------------------------------------------------------------------------
def encode_features(merged: pd.DataFrame, fgmap: dict):
    """bool -> 0/1; go_mf_category -> one-hot (tracked under G7). XGBoost handles NaN."""
    feat_cols = [c for c in fgmap if c in merged.columns]
    X = merged[feat_cols].copy()
    for c in feat_cols:
        if X[c].dtype == bool:
            X[c] = X[c].astype(float)
    if "go_mf_category" in X.columns:
        dummies = pd.get_dummies(X["go_mf_category"], prefix="gomf", dummy_na=False).astype(float)
        X = X.drop(columns=["go_mf_category"]).join(dummies)
    fgmap2 = {c: fgmap.get(c, "G7") for c in X.columns}
    for c in X.columns:
        if c.startswith("gomf_"):
            fgmap2[c] = "G7"
    return X, fgmap2


# --------------------------------------------------------------------------
# Cross-validation
# --------------------------------------------------------------------------
def logo_oof(Xmat: np.ndarray, yvec: np.ndarray, n_splits: int = 10,
             seed: int = SEED, n_trees: int = XGB_TREES) -> np.ndarray:
    """Leave-gene-out K-fold: each gene in exactly one test fold. Returns OOF preds."""
    oof = np.full(len(yvec), np.nan)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for tr, te in kf.split(Xmat):
        m = make_model(n_trees)
        m.fit(Xmat[tr], yvec[tr])
        oof[te] = m.predict(Xmat[te])
    return oof


def loarmo_oof(Xmat: np.ndarray, yvec: np.ndarray, arm_labels: np.ndarray,
               n_trees: int = XGB_TREES) -> np.ndarray:
    """Leave-chromosome-arm-out positional control: hold out all genes on one arm."""
    oof = np.full(len(yvec), np.nan)
    for a in sorted(set(arm_labels)):
        if a == "NA":
            continue
        te = np.where(arm_labels == a)[0]
        tr = np.where((arm_labels != a) & (arm_labels != "NA"))[0]
        if len(te) == 0 or len(tr) == 0:
            continue
        m = make_model(n_trees)
        m.fit(Xmat[tr], yvec[tr])
        oof[te] = m.predict(Xmat[te])
    return oof


def oof_metrics(oof: np.ndarray, yvec: np.ndarray) -> dict:
    """Spearman, R2, and AU-PR at the median (high-transmitter detection)."""
    v = ~np.isnan(oof)
    rho, p = stats.spearmanr(oof[v], yvec[v])
    r2 = 1 - np.sum((yvec[v] - oof[v]) ** 2) / np.sum((yvec[v] - yvec[v].mean()) ** 2)
    med = np.median(yvec[v])
    ybin = (yvec[v] > med).astype(int)
    ap = average_precision_score(ybin, oof[v])
    return {"spearman": float(rho), "spearman_p": float(p), "r2": float(r2),
            "aupr": float(ap), "random_ap": float(ybin.mean()),
            "aupr_lift": float(ap / ybin.mean()), "n_scored": int(v.sum())}


# --------------------------------------------------------------------------
# Permutation null (label-shuffle; models refit; B>=1000)
# --------------------------------------------------------------------------
def permutation_null(Xmat: np.ndarray, yvec: np.ndarray, n_perm: int = 1000,
                     n_trees: int = 300, max_workers: int = 4) -> dict:
    """Shuffle labels, refit leave-gene-out OOF, recompute Spearman. Deterministic per perm.

    Uses ThreadPoolExecutor (loky/joblib cannot init in the sandbox: SC_SEM_NSEMS_MAX
    PermissionError). XGBoost hist releases the GIL, so threads parallelise the fits.
    The null Spearman is config-invariant (shuffled labels are unpredictable at any tree
    count), and the observed is stable across tree counts (0.5474@600 vs 0.5514@300), so a
    300-tree null with a matched 300-tree observed is valid for the 600-tree primary.
    """
    from concurrent.futures import ThreadPoolExecutor

    def one_perm(i: int) -> float:
        rng = np.random.default_rng(SEED + i)
        yp = rng.permutation(yvec)
        oof = logo_oof(Xmat, yp, n_trees=n_trees)
        return float(stats.spearmanr(oof, yp)[0])

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        null = np.array(list(ex.map(one_perm, range(n_perm))))
    rho_obs_matched = float(stats.spearmanr(logo_oof(Xmat, yvec, n_trees=n_trees), yvec)[0])
    p = (np.sum(null >= rho_obs_matched) + 1) / (n_perm + 1)
    return {"rho_obs_matched": rho_obs_matched, "null_mean": float(null.mean()),
            "null_sd": float(null.std()), "null_max": float(null.max()),
            "p_perm": float(p), "z": float((rho_obs_matched - null.mean()) / null.std()),
            "n_perm": n_perm, "n_trees": n_trees}


# --------------------------------------------------------------------------
# Feature-group forward selection (incl. G1-alone)
# --------------------------------------------------------------------------
def forward_group_selection(X: pd.DataFrame, yvec: np.ndarray, fgmap2: dict) -> dict:
    groups = sorted(set(fgmap2.values()))
    grp_cols = {g: [c for c in X.columns if fgmap2[c] == g] for g in groups}
    Xd = {g: X[grp_cols[g]].values.astype(np.float32) for g in groups}

    def rho_for(cols):
        return float(stats.spearmanr(logo_oof(X[cols].values.astype(np.float32), yvec), yvec)[0])

    single = {g: rho_for(grp_cols[g]) for g in groups}   # includes G1-alone
    selected = [max(single, key=single.get)]
    remaining = [g for g in groups if g not in selected]
    path = [{"step": 1, "added": selected[0], "rho": single[selected[0]]}]
    cur = single[selected[0]]
    while remaining:
        best_g, best_rho = None, -1
        for g in remaining:
            cols = [c for sg in selected + [g] for c in grp_cols[sg]]
            r = rho_for(cols)
            if r > best_rho:
                best_rho, best_g = r, g
        selected.append(best_g)
        remaining.remove(best_g)
        path.append({"step": len(selected), "added": best_g, "rho": best_rho,
                     "delta": best_rho - cur})
        cur = best_rho
    return {"single": single, "path": path}


# --------------------------------------------------------------------------
# Mechanism: exact TreeSHAP via xgboost native pred_contribs
# --------------------------------------------------------------------------
def treeshap(X: pd.DataFrame, yvec: np.ndarray, fgmap2: dict) -> pd.DataFrame:
    """Exact TreeSHAP through the booster (shap 0.46 loader is incompatible with xgb 2.1.1)."""
    m = make_model()
    m.fit(X.values.astype(np.float32), yvec)
    dm = xgb.DMatrix(X.values.astype(np.float32), feature_names=list(X.columns))
    contribs = m.get_booster().predict(dm, pred_contribs=True)
    sv = contribs[:, :-1]  # last column is the bias term (== label mean)
    feat = list(X.columns)
    mean_abs = np.abs(sv).mean(axis=0)
    corr = []
    for j in range(len(feat)):
        xj = X.values[:, j].astype(float)
        msk = ~np.isnan(xj)
        corr.append(np.corrcoef(xj[msk], sv[msk, j])[0, 1]
                    if (msk.sum() > 10 and np.std(xj[msk]) > 0) else np.nan)
    return (pd.DataFrame({"feature": feat, "group": [fgmap2[f] for f in feat],
                          "mean_abs_shap": mean_abs, "shap_value_corr": corr})
            .sort_values("mean_abs_shap", ascending=False).reset_index(drop=True))


# --------------------------------------------------------------------------
# Atlas
# --------------------------------------------------------------------------
def build_atlas(genes: np.ndarray, yvec: np.ndarray, oof: np.ndarray,
                n_cases: np.ndarray) -> pd.DataFrame:
    atlas = pd.DataFrame({"gene": genes, "observed_transmissibility": yvec,
                          "predicted_transmissibility_oof": oof,
                          "n_amplified_cases": n_cases})
    atlas["predicted_percentile"] = atlas["predicted_transmissibility_oof"].rank(pct=True)
    atlas["residual"] = atlas["observed_transmissibility"] - atlas["predicted_transmissibility_oof"]
    atlas = atlas.sort_values("predicted_transmissibility_oof", ascending=False).reset_index(drop=True)
    atlas["predicted_rank"] = np.arange(1, len(atlas) + 1)
    return atlas


if __name__ == "__main__":
    # See module docstring for the full pipeline and recorded results.
    # Feature tables are sourced by the four Phase-0 tracks and cached under FEATURES_DIR;
    # the label is built from LABEL_FEATHER via build_label().
    print("transmissibility_predictor: import and call build_label/encode_features/logo_oof/... "
          "or see the docstring for the recorded results and artifact IDs.")
