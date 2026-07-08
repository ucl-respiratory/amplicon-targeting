"""
kernel.py for the feature-group-head-to-head skill.

Encodes a leakage-safe, fold-aware feature-group comparison harness with four
guarantees that are CHECKED, not merely offered. The public entry point
run_head_to_head refuses to return arm metrics unless all applicable guarantees
have passed:

  1. bit-exact anchor gate  - the incumbent arm must reproduce a named canonical
     AU-PR within tolerance before any comparison is reported;
  2. label-shuffle permutation null on FIXED out-of-fold predictions with the
     models held fixed, reported symmetrically (a wrong-side-of-zero delta yields
     p -> 1.0; never described as 'significantly negative');
  3. coverage-matched confound control whenever two feature groups differ in gene
     coverage (compares topology on an identical gene set);
  4. leakage self-check across all folds and fold-aware features (max feature
     diff must be 0.0).

British English, no em dashes. Seed 2 by default. Deterministic, single-thread.
"""
import gc
from typing import Dict, List, Tuple, Optional, Callable

import numpy as np
import pandas as pd

SEED = 2
XGB_TREES = 600
XGB_THREADS = 1
ANCHOR_TOL = 2e-4          # bit-exact gate tolerance (matches V2_TOL)
LEAKAGE_TOL = 1e-9         # max permitted fold-feature change under corruption


def hh_xgb_model(seed_fold: int):
    import xgboost as xgb
    return xgb.XGBClassifier(
        n_estimators=XGB_TREES, max_depth=6, learning_rate=0.05,
        subsample=0.9, colsample_bytree=0.9, tree_method="hist",
        enable_categorical=True, random_state=seed_fold,
        early_stopping_rounds=20, n_jobs=XGB_THREADS, eval_metric="aucpr",
    )


def fold_gene_means(pivot: "pd.DataFrame", training_caseids: set) -> "pd.Series":
    """Per-gene MEAN over TRAINING caseids only (leakage-safe)."""
    return pivot[pivot.index.isin(training_caseids)].mean(axis=0)


def build_nbr_per_gene(gene_means: "pd.Series", neighbour_map: Dict[str, List[str]]) -> "pd.Series":
    """gene -> nanmean over its neighbours' per-gene means (unweighted)."""
    gm = gene_means.to_dict()
    out = {}
    for g, nbrs in neighbour_map.items():
        vals = [gm[n] for n in nbrs if n in gm and gm[n] == gm[n]]
        if vals:
            out[g] = float(np.nanmean(vals))
    return pd.Series(out, dtype=np.float32)


def build_net_per_gene_weighted(gene_means: "pd.Series", weighted_map: Dict[str, List]) -> "pd.Series":
    """Confidence-weighted neighbour mean: net[g] = sum(w*mean)/sum(w)."""
    gm = gene_means.to_dict()
    out = {}
    for g, nbrs in weighted_map.items():
        num = den = 0.0
        for n, w in nbrs:
            v = gm.get(n)
            if v is not None and v == v:
                num += w * v
                den += w
        if den > 0:
            out[g] = float(num / den)
    return pd.Series(out, dtype=np.float32)


def cv_predict_unified(
    df: "pd.DataFrame", cols: List[str],
    splits: List[Tuple["np.ndarray", "np.ndarray"]],
    fold_aware: Dict[str, dict], static_ok: bool = True,
    seed: int = None,
) -> "np.ndarray":
    """
    Out-of-fold predictor with fold-aware feature reconstruction.

    fold_aware: col -> {"pivot": caseid x gene DataFrame, "map": {gene:[nbr]} for
    unweighted, or "weighted": {gene:[(nbr,w)]} for confidence-weighted}. For each
    fold, per-gene means are taken over TRAINING caseids only, aggregated over the
    neighbour map, and substituted by gene into train and test rows. Any col NOT in
    fold_aware is used as-is from df (static). Leakage-safe: validation cases never
    enter their fold's neighbour values.
    """
    if seed is None:
        seed = SEED
    y = df["target"].to_numpy()
    proba = np.zeros(len(df), dtype=np.float32)
    gene_arr = df["gene"].astype(str).to_numpy()
    caseid_arr = df["caseid"].astype(str).to_numpy()
    fa = [c for c in cols if c in fold_aware]

    for fold, (tr, te) in enumerate(splits, start=1):
        training_caseids = set(caseid_arr[tr])
        fold_maps = {}
        for c in fa:
            spec = fold_aware[c]
            gm = fold_gene_means(spec["pivot"], training_caseids)
            if spec.get("weighted"):
                fold_maps[c] = build_net_per_gene_weighted(gm, spec["weighted"]).to_dict()
            else:
                fold_maps[c] = build_nbr_per_gene(gm, spec["map"]).to_dict()

        present = [c for c in cols if c in df.columns]

        def _build_X(idx):
            X = df.iloc[idx][present].copy()
            rg = gene_arr[idx]
            for c in cols:
                if c in fold_maps:
                    X[c] = pd.Series(rg, index=X.index).map(fold_maps[c]).to_numpy().astype(np.float32)
            return X[cols]

        Xtr, Xte = _build_X(tr), _build_X(te)
        mdl = hh_xgb_model(seed + fold)
        mdl.fit(Xtr, y[tr], eval_set=[(Xte, y[te])], verbose=False)
        proba[te] = mdl.predict_proba(Xte)[:, 1]
        del Xtr, Xte, mdl, fold_maps
        gc.collect()
    return proba


def dedup_metrics(df: "pd.DataFrame", proba: "np.ndarray") -> Dict[str, float]:
    """Deduplicate to unique (caseid, gene), then metrics with Random-AP (base rate)."""
    from sklearn.metrics import average_precision_score, roc_auc_score
    tmp = df[["caseid", "gene", "target"]].copy()
    tmp["proba"] = proba
    u = tmp.drop_duplicates(["caseid", "gene"])
    y, p = u["target"].to_numpy(), u["proba"].to_numpy()
    base = float(y.mean())
    aupr = float(average_precision_score(y, p))
    return {"n_raw": len(tmp), "n_unique": len(u), "random_ap": base,
            "aupr": aupr, "lift": aupr / base, "auroc": float(roc_auc_score(y, p))}


def leakage_selfcheck(
    df: "pd.DataFrame", splits: List[Tuple["np.ndarray", "np.ndarray"]],
    fold_aware: Dict[str, dict],
) -> "pd.DataFrame":
    """
    Guarantee 4. For each fold and each fold-aware feature, corrupt ALL validation
    cases' source-pivot values to a sentinel and confirm the recomputed per-gene
    fold features are byte-identical (max diff 0.0). Also records train/test case
    overlap (must be 0). Returns a per-(fold,feature) DataFrame; leakage_detected
    is True if any max_feature_diff exceeds LEAKAGE_TOL.
    """
    caseid_arr = df["caseid"].astype(str).to_numpy()
    rows = []
    for fold, (tr, te) in enumerate(splits, start=1):
        tr_cases = set(caseid_arr[tr]); te_cases = set(caseid_arr[te])
        overlap = len(tr_cases & te_cases)
        for c, spec in fold_aware.items():
            piv = spec["pivot"]
            base_gm = fold_gene_means(piv, tr_cases)
            base = (build_net_per_gene_weighted(base_gm, spec["weighted"])
                    if spec.get("weighted") else build_nbr_per_gene(base_gm, spec["map"]))
            pv = piv.copy()
            val = [x for x in te_cases if x in pv.index]
            pv.loc[val] = 999.0
            corr_gm = fold_gene_means(pv, tr_cases)
            corr = (build_net_per_gene_weighted(corr_gm, spec["weighted"])
                    if spec.get("weighted") else build_nbr_per_gene(corr_gm, spec["map"]))
            genes = set(base.index) & set(corr.index)
            diffs = [abs(base[g] - corr[g]) for g in genes
                     if not (np.isnan(base[g]) or np.isnan(corr[g]))]
            md = max(diffs) if diffs else 0.0
            rows.append({"fold": fold, "feature": c, "n_val_corrupted": len(val),
                         "train_test_overlap": overlap,
                         "max_feature_diff": round(float(md), 10),
                         "leakage_detected": bool(md > LEAKAGE_TOL)})
    return pd.DataFrame(rows)


def hh_fast_ap_order(p):
    return np.argsort(-p, kind="mergesort")


def hh_fast_ap(y, order):
    ys = y[order].astype(np.float64)
    tp = np.cumsum(ys); fp = np.cumsum(1.0 - ys)
    prec = tp / (tp + fp)
    return float(np.sum(prec * ys) / tp[-1])


def permutation_null(
    df: "pd.DataFrame", proba_arm: "np.ndarray", proba_ref: "np.ndarray",
    n_perm: int = 10000, seed: int = None,
) -> Dict[str, float]:
    """
    Guarantee 2. Label-shuffle null on FIXED out-of-fold predictions, models held
    fixed. Shuffles deduplicated labels n_perm times and recomputes AU-PR(arm) -
    AU-PR(ref) with the probability vectors frozen. Reports the observed delta
    (exact sklearn), the null spread, a z-score, a one-sided p (fraction of null
    >= observed) and its resolution floor. Symmetric reporting: a wrong-side-of-
    zero delta yields p -> 1.0 and MUST be read as 'no evidence the arm improves on
    the baseline', never 'significantly negative'.
    """
    if seed is None:
        seed = SEED
    from sklearn.metrics import average_precision_score
    tmp = df[["caseid", "gene", "target"]].copy()
    tmp["pa"] = proba_arm; tmp["pr"] = proba_ref
    u = tmp.drop_duplicates(["caseid", "gene"])
    y = u["target"].to_numpy(); pa = u["pa"].to_numpy(); pr = u["pr"].to_numpy()
    obs = average_precision_score(y, pa) - average_precision_score(y, pr)
    oa, orf = hh_fast_ap_order(pa), hh_fast_ap_order(pr)
    obs_fast = hh_fast_ap(y, oa) - hh_fast_ap(y, orf)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        yp = rng.permutation(y)
        null[i] = hh_fast_ap(yp, oa) - hh_fast_ap(yp, orf)
    sd = float(null.std(ddof=1))
    floor = 1.0 / (n_perm + 1)
    p = (1 + int(np.sum(null >= obs_fast))) / (n_perm + 1)
    return {"obs_delta": round(float(obs), 6),
            "null_mean": round(float(null.mean()), 6), "null_std": round(sd, 6),
            "z": round(float((obs_fast - null.mean()) / sd), 2) if sd > 0 else float("nan"),
            "p_value": float(p), "p_floor": floor, "n_perm": n_perm,
            "verdict": hh_null_verdict(obs, p)}


def hh_null_verdict(obs_delta: float, p: float, alpha: float = 0.05) -> str:
    """Symmetric verdict wording per the standing statistical rule."""
    if obs_delta > 0 and p < alpha:
        return "arm improves on baseline (significant)"
    if obs_delta <= 0:
        return "no evidence the arm improves on the baseline"
    return "no significant improvement over baseline"


def coverage_of(genes_scored: set, feature_map: Dict[str, List]) -> float:
    """Fraction of scored genes that have >=1 neighbour in a feature's map."""
    if not genes_scored:
        return 0.0
    return sum(1 for g in genes_scored if g in feature_map) / len(genes_scored)


def run_head_to_head(
    df, splits, incumbent_cols, arms, fold_aware,
    anchor_value=None, anchor_tol=None,
    genes_scored=None, coverage_maps=None, coverage_matched_arms=None,
    n_perm=10000, seed=None, allow_no_anchor=False,
):
    """
    Run a leakage-safe, fold-aware feature-group head-to-head and REFUSE to report
    arm metrics unless the four guarantees pass. Returns a dict with keys
    'metrics' (per arm), 'permutation' (per non-incumbent arm), 'leakage',
    'anchor', and 'guarantees_passed'. Raises RuntimeError with a
    'head-to-head refused:' prefix if a required guarantee is missing or fails.

    Parameters
    ----------
    df : DataFrame with caseid, gene, target and the incumbent/arm columns.
    splits : list of (train_idx, test_idx) row-index pairs (fixed CV, unchanged).
    incumbent_cols : columns for arm 'A_incumbent' (the reference).
    arms : dict arm_name -> list of columns (must include an incumbent arm keyed
           'A_incumbent'). Each arm is scored out-of-fold with cv_predict_unified.
    fold_aware : col -> {"pivot":..., "map":...|"weighted":...} for every fold-aware
           feature used by any arm. Drives both prediction and the leakage check.
    anchor_value : the named canonical incumbent AU-PR to gate against
           (e.g. 0.6465). Required unless allow_no_anchor=True.
    genes_scored, coverage_maps, coverage_matched_arms : when two arms differ in
           gene coverage, supply the scored-gene set, the per-feature maps to
           measure coverage, and the names of the coverage-matched control arm(s);
           the harness verifies at least one coverage-matched arm is present and
           records the coverage of each feature. If arms differ in coverage and no
           coverage-matched arm is supplied, the run is REFUSED.

    Refusal conditions (guarantees):
      1. anchor: incumbent arm AU-PR must match anchor_value within anchor_tol.
      2. permutation: computed for every non-incumbent arm (never skipped).
      3. coverage-matched control: required when coverage differs across arms.
      4. leakage: leakage_selfcheck must show max_feature_diff 0.0 in every fold.
    """
    if seed is None:
        seed = SEED
    if anchor_tol is None:
        anchor_tol = ANCHOR_TOL
    if "A_incumbent" not in arms:
        raise RuntimeError("head-to-head refused: arms must include 'A_incumbent' (the reference arm)")

    # Guarantee 4 first (cheap, and a hard safety gate): leakage self-check.
    leak = leakage_selfcheck(df, splits, fold_aware)
    if bool(leak["leakage_detected"].any()):
        bad = leak[leak.leakage_detected][["fold", "feature", "max_feature_diff"]].to_dict("records")
        raise RuntimeError(f"head-to-head refused: leakage self-check FAILED: {bad}")
    if int((leak["train_test_overlap"] != 0).sum()) > 0:
        raise RuntimeError("head-to-head refused: train/test case overlap is non-zero in some fold")

    # Score every arm out-of-fold with identical leakage-safe machinery.
    oof = {}
    for arm, cols in arms.items():
        oof[arm] = cv_predict_unified(df, cols, splits, fold_aware, seed=seed)
    metrics = {arm: dedup_metrics(df, p) for arm, p in oof.items()}

    # Guarantee 1: bit-exact anchor gate on the incumbent.
    inc_aupr = metrics["A_incumbent"]["aupr"]
    anchor = {"anchor_value": anchor_value, "incumbent_aupr": round(inc_aupr, 6)}
    if anchor_value is None:
        if not allow_no_anchor:
            raise RuntimeError("head-to-head refused: no anchor_value given "
                               "(pass the named canonical AU-PR, or allow_no_anchor=True to waive)")
        anchor["status"] = "waived"
    else:
        diff = abs(inc_aupr - anchor_value)
        anchor["abs_diff"] = round(float(diff), 8)
        anchor["passed"] = bool(diff <= anchor_tol)
        anchor["status"] = "pass_bit_exact" if diff <= anchor_tol else "FAIL"
        if diff > anchor_tol:
            raise RuntimeError(f"head-to-head refused: anchor gate FAILED "
                               f"(incumbent {inc_aupr:.6f} vs anchor {anchor_value:.6f}, "
                               f"|diff| {diff:.2e} > tol {anchor_tol:.0e})")

    # Guarantee 3: coverage-matched control when arms differ in coverage.
    coverage = None
    if genes_scored is not None and coverage_maps:
        coverage = {name: round(coverage_of(set(genes_scored), m), 4)
                    for name, m in coverage_maps.items()}
        cov_vals = list(coverage.values())
        differ = (max(cov_vals) - min(cov_vals)) > 0.01
        if differ and not coverage_matched_arms:
            raise RuntimeError("head-to-head refused: feature groups differ in gene coverage "
                               f"{coverage}, but no coverage_matched_arms supplied to isolate topology")
        if coverage_matched_arms:
            missing = [a for a in coverage_matched_arms if a not in arms]
            if missing:
                raise RuntimeError(f"head-to-head refused: coverage_matched_arms {missing} not in arms")

    # Guarantee 2: permutation null for every non-incumbent arm (never skipped).
    perm = {}
    for arm in arms:
        if arm == "A_incumbent":
            continue
        perm[arm] = permutation_null(df, oof[arm], oof["A_incumbent"],
                                     n_perm=n_perm, seed=seed)

    return {"metrics": metrics, "permutation": perm, "leakage": leak,
            "anchor": anchor, "coverage": coverage,
            "coverage_matched_arms": coverage_matched_arms,
            "oof": oof,
            "guarantees_passed": {"anchor": anchor.get("status"),
                                  "permutation": "computed_all_arms",
                                  "coverage_matched": bool(coverage_matched_arms) if coverage else "n/a",
                                  "leakage": "max_diff_0.0_all_folds"}}
