---
name: feature-group-head-to-head
description: Leakage-safe, fold-aware head-to-head comparison of feature groups on fixed CV splits, for copy-number-to-protein transmission models (CN_targeting). Use when comparing whether one feature group (for example a STRING PPI neighbourhood) improves out-of-fold AU-PR over an incumbent (for example a KEGG neighbourhood), or for any leave-tumour-out or feature-swap comparison that must be reproducible and leakage-safe. Encodes four guarantees and refuses to report if any is skipped.
---

# feature-group-head-to-head

A reproducible harness for comparing feature groups in the CN_targeting
copy-number-to-protein transmission model. It does not merely enable good
practice; it ENCODES four guarantees and `run_head_to_head` raises
`RuntimeError("head-to-head refused: ...")` if any applicable guarantee is
skipped or fails.

British English, no em dashes. Deterministic: seed 2, XGBoost 600 trees,
single-thread. Metrics are always deduplicated to unique (caseid, gene) and
report Random-AP (base rate) beside AU-PR as lift.

## The four guarantees (refusal conditions)

1. **Bit-exact anchor gate.** The incumbent arm (`A_incumbent`) must reproduce a
   NAMED canonical AU-PR within `anchor_tol` (default 2e-4) before any arm is
   compared. Pass `anchor_value=` (for example 0.6465 for the ALL cohort, 0.7515
   for 3Q). Refuses if the incumbent does not reproduce it, or if no
   `anchor_value` is given (waive explicitly with `allow_no_anchor=True` only
   when there is genuinely no canonical number, for example a new LTO fold).
2. **Label-shuffle permutation null on FIXED out-of-fold predictions**, models
   held fixed. Computed for EVERY non-incumbent arm, never skipped. Reported
   symmetrically per the standing statistical rule: a wrong-side-of-zero delta
   yields p approaching 1.0 and the verdict "no evidence the arm improves on the
   baseline" - NEVER "significantly negative". Effect size (z) is reported beside
   the p-value.
3. **Coverage-matched confound control.** When two feature groups differ in gene
   coverage (measured from `coverage_maps` over `genes_scored`), a coverage-matched
   arm MUST be supplied via `coverage_matched_arms=` (an arm whose network feature
   is restricted to genes the incumbent also covers, so topology is compared on an
   identical gene set). Refuses if coverage differs by more than 1 percentage
   point and no coverage-matched arm is present.
4. **Leakage self-check across all folds and fold-aware features.** Before any arm
   is scored, every fold-aware feature is stress-tested: all validation cases'
   source-pivot values are corrupted to a sentinel and the recomputed per-gene
   fold features must be byte-identical (max feature diff 0.0). Refuses if any
   fold shows non-zero change, or if any fold has non-zero train/test case
   overlap.

## Contract

**Inputs**
- `df`: DataFrame with `caseid`, `gene`, `target`, and the incumbent columns.
- `splits`: list of (train_idx, test_idx) row-index pairs (the FIXED CV splits,
  unchanged; reconstruct from the frozen `cv_splits.json` by caseid).
- `arms`: dict `arm_name -> [columns]`. MUST include `A_incumbent` (the reference).
- `fold_aware`: dict `col -> {"pivot": caseid x gene DataFrame, "map": {gene:[nbr]}}`
  (unweighted) or `{"pivot":..., "weighted": {gene:[(nbr, weight)]}}` (confidence
  weighted). Every fold-aware column any arm uses must appear here; it drives both
  prediction and the leakage check. A column present in `df` and not in
  `fold_aware` is used static (matching the canonical treatment of `nbr_rna`).
- `anchor_value`, `anchor_tol`: the named canonical incumbent AU-PR and tolerance.
- `genes_scored`, `coverage_maps`, `coverage_matched_arms`: for guarantee 3.
- `n_perm` (default 10000), `seed` (default 2).

**Outputs** (dict)
- `metrics`: per arm, `{aupr, random_ap, lift, auroc, n_unique, n_raw}`.
- `permutation`: per non-incumbent arm, `{obs_delta, null_mean, null_std, z,
  p_value, p_floor, n_perm, verdict}`.
- `leakage`: per (fold, feature) DataFrame with `max_feature_diff` and
  `leakage_detected`.
- `anchor`: `{anchor_value, incumbent_aupr, abs_diff, passed, status}`.
- `coverage`: per-feature coverage fraction over `genes_scored`.
- `oof`: per-arm out-of-fold probability vectors (for downstream figures).
- `guarantees_passed`: a summary of which guarantee checks ran and passed.

## Usage

The kernel sidecar loads automatically. Build the fold-aware spec, define arms
including `A_incumbent`, then call `run_head_to_head`:

```python
fold_aware = {
    "nbr_cn_adjusted": {"pivot": cn_pivot, "map": kegg_map},         # incumbent KEGG nbr
    "net_cn":  {"pivot": cn_pivot,  "weighted": string_weighted},   # STRING net
    "net_rna": {"pivot": rna_pivot, "weighted": string_weighted},
    # coverage-matched variants restrict the neighbour map to KEGG-covered genes
    "net_cn_cm":  {"pivot": cn_pivot,  "weighted": string_weighted_kegg_covered},
    "net_rna_cm": {"pivot": rna_pivot, "weighted": string_weighted_kegg_covered},
}
arms = {
    "A_incumbent":            ["nbr_cn_adjusted", "rna", "nbr_rna"],
    "B_kegg_plus_string":     ["nbr_cn_adjusted", "rna", "nbr_rna", "net_cn", "net_rna"],
    "C_string_swaps_kegg":    ["net_cn", "rna", "net_rna"],
    "Cprime_coverage_matched":["net_cn_cm", "rna", "net_rna_cm"],
}
res = run_head_to_head(
    df, splits, incumbent_cols=arms["A_incumbent"], arms=arms,
    fold_aware=fold_aware, anchor_value=0.7515, anchor_tol=2e-4,
    genes_scored=set(df["gene"]), coverage_maps={"KEGG": kegg_map, "STRING": string_map},
    coverage_matched_arms=["Cprime_coverage_matched"], n_perm=10000, seed=2,
)
```

Individual helpers are also exposed for bespoke designs: `cv_predict_unified`,
`dedup_metrics`, `permutation_null`, `leakage_selfcheck`, `fold_gene_means`,
`build_nbr_per_gene`, `build_net_per_gene_weighted`, `coverage_of`. If you use
these directly, you are responsible for the guarantees; prefer `run_head_to_head`.

## Data locations (CN_targeting)

- Feathers: `Cache/fs_filtered_{all,3q}__*.feather`; apply the R4 filter
  (cn >= round(ploidy)+1 AND cn >= 3) and deduplicate to unique (caseid, gene).
- Frozen splits: `Cache/tournament_master/{all,3q}/cv_splits.json` (by caseid).
- Source pivots (caseid x gene): `Cache/pipeline_v2/_pivot_unfiltered_nbr_{cn_adjusted,rna,meth}.parquet`.
- KEGG map: `Data/kegg_top10_neighbours.json`. Canonical anchors: ALL 0.6465, 3Q 0.7515.

## Validation

Validated by re-running the completed network head-to-head (3Q cohort) through
`run_head_to_head`: anchor gate bit-exact (incumbent 0.751542 vs 0.7515,
|diff| 4.2e-5); B +0.0486, C +0.0149, C' -0.0145 (verdict "no evidence the arm
improves on the baseline", p approaching 1.0); leakage max diff 0.0 across all
folds and features; all four refusal conditions confirmed to fire. See
skill_card.md.
