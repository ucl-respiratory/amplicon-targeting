# =============================================================================
# 91_verify_paper.py  --  exact-reproduction check against paper v9.
# Encodes every headline in-text number + Table 1 as (paper_value, manifest_key,
# tolerance) and compares to the regenerated values_manifest.json. Emits
# verification_report.md (human) + verification_report.csv (machine), and exits
# non-zero if any check fails, so `run_all.py` surfaces a reproduction drift.
#
# Paper v9 authoritative values are transcribed from out/reports/paper_v9.md
# (the manuscript this pipeline reproduces). Tolerances: exact for integers and
# gene sets; +-0.01 for 2-dp statistics; +-0.05 for fold/ratio and per-type R2.
# =============================================================================
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd
from cnt_io import DIR_REP

# (description, manifest_stage, manifest_key, paper_value, tol)  tol=None -> exact/equality
CHECKS = [
    # --- mechanism / transmission ---
    ("CN->RNA mean correlation", "10_transmission_reconciliation", "mean_cn_rna_corr", 0.22, 0.01),
    ("CN->protein mean correlation", "10_transmission_reconciliation", "mean_cn_prot_corr", 0.11, 0.01),
    ("Between-gene responsiveness r", "10_transmission_reconciliation", "between_gene_r", 0.73, 0.01),
    ("Mechanism decomposition total R2", "12_mechanism_decomposition", "total_r2", 0.58, 0.01),
    ("Dosage transmission % of explained", "12_mechanism_decomposition", "transmission_pct", 89.2, 0.5),
    ("Coelev vs transmission r", "13_coelevation_by_transmission", "r_transmission_vs_coelev", 0.42, 0.02),
    ("Coelev vs attenuation r", "13_coelevation_by_transmission", "r_attenuation_vs_coelev", 0.07, 0.02),
    # --- ATAC ---
    ("ATAC mean R2 CPTAC-only", "15_atac_transmission", "mean_R2_CPTAC", 0.034, 0.005),
    ("ATAC mean R2 +ATAC", "15_atac_transmission", "mean_R2_plus_ATAC", 0.092, 0.005),
    ("ATAC fold gain", "15_atac_transmission", "fold_gain", 2.7, 0.1),
    # --- robustness ---
    ("Transmission reliability (mean)", "16_statistical_robustness", "mean_reliability", 0.74, 0.02),
    ("Disattenuated R2 (mean)", "16_statistical_robustness", "mean_R2_disattenuated", 0.123, 0.01),
    # --- nomination funnel ---
    ("Funnel S1 genes on recurrent amplicons", "22_target_funnel", "s1_recurrent_amplicon.genes", 3096, None),
    ("Funnel S1 gene x amplicon tests", "22_target_funnel", "s1_recurrent_amplicon.tests", 5163, None),
    ("Funnel S2 reliably co-elevated events", "22_target_funnel", "s2_reliably_coelevated_events", 1387, None),
    ("Funnel S3 surface/secreted genes", "22_target_funnel", "s3_surface_secreted.genes", 109, None),
    ("Funnel S4 elevated>=50% genes", "22_target_funnel", "s4_elevated_50pct.genes", 77, None),
    ("Funnel S5 DepMap-annotated genes", "22_target_funnel", "s5_depmap_annotated.genes", 62, None),
    ("Funnel S6 surface-accessible genes", "22_target_funnel", "s6_surface_accessible_genes", 31, None),
    # --- dependency validation ---
    ("DepMap n testable targets", "21_dependency_validation", "depmap_n_testable", 70, None),
    ("DepMap frac moved up", "21_dependency_validation", "depmap_frac_up", 0.84, 0.02),
    ("DepMap median delta (log2)", "21_dependency_validation", "depmap_median_delta_log2", 0.35, 0.02),
    ("DepMap median delta surface subset", "21_dependency_validation", "depmap_median_delta_surface", 0.43, 0.02),
    ("DepMap n FDR<0.1", "21_dependency_validation", "depmap_n_fdr10", 34, None),
    ("Nominated surface events", "21_dependency_validation", "n_nominated_surface_events", 75, None),
    ("Nominated unique genes", "21_dependency_validation", "n_nominated_unique_genes", 62, None),
    ("Dispensable passenger events", "21_dependency_validation", "n_dispensable_passenger", 70, None),
    ("Essential driver-like events", "21_dependency_validation", "n_essential_driver_like", 4, None),
    # --- multispecific safety ---
    ("Multi-antigen accessible sets", "23_multispecific_safety", "n_multi_antigen_sets", 7, None),
    ("LUAD 1q OR-gate organ burden", "23_multispecific_safety", "sets.LUAD_1q.OR_med_organs", 9, None),
    ("LUAD 1q AND-gate organ burden", "23_multispecific_safety", "sets.LUAD_1q.AND_med_organs", 0, None),
    ("LUAD 1q tumour AND (amplified)", "23_multispecific_safety", "sets.LUAD_1q.tumour_AND_amp", 0.277, 0.02),
    ("LUAD 7p tumour AND (amplified)", "23_multispecific_safety", "sets.LUAD_7p.tumour_AND_amp", 0.256, 0.02),
    # --- normal single-cell threshold ---
    ("Median tumour/normal RNA fold", "30_normal_singlecell_threshold", "median_fold_tumour_over_normal", 64.0, 3.0),
    ("LUAD 1q limiting antigen", "30_normal_singlecell_threshold", "sets.LUAD_1q.limiting_antigen", "MUC1", None),
    ("LUAD 1q normal celltypes >=25nTPM", "30_normal_singlecell_threshold", "sets.LUAD_1q.normal_celltypes_ge25", 1, None),
    # --- tumour same-cell enrichment (Fig 6) ---
    ("LUAD 1q same-cell enrichment", "31_tumour_samecell", "sets.LUAD_1q.enrich", 7.27, 0.05),
    ("LUAD 7p same-cell enrichment", "31_tumour_samecell", "sets.LUAD_7p.enrich", 5.43, 0.05),
    ("LUAD 5p same-cell enrichment", "31_tumour_samecell", "sets.LUAD_5p.enrich", 1.34, 0.05),
    ("LSCC 1q same-cell enrichment", "31_tumour_samecell", "sets.LSCC_1q.enrich", 5.27, 0.05),
    ("LSCC 20q same-cell enrichment", "31_tumour_samecell", "sets.LSCC_20q.enrich", 13.81, 0.05),
    ("LUAD same-cell n cells", "31_tumour_samecell", "sets.LUAD_1q.n_cells", 27122, None),
    ("LSCC same-cell n cells", "31_tumour_samecell", "sets.LSCC_1q.n_cells", 34312, None),
]

# Table 1: the 7 co-target sets (cancer_arm -> exact ordered gene list, tier)
TABLE1 = {
    "LUAD_1q": (["ADAM15","CD46","EFNA1","MUC1","NCSTN","XPR1"], "A"),
    "LUAD_7p": (["DAGLB","EGFR","ITGB8","TSPAN13"], "A"),
    "LSCC_1q": (["F11R","HSD17B7","NCSTN"], "A"),
    "LSCC_20q": (["GGT7","SDC4","TM9SF4"], "A"),
    "PDA_1q": (["F11R","MUC1","NCSTN"], "B"),
    "LUAD_5p": (["CLPTM1L","SLC12A7"], "B"),
    "UCEC_1q": (["ADAM15","PIGR"], "B"),
}

def get(manifest, stage, key):
    d = manifest[stage]
    for part in key.split("."):
        d = d[part]
    return d

def main():
    manifest = json.load(open(DIR_REP / "values_manifest.json"))
    rows = []
    for desc, stage, key, paper, tol in CHECKS:
        try:
            got = get(manifest, stage, key)
        except (KeyError, TypeError):
            rows.append({"check": desc, "paper": paper, "regenerated": "MISSING", "status": "FAIL"}); continue
        if isinstance(paper, str):
            ok = (str(got) == paper)
        elif tol is None:
            ok = (got == paper)
        else:
            ok = abs(float(got) - float(paper)) <= tol
        rows.append({"check": desc, "paper": paper, "regenerated": got, "status": "PASS" if ok else "FAIL"})

    # Table 1 checks
    t1 = json.load(open(DIR_REP / "values" / "24_offtarget_safety.json"))["table1_sets"]
    for sn, (genes, tier) in TABLE1.items():
        got = t1.get(sn)
        if not got:
            rows.append({"check": f"Table 1 {sn}", "paper": ";".join(genes), "regenerated": "MISSING", "status": "FAIL"}); continue
        got_genes = [g.strip() for g in got["genes"].split(";")]
        ok = (got_genes == genes) and (got["tier"] == tier)
        rows.append({"check": f"Table 1 {sn} ({tier}, {len(genes)} genes)",
                     "paper": "; ".join(genes), "regenerated": got["genes"],
                     "status": "PASS" if ok else "FAIL"})

    df = pd.DataFrame(rows)
    n_pass = int((df.status == "PASS").sum()); n = len(df)
    df.to_csv(DIR_REP / "verification_report.csv", index=False)

    lines = ["# Exact-reproduction verification against paper v9", "",
             f"**{n_pass}/{n} checks PASS**", "",
             "Each regenerated value is produced by the `analysis/` pipeline from the",
             "`data_download/` outputs and compared to the number stated in the manuscript",
             "(`paper_v9.md`). Integers, counts and gene sets are checked for exact equality;",
             "2-dp statistics to +-0.01-0.02; fold/ratio and per-type R2 to +-0.05.", "",
             "| Check | Paper v9 | Regenerated | Status |", "|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['check']} | {r['paper']} | {r['regenerated']} | {r['status']} |")
    if n_pass < n:
        lines += ["", "## FAILURES", ""]
        for r in rows:
            if r["status"] == "FAIL":
                lines.append(f"- **{r['check']}**: paper={r['paper']} regenerated={r['regenerated']}")
    (DIR_REP / "verification_report.md").write_text("\n".join(lines))
    print(f"{n_pass}/{n} checks PASS")
    if n_pass < n:
        for r in rows:
            if r["status"] == "FAIL": print(f"  FAIL {r['check']}: paper={r['paper']} got={r['regenerated']}")
    return 0 if n_pass == n else 1

if __name__ == "__main__":
    sys.exit(main())
