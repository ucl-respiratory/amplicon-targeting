# =============================================================================
# 24_offtarget_safety.py  --  off-target safety table + Table 1.
# For each nominated surface target: GTEx breadth across ADC-vulnerable normal
# organs, CPTAC matched tumour-vs-normal protein delta, HPA specificity, and a
# co-elevation "window". Emits adc_offtarget_safety.csv and table1_cotarget_sets.csv
# (the 7 multi-antigen accessible sets with confidence tiers A/B).
# =============================================================================
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import cnt_io, _depmap
from cnt_io import DIR_TAB, PATHS, record_values

VULN_GTEX = {
    "Heart": ["Heart_Left_Ventricle","Heart_Atrial_Appendage"], "Liver": ["Liver"], "Lung": ["Lung"],
    "Kidney": ["Kidney_Cortex","Kidney_Medulla"], "Brain": ["Brain_Cortex","Brain_Cerebellum","Brain_Hippocampus"],
    "GI_tract": ["Colon_Sigmoid","Colon_Transverse","Small_Intestine_Terminal_Ileum","Stomach","Esophagus_Mucosa"],
    "Bone_marrow_blood": ["Whole_Blood","Spleen"], "Pancreas": ["Pancreas"],
    "Skin": ["Skin_Sun_Exposed_Lower_leg"], "Nerve": ["Nerve_Tibial"]}
# Table 1 tiers: A = 3+ dispensable antigens with same-cell confirmation (lung);
# B = 2-antigen pairs or types lacking single-cell data (PDA, UCEC).
TIER = {"LUAD_1q":"A","LUAD_7p":"A","LSCC_1q":"A","LSCC_20q":"A","PDA_1q":"B","LUAD_5p":"B","UCEC_1q":"B"}
SET_ORDER = ["LUAD_1q","LUAD_7p","LSCC_1q","LSCC_20q","PDA_1q","LUAD_5p","UCEC_1q"]

def _load_gtex_matrix():
    """GTEx v8 median TPM (gene x tissueSiteDetailId) from the data_download
    stage-18 cache. The analysis pipeline never touches the network; subsetting
    to the nominated genes happens locally."""
    p = PATHS["gtex_bulk"]
    if not p.exists():
        raise FileNotFoundError(
            f"{p} missing. Run data_download stage 18 (18_gtex_bulk.py) to fetch "
            "and cache GTEx v8 bulk median TPM before running analysis.")
    return pd.read_parquet(p)

def main():
    gs = _depmap.cotarget_dependency_annotated()
    genes = sorted(gs.gene.unique())
    mat = _load_gtex_matrix()
    piv = mat.loc[mat.index.intersection(genes)]
    avail = {k: [t for t in v if t in piv.columns] for k, v in VULN_GTEX.items()}

    rec = []
    for g in piv.index:
        row = piv.loc[g]
        rec.append({"gene": g, "gtex_max_tissue": row.idxmax(), "gtex_max_TPM": round(row.max(),1),
                    "frac_tissues_expr_ge10": round((row >= 10).mean(),2),
                    **{f"gtex_{k}": round(row[v].max(),1) if v else np.nan for k,v in avail.items()}})
    gtx = pd.DataFrame(rec).set_index("gene")
    vuln_cols = [c for c in gtx.columns if c.startswith("gtex_") and c not in ("gtex_max_tissue","gtex_max_TPM")]
    gtx["max_vuln_organ_TPM"] = gtx[vuln_cols].max(axis=1)
    gtx["worst_vuln_organ"] = gtx[vuln_cols].idxmax(axis=1).str.replace("gtex_","")

    # CPTAC matched tumour-vs-normal protein
    prot = pd.read_csv(PATHS["cptac_matched_protein"], index_col=0)
    phe = pd.read_csv(PATHS["cptac_matched_pheno"]).set_index("case_id").loc[prot.columns]
    crec = []
    for g in gtx.index:
        if g not in prot.index: crec.append({"gene": g}); continue
        v = prot.loc[g]; tv = v[(phe.Group=="Tumor").values].dropna(); nv = v[(phe.Group=="Normal").values].dropna()
        try: p = mannwhitneyu(tv, nv, alternative="greater")[1]
        except Exception: p = np.nan
        crec.append({"gene": g, "cptac_tumor_med": round(tv.median(),2), "cptac_normal_med": round(nv.median(),2),
                     "cptac_tumor_vs_normal_delta": round(tv.median()-nv.median(),2), "cptac_tn_p": p})
    cp = pd.DataFrame(crec).set_index("gene")

    ax1 = gs.groupby("gene").agg(p_high_amp=("p_high_amp","max"), p_high_noamp=("p_high_noamp","min"),
                                 best_fdr=("fdr","min"), lift=("lift","max"))
    S = gtx.join(cp).join(ax1)
    S["window"] = (S.p_high_amp - S.p_high_noamp).round(2)

    def verdict(r):
        strong = (r.p_high_amp >= 0.6) and (r.window >= 0.25)
        if strong and r.max_vuln_organ_TPM < 50:  return "A_best_window_low_burden"
        if strong and r.max_vuln_organ_TPM < 150:  return "B_good_window_moderate_burden"
        if r.max_vuln_organ_TPM < 50:              return "C_low_burden_modest_window"
        return "D_high_normal_burden"
    S["adc_safety_tier"] = S.apply(verdict, axis=1)
    S.reset_index().rename(columns={"index":"gene"}).to_csv(DIR_TAB / "adc_offtarget_safety.csv", index=False)

    # ---- Table 1: the 7 multi-antigen accessible sets ----
    sets_csv = DIR_TAB / "multispecific_andgate_sets.csv"
    GA = pd.read_csv(sets_csv)
    t1 = []
    for sn in [s for s in SET_ORDER if s in GA.set.values]:
        row = GA[GA.set == sn].iloc[0]; tc, arm = sn.split("_")
        gl = row.genes.split(";")
        t1.append({"Cancer": tc, "Arm": arm, "n_accessible_antigens": len(gl),
                   "Co_target_set": "; ".join(gl), "Tier": TIER[sn]})
    table1 = pd.DataFrame(t1)
    table1.to_csv(DIR_TAB / "table1_cotarget_sets.csv", index=False)

    record_values("24_offtarget_safety", {
        "n_targets_with_gtex": int(len(gtx)),
        "safety_tier_counts": {k:int(v) for k,v in S.adc_safety_tier.value_counts().items()},
        "table1_n_sets": int(len(table1)),
        "table1_sets": {f"{r.Cancer}_{r.Arm}": {"n": int(r.n_accessible_antigens),
                        "genes": r.Co_target_set, "tier": r.Tier} for r in table1.itertuples()},
    })
    print(table1.to_string(index=False))
    print("\nsafety tiers:", dict(S.adc_safety_tier.value_counts()))

if __name__ == "__main__":
    main()
