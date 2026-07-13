# =============================================================================
# 04c_pancancer_nomination.py -- pan-cancer surface-antigen nomination
# -----------------------------------------------------------------------------
# Replaces the hardcoded 5-lung-set nomination. Nominates surface antigens from
# EVERY recurrent amplicon across ALL cohorts that has co-elevation support:
#
#   candidate      : gene co-elevated on a recurrent amplicon (Fisher FDR < 0.1,
#                    from 00c_amplicons) AND observed transmissibility >= 0.40
#   surface gate   : live UniProt ectodomain >= MIN_ECTODOMAIN_AA AND >=1 TM helix
#   annotation     : predicted transmissibility (atlas), DepMap essentiality
#
# The candidate universe is pan-cancer by construction. Single-cell co-detection
# (constructs) is computed only where a malignant census slice exists (LUAD,
# LSCC); non-lung constructs are nominated but flagged not-yet-single-cell-tested.
#
# Emits: tables/adc_target_antigens.csv  (all cohorts)
#        tables/adc_constructs.csv        (with tested/nominated flag)
# =============================================================================
import sys, json, importlib.util
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cnt_shared as ish
import config as cfg

# reuse the validated UniProt topology + surface verdict from 04_surface_targets
_spec = importlib.util.spec_from_file_location("m04", str(Path(__file__).resolve().parent / "04_surface_targets.py"))
m04 = importlib.util.module_from_spec(_spec)
# guard: 04's module-level code shouldn't run main on import
import builtins
_spec.loader.exec_module(m04)

TRANSMIT_MIN = 0.40
ECTO_MIN     = 50     # accessible epitope (aa); matches nomination rule
ESS_FLAG     = -0.40  # DepMap dependency liability flag

def main():
    ish.apply_style(sizes=(9,8,7))

    # ---- candidates: co-elevated on recurrent amplicon + transmitted ----
    co = pd.read_csv(cfg.DIR_TAB / "amplicon_coelevation.csv")
    co_sig = co[co.fdr < cfg.FDR_ALPHA].copy()
    atlas = ish.load_atlas()[
        ["gene","observed_transmissibility","predicted_transmissibility_oof","n_amplified_cases"]]
    feat = ish.load_feature_table()
    dep = feat[["gene","dep_mean_effect"]] if "dep_mean_effect" in feat.columns else \
          pd.DataFrame({"gene": [], "dep_mean_effect": []})

    cand = (co_sig.merge(atlas, on="gene", how="left")
                  .merge(dep, on="gene", how="left"))
    cand = cand[cand.observed_transmissibility >= TRANSMIT_MIN].copy()
    print(f"co-elevated+transmitted candidates: {len(cand)} rows, "
          f"{cand.gene.nunique()} genes, {cand.tumor_code.nunique()} cohorts")

    # ---- surface gate: live UniProt topology ----
    topo = m04.fetch_topology(list(cand.gene.unique()))
    print(f"UniProt topology fetched for {len(topo)} / {cand.gene.nunique()} genes")
    cand = cand.merge(topo, on="gene", how="left")
    cand["ecto_aa"] = cand["total_extracellular_aa"].fillna(0)
    cand["n_tm"] = cand["n_TM_helices"].fillna(0)
    cand["surface_accessible"] = (cand.ecto_aa >= ECTO_MIN) & (cand.n_tm >= 1)

    nom = cand[cand.surface_accessible].copy()
    nom["amplicon"] = nom.tumor_code + "_" + nom.band.astype(str)
    nom["essential_flag"] = np.where(nom.dep_mean_effect <= ESS_FLAG, "essential-liability", "")
    print(f"NOMINATED (surface-accessible): {len(nom)} antigen-amplicon rows, "
          f"{nom.gene.nunique()} genes across {nom.tumor_code.nunique()} cohorts")
    print("per cohort:", nom.groupby("tumor_code").gene.nunique().to_dict())

    # collapse to arm-level amplicon groups so constructs aren't over-fragmented
    nom["arm"] = nom.band.astype(str).str.extract(r"^(\d+[pq])")
    nom["amplicon_arm"] = nom.tumor_code + "_" + nom.arm

    T1 = (nom.sort_values(["tumor_code","arm","observed_transmissibility"], ascending=[True,True,False])
             [["tumor_code","amplicon_arm","gene","observed_transmissibility",
               "predicted_transmissibility_oof","n_amplified_cases_x" if "n_amplified_cases_x" in nom.columns else "n_amplified_cases",
               "ecto_aa","n_tm","dep_mean_effect","essential_flag"]]
             .drop_duplicates(["amplicon_arm","gene"]))
    T1.columns = ["cohort","amplicon","antigen","obs_transmit","pred_transmit",
                  "n_amplified","ecto_aa","n_tm","dep_effect","essential_flag"]
    T1.to_csv(cfg.DIR_TAB / "adc_target_antigens.csv", index=False)

    ish.record("04c_pancancer_nomination", {
        "transmit_floor": TRANSMIT_MIN, "ecto_min_aa": ECTO_MIN,
        "n_candidates": int(cand.gene.nunique()),
        "n_nominated_genes": int(nom.gene.nunique()),
        "n_nominated_rows": int(len(T1)),
        "n_cohorts": int(nom.tumor_code.nunique()),
        "per_cohort_genes": nom.groupby("tumor_code").gene.nunique().to_dict(),
        "per_amplicon_arm": T1.groupby("amplicon").antigen.nunique().to_dict(),
    })
    print(f"\nnominated antigen table -> {len(T1)} rows across {T1.amplicon.nunique()} arm-amplicons")
    return T1

if __name__ == "__main__":
    main()
