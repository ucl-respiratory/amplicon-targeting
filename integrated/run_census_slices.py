import os, sys, csv
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg
import cellxgene_census, tiledbsoma as soma
import pandas as pd, numpy as np, scipy.sparse as sp

CENSUS_VERSION = "2023-12-15"
PROXY = ("localhost", os.environ.get("HTTPS_PROXY","http://localhost:56768").split(":")[-1])
tcfg = {"vfs.s3.region":"us-west-2","vfs.s3.no_sign_request":"true",
        "vfs.s3.proxy_host":PROXY[0],"vfs.s3.proxy_port":PROXY[1],"vfs.s3.proxy_scheme":"http"}
ctx = soma.SOMATileDBContext(tiledb_config=tcfg)

# union of nominated antigens per cancer type, sourced from the pan-cancer
# nomination table (04c) so the census slices always cover the current antigen
# set rather than a stale hardcoded list.
# Malignant-cell slices are extracted for every cohort that (a) has nominated
# antigens and (b) has malignant primary cells in the CELLxGENE census. Of our
# six cohorts, LUAD/LSCC (measured nominations) and GBM (prediction-only, from
# 04c/M2) meet both; CCRCC, UCEC and PDA have no malignant primary cells in the
# census (renal cells are annotated as epithelial subtypes, not a malignant
# compartment) so cannot be tested and are reported as nominated-only.
_ant = pd.read_csv(cfg.DIR_TAB / "adc_target_antigens.csv")
SETS = {code: sorted(_ant.loc[_ant.cohort == code, "antigen"].unique())
        for code in ("LUAD", "LSCC")}
# GBM antigens are prediction-only (no CPTAC GBM proteome) and live in a
# separate table produced by the M2 analysis.
_gbm_csv = cfg.DIR_TAB / "m2_gbm_prediction_only.csv"
if _gbm_csv.exists():
    SETS["GBM"] = sorted(pd.read_csv(_gbm_csv)["antigen"].unique())
print("census gene sets:", {k: len(v) for k, v in SETS.items()})
DISEASE = {"LUAD":"lung adenocarcinoma","LSCC":"squamous cell lung carcinoma",
           "GBM":"glioblastoma"}

od = cfg.PATHS["cxg_luad"].parent; od.mkdir(parents=True, exist_ok=True)
census = cellxgene_census.open_soma(census_version=CENSUS_VERSION, context=ctx)
manifest=[]
with census:
    for code, genes in SETS.items():
        out = od / f"cellxgene_{code}_malignant.parquet"
        if out.exists() and out.stat().st_size>0:
            print(f"cached {out.name}"); continue
        flt = (f"disease == '{DISEASE[code]}' and cell_type == 'malignant cell' "
               "and is_primary_data == True")
        ad = cellxgene_census.get_anndata(
            census, "homo_sapiens", obs_value_filter=flt,
            var_value_filter="feature_name in "+str(genes),
            obs_column_names=["cell_type","donor_id","nnz"],
            var_column_names=["feature_name"])
        X = ad.X.toarray() if sp.issparse(ad.X) else np.asarray(ad.X)
        det = (X>0).astype("uint8")
        df = pd.DataFrame(det, columns=list(ad.var.feature_name))
        df.insert(0,"donor_id",ad.obs["donor_id"].values)
        df.insert(1,"nnz",ad.obs["nnz"].values)
        df.to_parquet(out, index=False)
        manifest.append({"cancer":code,"n_cells":int(ad.n_obs),"n_genes":int(ad.n_vars),
                         "n_donors":int(ad.obs['donor_id'].nunique())})
        print(f"{code}: {ad.n_obs} cells x {ad.n_vars} genes, {df.donor_id.nunique()} donors -> {out.name}")
if manifest:
    pd.DataFrame(manifest).to_csv(od/"cellxgene_manifest.csv", index=False)
    print("manifest written")
