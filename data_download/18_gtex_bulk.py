#!/usr/bin/env python3
# =============================================================================
# 18_gtex_bulk.py  --  GTEx v8 bulk median expression (gene x tissue)
# -----------------------------------------------------------------------------
# SOURCE-LEVEL DOWNLOAD (belongs in data_download, not analysis). The off-target
# safety analysis needs GTEx median TPM across normal tissues for the nominated
# antigens. Rather than the analysis pipeline hitting the GTEx portal API once
# per gene at run time, we download the single canonical bulk matrix ONCE here
# (all genes x 54 tissues) and cache it; analysis subsets its gene set locally.
#
# Source file (GTEx portal, GCS mirror):
#   GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
# Tissue columns are renamed to the portal `tissueSiteDetailId` convention
# ("Heart - Left Ventricle" -> "Heart_Left_Ventricle") so the analysis
# VULN_GTEX organ map keys match without further transformation.
#
# Emits: OUT/tables/gtex_v8_median_tpm.parquet
#          index gene (symbol); one column per tissueSiteDetailId; values = median TPM
#        + a source_manifest row.
# =============================================================================
import os, sys, time, csv
from pathlib import Path

def _out_tables():
    d = os.environ.get("CNT_DATA_OUT")
    if d: return Path(d)
    fs_out = os.environ.get("FS_OUT", os.path.join(os.getcwd(), "from_source", "out"))
    return Path(fs_out) / "tables"

DIR_TAB = _out_tables(); DIR_REP = DIR_TAB.parent / "reports"
DIR_RAW = DIR_TAB.parent.parent / "data" / "raw"
DIR_TAB.mkdir(parents=True, exist_ok=True); DIR_REP.mkdir(parents=True, exist_ok=True)
DIR_RAW.mkdir(parents=True, exist_ok=True)

GTEX_URL = ("https://adult-gtex.storage.googleapis.com/bulk-gex/v8/rna-seq/"
            "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")
GCT_GZ = DIR_RAW / "gtex_v8_gene_median_tpm.gct.gz"
OUT_PARQUET = DIR_TAB / "gtex_v8_median_tpm.parquet"

def _tissue_id(t):
    return t.replace(" - ", "_").replace(" ", "_").replace("(", "").replace(")", "")

def build():
    if OUT_PARQUET.exists() and OUT_PARQUET.stat().st_size > 0:
        print(f"[gtex] cached: {OUT_PARQUET.name}"); return
    import pandas as pd, requests
    if not GCT_GZ.exists() or GCT_GZ.stat().st_size == 0:
        print(f"[gtex] downloading {GCT_GZ.name}")
        r = requests.get(GTEX_URL, timeout=300)
        if r.status_code != 200:
            sys.stderr.write(f"[gtex] download failed ({r.status_code}); skipping.\n"); return
        GCT_GZ.write_bytes(r.content)
    g = pd.read_csv(GCT_GZ, sep="\t", skiprows=2)
    tissue_cols = [c for c in g.columns if c not in ("Name", "Description")]
    g = g.rename(columns={c: _tissue_id(c) for c in tissue_cols})
    ids = [_tissue_id(c) for c in tissue_cols]
    # collapse duplicate symbols (multiple Ensembl IDs) by max median TPM per tissue
    med = g.groupby("Description")[ids].max()
    med.index.name = "gene"
    med.to_parquet(OUT_PARQUET)
    print(f"[gtex] wrote {med.shape[0]} genes x {med.shape[1]} tissues -> {OUT_PARQUET.name}")
    with open(DIR_REP / "source_manifest.csv", "a", newline="") as fh:
        csv.writer(fh).writerow(
            ["18_gtex_bulk", "GTEx v8 bulk median TPM (gene x tissue)", GTEX_URL,
             str(OUT_PARQUET.resolve()), OUT_PARQUET.stat().st_size, "",
             time.strftime("%Y-%m-%d"), f"{med.shape[0]} genes x {med.shape[1]} tissues"])

if __name__ == "__main__":
    build()
