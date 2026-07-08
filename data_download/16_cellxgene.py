#!/usr/bin/env python3
# =============================================================================
# 16_cellxgene.py  --  CELLxGENE census malignant-cell slices (tumour scRNA)
# -----------------------------------------------------------------------------
# Downloads the on-tumour single-cell RNA slices used by paper Fig 6 (same-cell
# co-expression). For each testable cancer type we pull malignant cells (primary
# data) x the antigen genes of that type's co-target sets, and cache a compact
# per-type parquet of the detection matrix (cells x genes, plus donor_id + nnz).
#
# This is the ONLY Python stage of data_download (the census API is Python-only,
# via cellxgene-census / TileDB-SOMA). It is optional: if the census is
# unreachable the stage warns and exits 0, leaving Fig 6 to be skipped.
#
# REPRODUCIBILITY
#   - census_version is PINNED (CENSUS_VERSION below); the census is versioned
#     and immutable, so a pinned slice is byte-reproducible.
#   - TileDB's S3 VFS bypasses HTTP(S)_PROXY env vars. In a sandboxed network
#     the proxy must be passed explicitly via tiledb config (see PROXY_* env).
#
# ENV
#   CNT_DATA_OUT   output dir (default: $FS_OUT/tables or ./from_source/out/tables)
#   CENSUS_VERSION override the pinned census version
#   TILEDB_PROXY_HOST / TILEDB_PROXY_PORT   S3 VFS proxy (default localhost:49229)
#     In a sandboxed network, set these to the HTTP proxy that HTTP_PROXY points
#     at (e.g. host/port parsed from $HTTP_PROXY, commonly localhost:<port>);
#     TileDB's S3 VFS does NOT honour HTTP_PROXY automatically. On an open
#     network leave TILEDB_PROXY_HOST empty to connect to S3 directly.
#
# Emits (to OUT/tables/):
#   cellxgene_<TYPE>_malignant.parquet   one per {LUAD, LSCC}
#   cellxgene_manifest.csv               provenance (version, cell/donor counts)
# =============================================================================
import os, sys, json, csv

CENSUS_VERSION = os.environ.get("CENSUS_VERSION", "2025-11-08")

# The five testable co-target sets (LUAD + LSCC have malignant cells in census;
# PDA and UCEC have none). Genes are unioned per cancer type for the query.
SETS = {
    "LUAD": {  # LUAD_1q + LUAD_7p + LUAD_5p
        "genes": ["ADAM15","CD46","EFNA1","MUC1","NCSTN","XPR1",
                  "DAGLB","EGFR","ITGB8","TSPAN13","CLPTM1L","SLC12A7"],
        "disease": "lung adenocarcinoma",
        # restrict to malignant cells sampled from lung (excludes metastatic /
        # other-site LUAD cells) to match the paper's Fig 6 cohort exactly
        "tissue_general": "lung",
    },
    "LSCC": {  # LSCC_1q + LSCC_20q
        "genes": ["F11R","HSD17B7","NCSTN","GGT7","SDC4","TM9SF4"],
        "disease": "squamous cell lung carcinoma",
    },
}

def out_dir():
    d = os.environ.get("CNT_DATA_OUT")
    if d: return d
    fs_out = os.environ.get("FS_OUT", os.path.join(os.getcwd(), "from_source", "out"))
    return os.path.join(fs_out, "tables")

def build():
    try:
        import cellxgene_census
        import tiledbsoma as soma
        import pandas as pd, numpy as np, scipy.sparse as sp
    except ImportError as e:
        sys.stderr.write(f"[cellxgene] missing dependency ({e}); "
                         "install cellxgene-census to enable Fig 6. Skipping.\n")
        return 0

    od = out_dir(); os.makedirs(od, exist_ok=True)
    # TileDB S3 VFS proxy (TileDB ignores HTTP(S)_PROXY env vars).
    proxy_host = os.environ.get("TILEDB_PROXY_HOST", "localhost")
    proxy_port = os.environ.get("TILEDB_PROXY_PORT", "49229")
    tiledb_cfg = {
        "vfs.s3.region": "us-west-2",
        "vfs.s3.no_sign_request": "true",
    }
    if proxy_host:
        tiledb_cfg.update({"vfs.s3.proxy_host": proxy_host,
                           "vfs.s3.proxy_port": proxy_port,
                           "vfs.s3.proxy_scheme": "http"})
    ctx = soma.SOMATileDBContext(tiledb_config=tiledb_cfg)

    manifest = []
    try:
        census = cellxgene_census.open_soma(census_version=CENSUS_VERSION, context=ctx)
    except Exception as e:
        sys.stderr.write(f"[cellxgene] cannot open census {CENSUS_VERSION} ({e}); skipping.\n")
        return 0

    with census:
        for code, spec in SETS.items():
            out = os.path.join(od, f"cellxgene_{code}_malignant.parquet")
            if os.path.exists(out) and os.path.getsize(out) > 0:
                sys.stderr.write(f"[cellxgene] cached: {os.path.basename(out)}\n"); continue
            genes = spec["genes"]
            flt = (f"disease == '{spec['disease']}' and cell_type == 'malignant cell' "
                   "and is_primary_data == True")
            if spec.get("tissue_general"):
                flt += f" and tissue_general == '{spec['tissue_general']}'"
            try:
                ad = cellxgene_census.get_anndata(
                    census, "homo_sapiens",
                    obs_value_filter=flt,
                    var_value_filter="feature_name in " + str(genes),
                    obs_column_names=["cell_type","donor_id","nnz"],
                    var_column_names=["feature_name"])
            except Exception as e:
                sys.stderr.write(f"[cellxgene] {code} query failed ({e}); skipping type.\n"); continue
            X = ad.X.toarray() if sp.issparse(ad.X) else np.asarray(ad.X)
            det = (X > 0).astype("uint8")
            df = pd.DataFrame(det, columns=list(ad.var.feature_name))
            df.insert(0, "donor_id", ad.obs["donor_id"].values)
            df.insert(1, "nnz", ad.obs["nnz"].values)
            df.to_parquet(out, index=False)
            manifest.append({"cancer": code, "disease": spec["disease"],
                             "census_version": CENSUS_VERSION,
                             "n_cells": int(ad.n_obs), "n_genes": int(ad.n_vars),
                             "n_donors": int(ad.obs["donor_id"].nunique()),
                             "file": os.path.basename(out)})
            sys.stderr.write(f"[cellxgene] {code}: {ad.n_obs} cells x {ad.n_vars} genes "
                             f"({df['donor_id'].nunique()} donors) -> {os.path.basename(out)}\n")
    if manifest:
        mpath = os.path.join(od, "cellxgene_manifest.csv")
        with open(mpath, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(manifest[0].keys()))
            w.writeheader(); w.writerows(manifest)
        sys.stderr.write(f"[cellxgene] manifest -> {mpath}\n")
    return 0

if __name__ == "__main__":
    sys.exit(build())
