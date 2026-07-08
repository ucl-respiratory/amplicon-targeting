#!/usr/bin/env python3
# =============================================================================
# 17_ensembl_coords.py  --  Ensembl gene coordinates for the CPTAC gene universe
# -----------------------------------------------------------------------------
# SOURCE-LEVEL DOWNLOAD (belongs in data_download, not analysis). Genomic
# coordinates (chrom, start, end, strand) are a reference annotation, so they
# are fetched ONCE here and cached; the analysis pipeline reads the cache and
# never touches the network.
#
# For every gene symbol in the assembled feature table (stage 13
# str_omic_table_rebuilt.csv) we POST to the Ensembl REST symbol-lookup endpoint
# in chunks, with 429 back-off and a failed-chunk retry loop so the table is
# COMPLETE and deterministic (a truncated coordinate set silently collapses the
# analysis ATAC / gene-density joins). Genes without a current Ensembl symbol
# lookup are simply absent (expected for a small tail of aliases/withdrawn IDs).
#
# Emits: OUT/tables/ensembl_gene_coords.parquet
#          columns: gene, chrom, start, end, strand
#        + a source_manifest row. GRCh38 assembly, Ensembl REST (current release).
# =============================================================================
import os, sys, json, time, csv
from pathlib import Path

def _out_tables():
    d = os.environ.get("CNT_DATA_OUT")
    if d: return Path(d)
    fs_out = os.environ.get("FS_OUT", os.path.join(os.getcwd(), "from_source", "out"))
    return Path(fs_out) / "tables"

DIR_TAB = _out_tables(); DIR_REP = DIR_TAB.parent / "reports"
DIR_TAB.mkdir(parents=True, exist_ok=True); DIR_REP.mkdir(parents=True, exist_ok=True)

REST = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
OUT_PARQUET = DIR_TAB / "ensembl_gene_coords.parquet"
STR_OMIC = DIR_TAB / "str_omic_table_rebuilt.csv"

def _batch(session, symbols):
    return session.post(REST, headers={"Content-Type": "application/json",
                                        "Accept": "application/json"},
                        data=json.dumps({"symbols": symbols}), timeout=120)

def build():
    if OUT_PARQUET.exists() and OUT_PARQUET.stat().st_size > 0:
        print(f"[ensembl] cached: {OUT_PARQUET.name}"); return
    import pandas as pd, requests
    if not STR_OMIC.exists():
        sys.stderr.write(f"[ensembl] {STR_OMIC} missing; run stage 13 first. Skipping.\n"); return
    genes = sorted(pd.read_csv(STR_OMIC, usecols=["gene"]).gene.dropna().unique())
    print(f"[ensembl] resolving coordinates for {len(genes)} genes")
    s = requests.Session()
    coords, pending, rounds, CHUNK = {}, list(genes), 0, 300
    while pending and rounds < 8:
        rounds += 1; fails = []
        for i in range(0, len(pending), CHUNK):
            chunk = pending[i:i+CHUNK]; ok = False
            for attempt in range(4):
                try:
                    r = _batch(s, chunk)
                    if r.status_code == 200:
                        coords.update({k: v for k, v in r.json().items() if v}); ok = True; break
                    elif r.status_code == 429:
                        time.sleep(float(r.headers.get("Retry-After", 2)) + 0.5); continue
                    else:
                        time.sleep(1); continue
                except Exception:
                    time.sleep(1.5)
            if not ok: fails.extend(chunk)
            time.sleep(0.25)
        print(f"[ensembl] round {rounds}: {len(coords)} coords, {len(fails)} pending")
        pending = fails
    rows = [{"gene": k, "chrom": v["seq_region_name"], "start": int(v["start"]),
             "end": int(v["end"]), "strand": int(v["strand"])}
            for k, v in coords.items() if v.get("start") and v.get("seq_region_name") in CHROMS]
    df = pd.DataFrame(rows).sort_values("gene").reset_index(drop=True)
    df.to_parquet(OUT_PARQUET)
    print(f"[ensembl] wrote {len(df)}/{len(genes)} genes -> {OUT_PARQUET.name}")
    # manifest row
    with open(DIR_REP / "source_manifest.csv", "a", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["17_ensembl_coords", "Ensembl REST gene coordinates (GRCh38)",
                    REST, str(OUT_PARQUET.resolve()), OUT_PARQUET.stat().st_size, "",
                    time.strftime("%Y-%m-%d"), f"{len(df)} of {len(genes)} CPTAC genes"])

if __name__ == "__main__":
    build()
