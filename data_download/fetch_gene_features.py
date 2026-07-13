#!/usr/bin/env python3
"""fetch_gene_features.py -- download the PUBLIC gene-property sources used by the
integrated predictor's feature table, into data_download/from_source so the
feature table can be rebuilt from source (no committed gene_intrinsic tables).

Sources (all public, all gene-intrinsic -- none touches CPTAC tumour data):
  - gnomAD v4.1 gene constraint (LOEUF, pLI, mis_z)          -> dosage sensitivity
  - DescribePROT 9606 (length/MW/pI/GRAVY/TM/disorder)       -> biophysics
  - STRING v12 human protein links (degree/centrality)        -> network
  - Ensembl Compara dN/dS + gene-age (orthologue depth)       -> evolution
    (fetched separately via the Ensembl REST helper; see build step)

Writes raw files to RAW/features/ and appends rows to a fetch manifest.
Idempotent: skips a file that already exists with non-zero size unless FORCE=1.
"""
import os, sys, hashlib, time, gzip, io, csv
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parent / "from_source"
RAW  = Path(os.environ.get("CNT_RAW", ROOT / "data" / "raw")) / "features"
RAW.mkdir(parents=True, exist_ok=True)
MANIFEST = ROOT / "out" / "tables" / "feature_source_manifest.csv"
FORCE = os.environ.get("FORCE", "0") == "1"

SOURCES = {
    "gnomad_constraint": {
        "url": "https://gcp-public-data--gnomad.storage.googleapis.com/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv",
        "dst": "gnomad.v4.1.constraint_metrics.tsv",
        "note": "gnomAD v4.1 gene constraint metrics (LOEUF/pLI/mis_z)",
    },
    # NOTE: DescribePROT (disorder/SS descriptors) and Ensembl-Compara dN/dS are
    # NOT live-fetchable: the DescribePROT bulk endpoint now serves a JS single-page
    # app for every path (no direct CSV/JSON), and current Ensembl BioMart no longer
    # exposes dN/dS attributes. Those features (vsl2_disorder, psipred_*, asaquick_buried,
    # dn_ds, gene_age_proxy) are carried as a PINNED SNAPSHOT in
    # data/raw/features/external_predicted_descriptors_snapshot.csv, logged in the
    # manifest with provenance='pinned_snapshot'. They are gene-intrinsic external
    # annotations (no CPTAC content). All other predictor features are rebuilt from
    # the live-fetched sources below plus UniProt (stage 15) and the from_source omics.
    "string_links": {
        "url": "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz",
        "dst": "9606.protein.links.v12.0.txt.gz",
        "note": "STRING v12 human PPI links (network centrality)",
    },
    "string_aliases": {
        "url": "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz",
        "dst": "9606.protein.aliases.v12.0.txt.gz",
        "note": "STRING v12 protein->gene-symbol alias map",
    },
}


def _md5(p, block=1 << 20):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(block), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(key):
    s = SOURCES[key]
    dst = RAW / s["dst"]
    if dst.exists() and dst.stat().st_size > 0 and not FORCE:
        print(f"[skip] {key}: {dst.name} present ({dst.stat().st_size} bytes)")
        return dst
    print(f"[get ] {key}: {s['url']}")
    req = urllib.request.Request(s["url"], headers={"User-Agent": "python-urllib"})
    with urllib.request.urlopen(req, timeout=300) as r, open(dst, "wb") as out:
        while True:
            buf = r.read(1 << 20)
            if not buf:
                break
            out.write(buf)
    print(f"[ok  ] {key}: {dst.stat().st_size} bytes -> {dst}")
    return dst


def write_manifest(rows):
    hdr = ["key", "source", "url", "file", "bytes", "md5", "retrieved_utc"]
    with open(MANIFEST, "w", newline="") as f:
        w = csv.writer(f); w.writerow(hdr)
        for row in rows:
            w.writerow(row)
    print(f"[manifest] {MANIFEST} ({len(rows)} rows)")


def main():
    rows = []
    for key in SOURCES:
        try:
            p = fetch(key)
            rows.append([key, SOURCES[key]["note"], SOURCES[key]["url"],
                         str(p.relative_to(ROOT)), p.stat().st_size, _md5(p),
                         time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())])
        except Exception as e:
            print(f"[FAIL] {key}: {type(e).__name__}: {e}")
            rows.append([key, SOURCES[key]["note"], SOURCES[key]["url"], "MISSING", 0, "", ""])
    write_manifest(rows)


if __name__ == "__main__":
    main()
