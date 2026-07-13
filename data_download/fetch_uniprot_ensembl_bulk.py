#!/usr/bin/env python3
"""Bulk-fetch UniProt (sequences + annotations) and Ensembl (mRNA features) for
the integrated predictor feature table. One stream query each -- no per-gene calls.

UniProt: all human reviewed proteins with sequence, keywords, subcellular location,
transmembrane + signal features, protein names -> biophysics (via ProtParam),
is_tf/kinase/receptor/enzyme, tm_domain_count, signal_peptide.

Ensembl BioMart: per-gene canonical transcript length, GC%, transcript count
-> transcript_length, gc_content, n_isoforms. (UTR lengths + codon optimality are
derived downstream from the canonical CDS where available; otherwise left NaN.)
"""
import os, gzip, io, time
from pathlib import Path
import urllib.request, urllib.parse

ROOT = Path(__file__).resolve().parent / "from_source"
OUT = ROOT / "data" / "raw" / "features"
OUT.mkdir(parents=True, exist_ok=True)


def fetch_uniprot():
    dst = OUT / "uniprot_human_reviewed.tsv.gz"
    if dst.exists() and dst.stat().st_size > 0:
        print(f"[skip] uniprot: {dst.stat().st_size} bytes"); return
    fields = ["accession", "gene_primary", "sequence", "keyword",
              "cc_subcellular_location", "ft_transmem", "ft_signal",
              "protein_name", "go_f", "ec"]
    base = "https://rest.uniprot.org/uniprotkb/stream"
    params = {"query": "organism_id:9606 AND reviewed:true",
              "fields": ",".join(fields), "format": "tsv", "compressed": "true"}
    url = base + "?" + urllib.parse.urlencode(params)
    print(f"[get ] uniprot stream: {url[:90]}...")
    req = urllib.request.Request(url, headers={"User-Agent": "python-urllib"})
    with urllib.request.urlopen(req, timeout=600) as r, open(dst, "wb") as f:
        while True:
            b = r.read(1 << 20)
            if not b: break
            f.write(b)
    print(f"[ok  ] uniprot: {dst.stat().st_size} bytes -> {dst}")


def fetch_ensembl_biomart():
    dst = OUT / "ensembl_mrna_features.tsv"
    if dst.exists() and dst.stat().st_size > 0:
        print(f"[skip] ensembl: {dst.stat().st_size} bytes"); return
    # BioMart XML: gene symbol, transcript length, GC%, transcript count-per-gene proxy
    q = ('<?xml version="1.0" encoding="UTF-8"?>'
         '<!DOCTYPE Query>'
         '<Query virtualSchemaName="default" formatter="TSV" header="1" '
         'uniqueRows="0" count="0" datasetConfigVersion="0.6">'
         '<Dataset name="hsapiens_gene_ensembl" interface="default">'
         '<Filter name="transcript_is_canonical" value="only"/>'
         '<Attribute name="external_gene_name"/>'
         '<Attribute name="transcript_length"/>'
         '<Attribute name="percentage_gene_gc_content"/>'
         '<Attribute name="ensembl_transcript_id"/>'
         '</Dataset></Query>')
    url = "https://www.ensembl.org/biomart/martservice?query=" + urllib.parse.quote(q)
    print(f"[get ] ensembl biomart canonical transcripts...")
    req = urllib.request.Request(url, headers={"User-Agent": "python-urllib"})
    with urllib.request.urlopen(req, timeout=600) as r, open(dst, "wb") as f:
        f.write(r.read())
    print(f"[ok  ] ensembl: {dst.stat().st_size} bytes -> {dst}")
    # also n_isoforms: count of all transcripts per gene
    dst2 = OUT / "ensembl_transcript_counts.tsv"
    if not (dst2.exists() and dst2.stat().st_size > 0):
        q2 = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
              '<Query virtualSchemaName="default" formatter="TSV" header="1" '
              'uniqueRows="0" count="0" datasetConfigVersion="0.6">'
              '<Dataset name="hsapiens_gene_ensembl" interface="default">'
              '<Attribute name="external_gene_name"/>'
              '<Attribute name="ensembl_transcript_id"/>'
              '</Dataset></Query>')
        url2 = "https://www.ensembl.org/biomart/martservice?query=" + urllib.parse.quote(q2)
        print(f"[get ] ensembl transcript counts...")
        req2 = urllib.request.Request(url2, headers={"User-Agent": "python-urllib"})
        with urllib.request.urlopen(req2, timeout=600) as r, open(dst2, "wb") as f:
            f.write(r.read())
        print(f"[ok  ] ensembl counts: {dst2.stat().st_size} bytes")


if __name__ == "__main__":
    fetch_uniprot()
    fetch_ensembl_biomart()
    print("[done] bulk fetch complete")
