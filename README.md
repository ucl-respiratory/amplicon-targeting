# CN-targeting: copy-number–driven multi-antigen co-target nomination

Reproducible code for the manuscript *"Dosage transmission, not complex
buffering, is the primary determinant of protein-level co-elevation on cancer
amplicons"* — a pipeline that predicts which co-amplified genes present as
elevated protein and nominates multi-antigen ADC co-target sets.

The repository is two pipelines you run in order:

```
data_download/   →   analysis/
(fetch + cache        (reproduce every figure
 all source data)      and in-text value)
```

## 1. `data_download/` — everything from source
An 18-stage pipeline (R + two Python stages) that downloads, parses and caches
**every** input from its primary source: CPTAC proteogenomics (PDC + GDC), CORUM
complexes, DepMap 23Q4 CRISPR dependency, the TCGA ATAC-seq atlas, UniProt
membrane topology, HPA, GTEx v8, Ensembl gene coordinates, the ADC Atlas, and CZ
CELLxGENE malignant-cell slices. All network access lives here; releases are
pinned (GTEx v8, CELLxGENE census `2025-11-08`, HPA v24, DepMap 23Q4) and every
source is logged to `source_manifest.csv`.

```bash
cd data_download
Rscript run_all.R                 # from scratch: download + assemble
# outputs land in data_download/from_source/{tables,data,reports}
```
See [`data_download/DOWNLOAD_PIPELINE.md`](data_download/DOWNLOAD_PIPELINE.md)
for the stage-by-stage source map and reproducibility notes.

## 2. `analysis/` — reproduce the paper
A 16-stage Python pipeline that assumes `data_download` has run. It **reads only
the cached tables** — no analysis stage touches the network — and regenerates
all 13 figures, Table 1, every in-text number, and a verification report.

```bash
cd analysis
conda env create -f environment.yml && conda activate cnt-analysis
export CNT_DATA=$PWD/../data_download/from_source   # the dir with tables/, data/
python run_all.py
```
Outputs land in `analysis/out/{figures,tables,reports}`, ending with
`out/reports/verification_report.md` — **50/50 checks PASS** against the
manuscript (in-text values + Table 1 gene sets).
See [`analysis/README.md`](analysis/README.md) for the stage→figure map.

## Role partitioning (a hard rule)
- **`data_download` fetches; `analysis` computes.** Every source download and
  its cache lives in `data_download`. `analysis` is deterministic and offline:
  given the pinned tables, fixed seeds (`analysis/00_config.py`) reproduce every
  value exactly. A *fresh* `data_download` run drifts only within the tolerances
  documented in `DOWNLOAD_PIPELINE.md` (aligner/annotation/release drift);
  structure and conclusions are unchanged.

## Layout
```
hackathon/
├── data_download/          # 18 source stages + run_all.R + DOWNLOAD_PIPELINE.md
│   ├── 00_config.R … 18_gtex_bulk.py
│   ├── source_manifest.csv # every source, URL, output
│   └── from_source/        # (generated) cached tables the analysis reads
└── analysis/               # 16 reproduction stages + run_all.py
    ├── 00_config.py         # input contract, pinned params + seeds
    ├── cnt_io.py            # cached loaders + value recording
    ├── _mechanism/_features/_amplicon/_depmap.py   # shared helpers
    ├── 10_… 31_…            # figure/value stages
    ├── 90_values_manifest.py, 91_verify_paper.py
    └── out/                 # (generated) figures, tables, reports
```

## Requirements
- **data_download**: R ≥ 4.4 (data.table, arrow, httr, jsonlite, org.Hs.eg.db)
  + Python ≥ 3.10 for the two `.py` stages (`cellxgene-census`, `pandas`,
  `pyarrow`, `requests`).
- **analysis**: the `cnt-analysis` conda env (`analysis/environment.yml`).
