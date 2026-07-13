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
**v9 dosage-transmission manuscript** (in-text values + Table 1 gene sets).
This verifier checks the v9 paper only; the integrated bioRxiv preprint has its
own verifier (see §3).
See [`analysis/README.md`](analysis/README.md) for the stage→figure map.

## 3. `integrated/` — the integrated bioRxiv preprint
A Python pipeline for the integrated manuscript *"From amplicon to antigen"*
(passenger-gene multispecific ADC targeting). Like `analysis/`, it reads only
the cached `data_download/from_source` tables and computes **every** figure,
table and in-text number — including the transmissibility atlas (`00a`), the
40-feature predictor table (`00d`), the gene-property predictor retrain (`02`),
empirical-Bayes shrinkage (`03`), the nomination funnel (`00b`/`00c`/`04*`),
single-cell AND-gate co-detection (`05`), and the passenger-only driver-excluded
construct (`04d`). The transmissibility label and the nomination pipeline are
rebuilt entirely from source; of the 40 predictor features, provenance is
three-way and fully itemised in `tables/feature_coverage_report.csv`:

- **29 features rebuilt live from source** (gnomAD v4.1, DepMap, CORUM, UniProt
  sequence/topology/keywords, Ensembl mRNA, GTEx v8, STRING v12). 24 of the 29
  reproduce the previous committed values at correlation/agreement ≥ 0.80; the
  five that lag (gnomAD pLI 0.80, two CORUM complex metrics, an aggregation
  proxy, the GO category) are all low-importance features.
- **10 features carried as a pinned snapshot copied from the previous
  `gene_intrinsic` feature table** (protein disorder and secondary-structure
  predictions vsl2/psipred/asaquick, dN/dS, gene age, UTR5/UTR3 lengths, codon
  optimality). These are *not* independently reproduced — their source servers
  (DescribePROT) or bulk endpoints (Ensembl UTR/CDS) are no longer reliably
  fetchable, so the committed values are pinned as an external reference. They
  are gene-intrinsic external annotations (no CPTAC/tumour data), and together
  carry ~15% of predictor SHAP importance (the disorder group alone ~12%).
- **1 feature (phylop_mean) is an all-null stub** — unavailable in both the
  committed table and this build; it contributes nothing to the predictor.

So the transmissibility label, atlas, predictor (including the cross-lineage
Kendall-W transfer metric, computed by leave-one-lineage-out refit in `02`),
funnel, nomination and constructs are end-to-end from source; the predictor is
trained on a feature matrix that still inherits 10 pinned columns from the
previous project's table (carried as a cached snapshot in `from_source`); this
is documented rather than hidden.

**The pipeline no longer reads the old `gene_intrinsic/` folder at runtime.** It
was verified by removing that folder and re-running the whole chain
(`00a`→`00d`→`02`→`03`→`00b`/`00c`→`04*`→`05`→`verify_paper.py`): every stage
completes and the verifier passes 20/20. The only residual reference is an
optional diagnostic in `00a`/`00d` that prints a reproduction correlation
against the previous atlas *if that folder is present*, and tolerates its
absence. `gene_intrinsic/` and the v9-only `analysis/` directory can therefore
be removed from the repository without affecting the integrated paper.

```bash
cd integrated
export CNT_DATA=$PWD/../data_download/from_source
python 00a_transmissibility.py && python 00d_gene_features.py   # atlas + features
python 02_predictor.py && python 03_empirical_bayes.py          # predictor + EB
python 00b_target_funnel.py && python 00c_amplicons.py          # funnel + amplicons
python 04c_pancancer_nomination.py && python 04d_constructs.py  # antigens + constructs
python 05_single_cell_andgate.py                                # single-cell AND-gate
python d1_threshold_sensitivity.py && python d2_cross_cohort.py # Tables S3, S5
python verify_paper.py                                          # reproduction check
```
`verify_paper.py` writes `reports/verification_report_integrated.md` —
**20/20 checks PASS** against the integrated preprint, every checked value read
from a from-source pipeline artifact (none hand-transcribed). The predictor
feature-matrix provenance caveat (10 pinned columns) is itemised above and in
`tables/feature_coverage_report.csv`.

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
├── analysis/               # 16 reproduction stages (v9 paper) + run_all.py
│   ├── 00_config.py         # input contract, pinned params + seeds
│   ├── cnt_io.py            # cached loaders + value recording
│   ├── _mechanism/_features/_amplicon/_depmap.py   # shared helpers
│   ├── 10_… 31_…            # figure/value stages
│   ├── 90_values_manifest.py, 91_verify_paper.py
│   └── out/                 # (generated) figures, tables, reports
└── integrated/             # integrated bioRxiv preprint pipeline
    ├── config.py            # input contract, pinned params + seeds
    ├── cnt_shared.py        # cached loaders (integrated atlas/features)
    ├── 00a_transmissibility.py, 00d_gene_features.py   # atlas + 40-feature table
    ├── 00b_…04d_…, 05_…     # funnel, nomination, constructs, single-cell
    ├── d1_threshold_sensitivity.py, d2_cross_cohort.py # Tables S3, S5
    ├── 07_manuscript.py, verify_paper.py               # docx + reproduction check
    └── tables/ figures/ reports/   # (generated) all from data_download/from_source
```

## Requirements
- **data_download**: R ≥ 4.4 (data.table, arrow, httr, jsonlite, org.Hs.eg.db)
  + Python ≥ 3.10 for the two `.py` stages (`cellxgene-census`, `pandas`,
  `pyarrow`, `requests`).
- **analysis**: the `cnt-analysis` conda env (`analysis/environment.yml`).
