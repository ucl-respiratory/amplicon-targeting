# CN-targeting: copy-number–driven multi-antigen co-target nomination

Reproducible code for the bioRxiv preprint *"From amplicon to antigen: a
quantified transmission map that nominates multi-antigen antibody–drug-conjugate
co-target sets across cancer types"* (Lam, Walker-Samuel & Pennycuick, 2026).

The pipeline maps how a copy-number amplification reaches the cell surface —
copy number → mRNA → protein → single-cell co-detection — and uses that map to
nominate multi-antigen ADC co-target sets built from **co-amplified passenger**
genes. Every figure, table and in-text number in the preprint is computed from
primary source data by the code here; nothing is hand-transcribed or cached from
an earlier run.

The repository is two pipelines you run in order:

```
data_download/   →   integrated/
(fetch + cache        (reproduce every figure, table
 all source data)      and in-text value + verifier)
```

## 1. `data_download/` — everything from source
An R + Python pipeline that downloads, parses and caches **every** input from
its primary source: CPTAC proteogenomics (PDC + GDC), CORUM complexes, DepMap
23Q4 (CRISPR dependency probability *and* Chronos gene-effect, plus CN and
expression), the TCGA ATAC-seq atlas, UniProt membrane topology, HPA, GTEx v8,
Ensembl gene coordinates, gnomAD v4.1 constraint, STRING v12, the ADC Atlas, and
CZ CELLxGENE malignant-cell slices. All network access lives here; releases are
pinned (GTEx v8, CELLxGENE census `2025-11-08`, HPA v24, DepMap 23Q4, gnomAD
v4.1) and every source is logged to `source_manifest.csv`.

```bash
cd data_download
Rscript run_all.R                 # from scratch: download + assemble
# outputs land in data_download/from_source/{tables,data,reports}
```
See [`data_download/DOWNLOAD_PIPELINE.md`](data_download/DOWNLOAD_PIPELINE.md)
for the stage-by-stage source map and reproducibility notes.

## 2. `integrated/` — reproduce the preprint
A Python pipeline that assumes `data_download` has run. It **reads only the
cached `data_download/from_source` tables** — no stage touches the network — and
computes **every** figure, table and in-text number: the transmissibility atlas
(`00a`), the gene-property feature table (`00d`), the cascade/gates figure
(`01`), the gene-property predictor and its leave-arm-out / leave-lineage-out
controls (`02`), empirical-Bayes shrinkage (`03`), the nomination funnel
(`00b`/`00c`/`04*`), the passenger-only driver-excluded construct (`04d`),
single-cell AND-gate co-detection (`05`), and Supplementary Tables S3 and S5
(`d1`/`d2`).

```bash
cd integrated
export CNT_DATA=$PWD/../data_download/from_source   # the dir with tables/, data/
python 00a_transmissibility.py && python 00d_gene_features.py   # atlas + features
python 01_transmission_gates.py                                 # cascade + gates (Fig 1)
python 02_predictor.py && python 03_empirical_bayes.py          # predictor + EB (Fig 2)
python 00b_target_funnel.py && python 00c_amplicons.py          # funnel + amplicons (Fig 3)
python 04c_pancancer_nomination.py && python 04d_constructs.py  # Table 1, constructs (Fig 4/5)
python 05_single_cell_andgate.py                                # single-cell AND-gate (Fig 6)
python d1_threshold_sensitivity.py && python d2_cross_cohort.py # Tables S3, S5
python verify_paper.py                                          # reproduction check
```
Environments: most stages run in the `cnt-census` conda env; `02_predictor.py`
needs XGBoost (the `james-ml` env). `04d`/`05` run a 2,000-resample donor-block
bootstrap over the single-cell slices and take ~20 min each; every other stage
is seconds-to-minutes. Thresholds and seeds are fixed in `integrated/config.py`
(seed 2 throughout), so outputs are deterministic.

`verify_paper.py` writes `reports/verification_report_integrated.md` —
**20/20 checks PASS** against the preprint, every checked value read from a
from-source pipeline artifact (none hand-transcribed). Outputs land in
`integrated/{tables,figures,reports}`.

### Predictor feature provenance
The predictor uses **39 gene-property features** (no protein-derived input);
provenance is itemised in `integrated/tables/feature_coverage_report.csv`:

- **29 features rebuilt live from source** (gnomAD v4.1, DepMap 23Q4 Chronos
  gene-effect + dependency, CORUM, UniProt sequence/topology/keywords, Ensembl
  mRNA, GTEx v8, STRING v12).
- **10 features carried as a pinned snapshot**
  (`data_download/from_source/data/raw/features/external_predicted_descriptors_snapshot.csv`)
  — protein disorder and secondary-structure predictions (vsl2/psipred/asaquick),
  dN/dS, gene age, UTR5/UTR3 lengths, codon optimality. Their source servers
  (DescribePROT) or bulk endpoints (Ensembl UTR/CDS) are no longer reliably
  fetchable, so the values are committed as a pinned external reference. They are
  gene-intrinsic annotations (no CPTAC/tumour data) and carry ~15% of predictor
  SHAP importance. This is the one documented provenance caveat; it is a *feature
  matrix* inheritance, not a result read from an old run — the transmissibility
  label, atlas, predictor, funnel, nomination and constructs are all end-to-end
  from source.

The pipeline reads nothing outside `data_download/from_source` at runtime. The
two predecessor projects it was fused from have been removed from the repo; the
full chain (`00a`→`00d`→`02`→`03`→`00b`/`00c`→`04*`→`05`→`verify_paper.py`) runs
and passes the verifier 20/20 with only `data_download/` and `integrated/`
present.

## Role partitioning (a hard rule)
- **`data_download` fetches; `integrated` computes.** Every source download and
  its cache lives in `data_download`. `integrated` is deterministic and offline:
  given the pinned tables, fixed seeds (`integrated/config.py`) reproduce every
  value exactly. A *fresh* `data_download` run drifts only within the tolerances
  documented in `DOWNLOAD_PIPELINE.md` (aligner/annotation/release drift);
  structure and conclusions are unchanged.

## Layout
```
hackathon/
├── data_download/          # source stages + run_all.R + DOWNLOAD_PIPELINE.md
│   ├── 00_config.R … 18_gtex_bulk.py
│   ├── source_manifest.csv # every source, URL, output
│   └── from_source/        # (generated) cached tables the pipeline reads
├── integrated/             # ★ the preprint pipeline — run this to reproduce
│   ├── config.py            # input contract, pinned params + seeds
│   ├── cnt_shared.py        # cached loaders (atlas/features)
│   ├── 00a_transmissibility.py, 00d_gene_features.py   # atlas + feature table
│   ├── 01_transmission_gates.py                        # cascade + gates (Fig 1)
│   ├── 02_predictor.py, 03_empirical_bayes.py          # predictor + EB (Fig 2)
│   ├── 00b_…04d_…, 05_…     # funnel, nomination, constructs, single-cell
│   ├── d1_threshold_sensitivity.py, d2_cross_cohort.py # Tables S3, S5
│   ├── 07_manuscript.py, verify_paper.py               # docx + reproduction check
│   └── tables/ figures/ reports/   # (generated) all from data_download/from_source
└── biorxiv/                # LaTeX source + rendered PDF of the preprint
    ├── manuscript.tex, manuscript.pdf
    └── Figures/            # figures copied from integrated/figures
```

## Requirements
- **data_download**: R ≥ 4.4 (data.table, arrow, httr, jsonlite, org.Hs.eg.db)
  + Python ≥ 3.10 for the `.py` stages (`cellxgene-census`, `pandas`,
  `pyarrow`, `requests`).
- **integrated**: the `cnt-census` conda env for most stages; `02_predictor.py`
  additionally needs XGBoost (the `james-ml` env). Both are standard
  scientific-Python stacks (pandas, numpy, scipy, scikit-learn, matplotlib).
