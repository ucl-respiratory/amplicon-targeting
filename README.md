# CN-targeting: copy-number–driven multi-antigen co-target nomination

Code for the bioRxiv preprint _"From amplicon to antigen: a quantified
transmission map that nominates multi-antigen antibody–drug-conjugate co-target
sets across cancer types"_ (Lam, Walker-Samuel & Pennycuick, 2026).

Created as part of the Claude Science Hackathon, July 2026. We thank Anthropic and Glastone Institutes for their support.

Cancer amplicons carry many genes, most of them innocent passengers. This
codebase traces how a copy-number amplification reaches the cell surface —
copy number → mRNA → protein → single-cell co-detection — and uses that map to
nominate sets of co-amplified **passenger** surface antigens that a multispecific
antibody–drug conjugate could target together.

The work is in two stages: `data_download/` fetches and caches every input from
its primary source, and `integrated/` computes every figure, table and number in
the preprint from those cached tables.

## `data_download/` — assemble the source data

An R + Python pipeline that downloads and parses each input into
`from_source/{tables,data}/`: CPTAC proteogenomics (PDC + GDC), CORUM complexes,
DepMap 23Q4 (CRISPR dependency and Chronos gene-effect, copy number, expression),
the TCGA ATAC-seq atlas, UniProt membrane topology, HPA, GTEx v8, Ensembl gene
coordinates, gnomAD v4.1 constraint, STRING v12, the ADC Atlas, and CZ CELLxGENE
malignant-cell slices. Releases are pinned (CELLxGENE census `2025-11-08`, HPA
v24, GTEx v8, DepMap 23Q4, gnomAD v4.1) and every source is logged to
`source_manifest.csv`.

```bash
cd data_download
Rscript run_all.R          # download + assemble → from_source/
```

See [`data_download/DOWNLOAD_PIPELINE.md`](data_download/DOWNLOAD_PIPELINE.md)
for the stage-by-stage source map.

## `integrated/` — reproduce the preprint

A Python pipeline over the cached `from_source/` tables. Each stage produces part
of the paper:

| Stage                       | Produces                                                           |
| --------------------------- | ------------------------------------------------------------------ |
| `00a_transmissibility.py`   | transmissibility atlas (the predicted quantity)                    |
| `00d_gene_features.py`      | 39-feature gene-property table                                     |
| `01_transmission_gates.py`  | cascade + regulatory gates — Figure 1                              |
| `02_predictor.py`           | gene-property predictor + positional/lineage controls — Figure 2   |
| `03_empirical_bayes.py`     | measurement–prediction shrinkage                                   |
| `00b`/`00c`/`04c`           | recurrent amplicons, co-elevation, nomination funnel — Figure 3    |
| `04d_constructs.py`         | nominated antigens + multivalent constructs — Table 1, Figures 4–5 |
| `05_single_cell_andgate.py` | single-cell co-detection + normal-tissue burden — Figure 6         |
| `d1`/`d2`                   | Supplementary Tables S3, S5                                        |
| `verify_paper.py`           | checks every reported value against the computed artifacts         |

```bash
cd integrated
export CNT_DATA=$PWD/../data_download/from_source
python 00a_transmissibility.py && python 00d_gene_features.py
python 01_transmission_gates.py
python 02_predictor.py && python 03_empirical_bayes.py
python 00b_target_funnel.py && python 00c_amplicons.py
python 04c_pancancer_nomination.py && python 04d_constructs.py
python 05_single_cell_andgate.py
python d1_threshold_sensitivity.py && python d2_cross_cohort.py
python verify_paper.py
```

Outputs land in `integrated/{tables,figures,reports}`. Thresholds and seeds are
fixed in `config.py` (seed 2), so runs are deterministic. `04d` and `05` run a
2,000-resample donor-block bootstrap over the single-cell slices (~20 min each);
the rest are seconds-to-minutes.

The predictor uses 39 gene-property features and no protein-derived input;
`tables/feature_coverage_report.csv` lists each feature and its source. Most are
rebuilt live from source; ten protein-descriptor features (disorder,
secondary-structure, dN/dS, gene age, UTR lengths, codon optimality) are served
from a committed snapshot at
`data_download/from_source/data/raw/features/external_predicted_descriptors_snapshot.csv`
because their upstream servers are no longer reliably reachable.

## `biorxiv/` — the manuscript

LaTeX source and rendered PDF (`manuscript.tex`, `manuscript.pdf`), with figures
copied from `integrated/figures/`.

## Requirements

- **data_download**: R ≥ 4.4 (data.table, arrow, httr, jsonlite, org.Hs.eg.db)
  and Python ≥ 3.10 (`cellxgene-census`, pandas, pyarrow, requests) for the two
  Python stages.
- **integrated**: a scientific-Python stack (pandas, numpy, scipy, scikit-learn,
  matplotlib); `02_predictor.py` additionally needs XGBoost.
