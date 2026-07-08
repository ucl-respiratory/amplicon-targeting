# analysis

Reproduces **every figure and in-text value** of the CN-targeting manuscript
from the tables produced by [`../data_download`](../data_download). Run
`data_download` first, then point this pipeline at its output.

## Quick start
```bash
# 1. build the data (once) — see ../data_download
Rscript ../data_download/run_all.R

# 2. tell analysis where that output landed
export CNT_DATA=$PWD/../data_download/from_source     # dir containing tables/, data/

# 3. reproduce everything
conda env create -f environment.yml && conda activate cnt-analysis
python run_all.py
```
Outputs land in `analysis/out/` (`figures/`, `tables/`, `reports/`). The
verification report `out/reports/reproduction_report.md` compares every
regenerated number against the manuscript.

## What it produces
- **13 figures** (`out/figures/`) — the paper's Figures 1–13.
- **Table 1** (`out/tables/table1_cotarget_sets.csv`) — 7 tiered co-target sets.
- **paper_values.json** + **values_manifest.md** — every in-text number, keyed to
  the claim it supports.
- **reproduction_report.md** — regenerated-vs-paper comparison with pass/drift.

## Stages
| Stage | Reproduces | Fig |
|---|---|---|
| `10_transmission_reconciliation` | CN→mRNA vs CN→protein magnitude/ranking | 7 |
| `11_per_tissue_attenuation` | per-tissue dosage attenuation | 8 |
| `12_mechanism_decomposition` | responsiveness r=0.73 / 89% transcriptional | 9 |
| `13_coelevation_by_transmission` | co-elevation as a function of transmission | 10 |
| `14_regulatory_decomposition` | methylation + genomic variance | 11 |
| `15_atac_transmission` | ATAC promoter accessibility nested R² (0.034→0.092) | 12 |
| `16_statistical_robustness` | reliability / arm-CV / disattenuation | 13 |
| `20_tissue_amplicons` | recurrent amplicon landscape | 2 |
| `21_coelevation_nomination` | Fisher-FDR co-target nomination | — |
| `22_dependency_essentiality` | DepMap dependency + EGFR context reclassification | 3 |
| `23_accessibility_funnel` | 62→31→8→7 funnel + Table 1 | 1 |
| `24_multispecific_andgate` | AND-gate vs OR-gate normal burden | 4 |
| `25_offtarget_safety` | GTEx/HPA off-target screen | — |
| `30_singlecell_magnitude` | normal single-cell threshold collapse, 64× fold | 5 |
| `31_tumour_singlecell` | on-tumour same-cell co-expression + bootstrap CIs | 6 |
| `90_values_manifest` | aggregate all values | — |
| `91_verify` | compare to paper, write reproduction_report.md | — |

## Reproducibility boundary
Given the **pinned** data tables the manuscript used, this pipeline is
deterministic (fixed seeds in `00_config.py`) and reproduces every value exactly.
A *fresh* `data_download` run drifts within the tolerances documented in
`../data_download/DOWNLOAD_PIPELINE.md` (STAR-vs-HTSeq RNA, org.Hs.eg.db and
DepMap release drift); the pipeline structure and conclusions are unchanged.
`31_tumour_singlecell` depends on the CELLxGENE census slice from
`data_download` stage 16 (or re-pulls it); skip with `--skip-cellxgene`.
