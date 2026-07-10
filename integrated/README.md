# integrated/ — "From amplicon to antigen"

A single, integrated pipeline for the combined paper. It unifies the two source
analyses (`analysis/` = measured co-elevation catalogue; `gene_intrinsic/` =
gene-property transmissibility predictor) into one narrative and one codebase:

> **Which co-amplified passenger genes reach the cell surface — measured where
> proteomics exists, predicted where it does not — and how to target them as
> AND-gated multi-antigen ADC sets.**

This folder writes *new* code for the whole pipeline. It reuses the two source
pipelines' **data contract** (it reads the same `data_download/from_source`
tables) but not their code paths, and it does not modify `analysis/`,
`gene_intrinsic/`, or `data_download/`.

## Paper flow (and the module that produces each part)

| # | Paper section | Module | Runs on |
|---|---|---|---|
| 1 | Clinical intro: multispecific ADCs for passenger targeting | (manuscript) | — |
| 2 | Transmission cascade & where it is gated | `01_transmission_gates.py` | data_download tables |
| 3 | Predicting transmissibility without protein data | `02_predictor.py` | gene-property features |
| 4 | Empirical-Bayes combine (measured ⊕ predicted) | `03_empirical_bayes.py` | atlas + predictor |
| 5 | Surface target identification (empirical **and** predicted) | `04_surface_targets.py` | transmissibility × UniProt topology × amplification |
| 6 | Single-cell AND-gate demonstration | `05_single_cell_andgate.py` | CELLxGENE malignant cells |
| 7 | Experimental plan | `06_experimental_plan.py` | (schematic) |
| 8 | Assemble the integrated manuscript (Word) | `07_manuscript.py` | computed values + figures |

## Data

All inputs come from the **from-source** download pipeline, cached in the
Dropbox repo at `data_download/from_source/` (visible, not hidden in a private
cache). Build it with:

```bash
cd data_download
Rscript ../integrated/data_prep.R      # runs every reachable + load-bearing layer
```

`data_prep.R` drives `data_download`'s own stage scripts, caches every download
under `from_source/` (resumable — GDC downloads and parsed `.RData` are
skip-if-exists), and records per-layer status in `reports/data_prep_layers.csv`.

**Two layers are carried as documented gaps**, because their R packages
(`ChAMP`, `maftools`) do not solve against this environment's R 4.5:
- **methylation** — a minor regulatory proxy; the headline chromatin gate is
  TCGA ATAC-seq accessibility, which builds on real data. Methylation features
  resolve to `NA` and are reported as absent, never imputed.
- **somatic SNV** — not consumed by any downstream analysis in either source
  pipeline (verified); carried as an empty layer.

## Configuration (`config.py`)

- **Amplification:** ploidy-adjusted CN ≥ 1.4 (the clinical-output basis of the
  `analysis/` pipeline), applied uniformly to both the measured and predicted arms.
- **Seed:** 2 (predictor + empirical-Bayes resamples).
- **Confidence tiers:** every nominated antigen is tagged `measured_high`,
  `measured_pred`, or `predicted_only` (see `CONF_TIERS`).
- Paths resolve `CNT_DATA` → `data_download/from_source` by default; override via
  `CNT_DATA` / `CNT_<NAME>` env vars.

## Outputs

- `figures/` — the paper figures (real, computed).
- `tables/` — the integrated surface-target table (empirical + predicted, tiered)
  and supporting tables.
- `reports/` — layer-status and value manifests for reproducibility.

## Run order

```bash
python 01_transmission_gates.py     # needs data_download CN+RNA+ATAC tables
python 02_predictor.py              # runs standalone (gene-property features)
python 03_empirical_bayes.py        # runs standalone (transmissibility atlas)
python 04_surface_targets.py        # UniProt topology fetched live if stage 15 absent
python 05_single_cell_andgate.py    # needs CELLxGENE slices (see below)
python 06_experimental_plan.py      # schematic, no data
python 07_manuscript.py             # assembles reports/integrated_manuscript.docx
```

Each module is standalone, reads `config.py`, and asserts its required inputs
exist (pointing back to `data_prep.R` if not). `07_manuscript.py` embeds
whichever figures exist and inserts a labelled placeholder for any not yet
produced, so the manuscript is always complete and re-runs cleanly.

### Single-cell layer (CELLxGENE census)

`05_single_cell_andgate.py` reads malignant-cell detection slices extracted from
the CELLxGENE census. Extract them with `run_census_slices.py` in a Python
environment carrying `cellxgene-census` + `tiledbsoma`:

```bash
python run_census_slices.py   # writes from_source/tables/cellxgene_{LUAD,LSCC}_malignant.parquet
```

TileDB's S3 VFS ignores `HTTP(S)_PROXY`; in a sandboxed network the proxy must
be passed explicitly via `vfs.s3.proxy_host` / `vfs.s3.proxy_port` (the script
reads it from `HTTPS_PROXY`). Normal single-cell reference is the HPA
`rna_single_cell_type.tsv` (fetched into `from_source/data/raw/`).

## Run state at first build

- Modules 02, 03, 05, 06, 07 run on real data now: the predictor reproduces
  leave-gene-out ρ = 0.54, R² = 0.32; empirical-Bayes gives the prior-only
  floor ρ = 0.53; the single-cell AND-gate shows 6–14× same-cell co-detection
  for four of five sets (the two-antigen LUAD 5p set is only 1.5×).
- Modules 01 and the genome-wide arm of 04 finalise when the copy-number layer
  (`data_prep.R`, a multi-GB GDC download) completes; the manuscript's Figure 1
  is a placeholder until then.
