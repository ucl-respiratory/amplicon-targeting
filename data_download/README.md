# data_download

Self-contained, from-source download & assembly pipeline for the CN-targeting
multi-omic analysis. Downloads every input from its canonical public repository
(PDC, GDC, GTEx, TCSA, HPA, CORUM, DepMap, UniProt, CELLxGENE) and produces the
analysis tables the downstream `analysis/` pipeline consumes. No
cached/pre-downloaded data required.

**Full documentation:** [DOWNLOAD_PIPELINE.md](DOWNLOAD_PIPELINE.md)

## Quick start
```bash
Rscript run_all.R                 # full from-scratch run
Rscript run_all.R --from-cache    # assemble only, from existing parsed .RData
Rscript run_all.R --stages 06:13  # run a contiguous subset
```

## Layout
- `00_config.R` — paths, endpoints, constants (override via `FS_*` env vars)
- `01`–`05`, `14`–`16` — download stages (annotation, proteome, GDC genomics, auxiliary, corum+depmap, TCGA ATAC-seq accessibility, UniProt membrane topology, CELLxGENE malignant-cell scRNA)
- `06`–`13` — assembly chain → `omic_table_annotated.parquet`, `omic_table_protein_core.parquet`, `str_omic_table_rebuilt.csv`
- stage `16_cellxgene.py` is the one **Python** stage (census API is Python-only); `run_all.R` shells out to it. Everything else is R.
- `_load_inputs.R` — helper (loads parsed objects for the assembly stages)
- `run_all.R` — master runner
- `source_manifest.csv` — provenance: every download source + URL
- `qa_dryrun.txt` + `qa_figures/` — validation record

Outputs land under `$FS_ROOT/out/` (git-ignored). See DOWNLOAD_PIPELINE.md for
the environment recipe and the documented differences from the 2022 cached build.
