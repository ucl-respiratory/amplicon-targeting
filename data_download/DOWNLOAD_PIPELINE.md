# CN-targeting — from-source download & assembly pipeline

A **self-contained, publication-grade** rebuild of the CN-targeting multi-omic
data pipeline. Every input is downloaded from its canonical public repository;
**nothing** depends on the project's pre-downloaded / cached `.RData` files. The
pipeline ends with the same final analysis tables the project has been using.

This pipeline lives in `from_source/` and is numbered `01 → 13` with a master
runner `run_all.R`. It is independent of the legacy `pipeline/` (00–09) tree.

---

## 1. What it produces

| Output (in `OUT/tables/`) | Shape | Description |
|---|---|---|
| `omic_table_annotated.parquet` | 18,788,020 rows × 27 cols | gene × caseid tidy table, all layers + annotation |
| `omic_table_protein_core.parquet` | 4,976,680 rows | subset with a measured/imputed protein value |
| `str_omic_table_rebuilt.csv` | wide | James-compatible feature matrix (ML adapter) |
| `protein_relative.parquet`, `cn_table.parquet`, `rna_table.parquet`, `meth_table.parquet`, `snv_table.parquet` | tidy | per-layer intermediates |
| `describeprot_gene_features.parquet` | gene-level | structural / disorder features |

Plus 8 QC figures in `OUT/figures/`, provenance in `OUT/reports/source_manifest.csv`,
and a run log in `OUT/reports/run_log.csv`.

The `(gene, caseid)` pair is **unique** in both final tables (asserted with
`stopifnot`).

---

## 2. Pipeline stages

| # | Script | Source(s) | Emits |
|---|---|---|---|
| 01 | `01_annotation.R` | Bioconductor `org.Hs.eg.db` | `genemap` (gene → chr/arm/cytoband/ensembl) |
| 02 | `02_proteome_pdc.R` | **PDC GraphQL** (7 CPTAC Discovery Proteome studies) | `cptac.data` (TMT log-ratios), `cptac.pheno`, matched T–N |
| 03 | `03_gdc_genomics.R` | **GDC** (`GenomicDataCommons`) CPTAC-3 | CN (ASCAT), RNA (STAR-Counts), methylation (SeSAMe), SNV (MAF) |
| 04 | `04_auxiliary.R` | GTEx v8, TCSA, HPA, DescribePROT, ADC Atlas | `gtex`, `surface.genes`, `secreted.genes`, describeprot/adc raws |
| 05 | `05_corum_depmap.R` | CORUM 5.x, DepMap 23Q4 (figshare), CCLE | `corum_humanComplexes.txt`, DepMap cell-line matrices |
| 14 | `14_atac_gdc.R` | **TCGA ATAC-seq atlas** (Corces 2018, GDC) | `atac_gene_promoter_accessibility.parquet` (chromatin accessibility) |
| 06 | `06_harmonize.R` | (assembly) | `master_samples.csv` (aliquot → caseid/type/group) |
| 07 | `07_proteome_relative.R` | (assembly) | tumour-vs-normal relative protein + measured/ref flags |
| 08 | `08_copynumber.R` | (assembly) | tidy CN + **continuous ploidy** + ploidy-adjusted CN |
| 09 | `09_rna_meth_snv.R` | (assembly) | tidy RNA / methylation / SNV tables |
| 10 | `10_assemble.R` | (assembly) | `omic_table_clean.parquet` (join all layers) |
| 11 | `11_annotate.R` | (assembly) | `omic_table_annotated.parquet` (+ surface/secreted/gtex flags) |
| 12 | `12_protein_core.R` | (assembly) | `omic_table_protein_core.parquet` |
| 13 | `13_feature_adapter.R` | (assembly) | `str_omic_table_rebuilt.csv` + describeprot features |

`_load_inputs.R` is a helper (not a numbered stage) that hoists the parsed
objects into an `IN` environment for the assembly stages.

---

## 3. Running it

```bash
# full from-scratch run (downloads everything, then assembles)
Rscript from_source/run_all.R

# assemble only, from an existing set of parsed .RData (skips downloads)
Rscript from_source/run_all.R --from-cache

# run a contiguous subset of stages
Rscript from_source/run_all.R --stages 06:13
```

Every stage is **resumable** (skip-if-exists): re-running after an interruption
skips downloads/outputs that already exist. Configuration is via environment
variables (see `00_config.R`): `FS_ROOT`, `FS_RAW`, `FS_PARSED`, `FS_GDC_CACHE`,
`FS_OUT`, plus optional `FS_OLD_PLOIDY_CSV` and `FS_DEPMAP_ARM_URL`.

### Environment

R ≥ 4.3 with: `data.table`, `arrow`, `jsonlite`, `httr`, `tidyverse`, `Hmisc`,
`readxl`, `R.utils`; Bioconductor: `org.Hs.eg.db`, `limma`, `GenomicDataCommons`,
`DESeq2`, `maftools`, `ChAMP`. A ready conda recipe:

```bash
conda create -n cnt-source -c bioconda -c conda-forge \
  r-base r-data.table r-arrow r-jsonlite r-httr r-tidyverse r-hmisc \
  r-readxl r-r.utils bioconductor-org.hs.eg.db bioconductor-limma \
  bioconductor-genomicdatacommons bioconductor-deseq2 \
  bioconductor-maftools bioconductor-champ
```

---

## 4. Modernized endpoints & differences from the 2022 cached build

This rebuild uses **current** repository APIs. Where an original source was
retired or changed, the modern equivalent is used and the resulting numerical
drift is documented here. **The pipeline structure and final table shapes are
identical; individual values may differ within the tolerances below.**

### 4.1 Proteome — PDC GraphQL replaces hand-curated TMT + xlsx
The original parsed raw TMT `.tsv` files plus **seven** per-study, differently
formatted clinical/specimen `.xlsx` files by hand. This rebuild pulls the same
data programmatically from the **PDC GraphQL API**:
- `quantDataMatrix(data_type: "unshared_log2_ratio")` → the gene × aliquot
  matrix (the original used `Unshared.Log.Ratio`; `unshared_log2_ratio` is the
  current field name and matches to ~1e-3, e.g. A1BG −0.3517).
- `biospecimenPerStudy` → aliquot → case → `sample_type` (Primary Tumor / Solid
  Tissue Normal), replacing the hand-parsed specimen sheets.
- 7 CPTAC Discovery Proteome studies: PDC000125 (UCEC), PDC000127 (CCRCC),
  PDC000153 (LUAD), PDC000204 (GBM), PDC000221 (HNSCC), PDC000234 (LSCC),
  PDC000270 (PDA).

**Clinical columns:** the original `cptac.pheno` carried ~14 columns including
age/BMI/histology/smoking. None of these are read by any assembly stage or
appear in any final table, so the rebuild emits the 5 identity/group columns the
pipeline actually consumes. Full clinical annotation is available via PDC
`clinicalPerStudy` if ever needed.

### 4.2 RNA — GDC `STAR - Counts` replaces the **retired** `HTSeq - Counts`
GDC has **retired HTSeq**; the only current bulk RNA workflow for CPTAC-3 is
`STAR - Counts` (GENCODE v36). This rebuild reads the `unstranded` column of the
`augmented_star_gene_counts.tsv` file (the HTSeq-count analogue), skipping the 4
STAR summary rows (`N_unmapped`, `N_multimapping`, `N_noFeature`, `N_ambiguous`),
then ENSEMBL→symbol sums and DESeq2 size-factor normalizes as before.
**Drift:** different aligner + gene model → absolute counts and normalized values
**will not** match the 2022 HTSeq build to the digit. The CN→RNA dosage
relationship is preserved (median Spearman ρ = 0.16, 80% of genes positive).

### 4.3 Annotation — `org.Hs.eg.db` version drift
`genemap` is rebuilt from the installed `org.Hs.eg.db` (3.22.0 here vs the older
release used in 2022). Schema is identical (8 columns, autosomes only). Gene
coordinates shift slightly with the annotation release (e.g. A1BG start
58,345,182 vs 58,346,805). ~24.9k unique symbols (vs 24.6k), 23.2k overlapping.
The autosomes-only restriction means chrX/Y genes (e.g. TSPAN6) map to NA — this
is **faithful** to the original.

### 4.4 GTEx — identical release, **zero drift**
GTEx v8 median-TPM gct is pinned (`GTEx_Analysis_2017-06-05_v8_...gene_median_tpm`).
The 6-tissue → cancer-type collapse reproduces the cached `gtex` object with
**correlation 1.0000 on all 6 cancer types** and byte-identical values on shared
genes. **LUAD and LSCC both map to GTEx `Lung`** — an intentional, documented
identical-column collision.

### 4.5 CORUM 5.x — new API and schema
The old `/download/releases/current/*.zip` paths are gone; CORUM is now a SPA
backed by a FastAPI service at `mips.helmholtz-muenchen.de/fastapi-corum`. This
rebuild calls
`GET /public/file/download_current_file?file_id=human&file_format=txt`.
**The host serves an incomplete TLS chain**, so the fetch retries with peer
verification relaxed *for that host only* if the verified request fails. CORUM
5.x (7,866 human complexes) already exposes `complex_id` + `subunits_gene_name`
(semicolon-joined), matching the downstream `corum_data_cleaned.csv` contract.

### 4.6 DepMap 23Q4 — figshare rename **and gene-header reformat**
DepMap 23Q4 Public files come from figshare article **24667905**. Beyond a
filename rename, the file layout changed:
- `OmicsCNGene.csv` → `Copy_Number_Public_23Q4.csv`
- `OmicsExpressionProteinCodingGenesTPMLogp1.csv` → `Expression_Public_23Q4.csv`
- `Model.csv` → `Model.csv` (verbatim)

The cached CN/Expression files use **bare** gene-symbol headers (`FAM87B`,
`TSPAN6`); the 23Q4 files use `SYMBOL (EntrezID)` (`FAM87B (400728)`). Same genes,
same order, so `05_corum_depmap.R` strips the ` (EntrezID)` suffix to reproduce
the exact cached layout (cell lines as rows, genes as columns, readable with
`read.csv(row.names = 1)`). `Model.csv` is written **verbatim** (the analysis
reads it with plain `read.csv`, no `row.names`).
- **`Proteomics.csv`** is the CCLE/Gygi proteomics table (`GENE (UniProt)`
  headers), **not** part of the DepMap omics release — sourced separately.
- **`Arm-level_CNAs.csv`** is derived arm-level CN, downloaded only if
  `FS_DEPMAP_ARM_URL` is configured; it is used solely by empirical validation.

### 4.7 Chromatin accessibility — TCGA ATAC-seq atlas (new layer)
Stage 14 adds the chromatin-accessibility layer that CPTAC lacks, from the
pan-cancer TCGA ATAC-seq atlas (Corces 2018, Science). The GDC ATAC-seq AWG page
hosts "All cancer type-specific count matrices in normalized counts" as a single
ZIP (`api.gdc.cancer.gov/data/38b8f311-...`, ~630 MB) of 23 per-type
`<TYPE>_log2norm.txt` matrices (hg38 peaks). Five of the six CPTAC types are
matched — **LUAD, LUSC→LSCC, KIRC→CCRCC, GBM, UCEC** — while **PDA has no TCGA
ATAC**. Per gene, promoter accessibility = the maximum peak accessibility (mean
log2norm across that type's tumour samples) among peaks overlapping the TSS ±2 kb
window. TSS comes from the pipeline's own `genemap` (org.Hs.eg.db CHRLOC is
GRCh38: EGFR TSS 55,019,016 vs the hg38 reference 55,018,820 — within the 2 kb
window). Output: `atac_gene_promoter_accessibility.parquet`
(gene × {LUAD, LSCC, CCRCC, GBM, UCEC}).

This layer is the strongest predictor of CN→mRNA dosage transmission the project
has found: within each tumour type, promoter accessibility correlates positively
with per-gene transmission (LUAD r≈+0.31, LSCC r≈+0.33, CCRCC r≈+0.19, GBM
r≈+0.21, UCEC r≈+0.24), and adding it roughly triples the cross-validated R² for
predicting transmission over CPTAC-only features. **Caveats:** bulk tumour ATAC
(purity/cell-mixture confounding); ATAC and CPTAC are different cohorts of the
same cancer types (gene/type-level join, not sample-matched); promoter-only
(distal enhancers not yet included); PDA has no accessibility data.

---

## 5. Methodological fixes carried over from the clean rebuild

These corrections (introduced in the project's `pipeline/` rebuild) are retained:

- **Genemap multi-mapping** — symbols with >1 Ensembl id are collapsed to one row
  before the join; `(gene, caseid)` uniqueness re-asserted with `stopifnot`.
- **Continuous ploidy** — the old integer ploidy (np.bincount quantization) put
  129 tumours at an implausible ploidy = 1. Ploidy is now the median over
  chromosomes of per-`(caseid, chromosome)` median gene CN (correlates 0.85 with
  the old estimate where both exist). The old-vs-new comparison figure is
  optional (`FS_OLD_PLOIDY_CSV`).
- **Methylation replicate leak** — `make.names` replicate suffixes (`.1`, `.2`)
  that created phantom caseids are stripped; beta averaged per case.
- **No silent imputation** — `protein_measured`, `has_normal_ref`, `n_normal_ref`
  flags record what was actually measured vs upstream-imputed. `MIN_NORMAL_REF =
  3`. GBM has zero adjacent normals → `prot_relative` is left NA, not imputed.

---

## 6. Validation

`qa_dryrun.txt` records an end-to-end assembly run (`run_all.R --from-cache
--stages 06:12`) that reproduces the documented target tables **exactly** on
identical inputs — isolating assembly correctness from download drift:

- `omic_table_annotated`: 18,788,020 rows / 29,313 genes / 666 caseids ✓
- `omic_table_protein_core`: 4,976,680 rows / 7,495 genes / 664 caseids ✓
- `(gene, caseid)` uniqueness holds in both ✓

The download stages (01–05) are unit-tested individually against the live public
endpoints (recorded per stage). A full from-scratch run differs from the 2022
build only by the documented source drift above.
