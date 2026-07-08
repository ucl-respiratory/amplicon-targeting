# =============================================================================
# 00_config.R  --  Configuration for the FROM-SOURCE CN-targeting pipeline
# -----------------------------------------------------------------------------
# This is a SELF-CONTAINED, publication-grade rebuild of the CN-targeting
# multi-omic ingestion pipeline. Unlike the companion `pipeline/` (which reads
# pre-downloaded .RData matrices cached in old/data/parsed/), this pipeline
# DOWNLOADS EVERYTHING FROM SOURCE REPOSITORIES and reproduces the same final
# analysis tables:
#     omic_table_annotated.parquet        (full-union multi-omic table)
#     omic_table_protein_core.parquet     (protein-scoped modelling table)
#     str_omic_table_rebuilt.csv          (James-compatible wide modelling CSV)
#     + CORUM complexes and DepMap 23Q4 cell-line matrices
#
# Data-flow overview (stage -> parsed object -> final table):
#   01_annotation      org.Hs.eg.db            -> genemap
#   02_proteome_pdc    PDC GraphQL API         -> cptac.data/pheno (+ matched)
#   03_gdc_genomics    GDC (GenomicDataCommons) -> cndata, rnadata, methdata, snvdata
#   04_auxiliary       GTEx / TCSA / HPA / DescribePROT / ADC Atlas
#   05_corum_depmap    CORUM 5.0 / DepMap 23Q4 (figshare)
#   06..13             assembly chain (harmonize -> ... -> annotate -> adapter)
#
# All source identifiers, URLs and versions are recorded here and re-emitted to
# source_manifest.csv by each stage. No hard-coded user home directories.
# =============================================================================

suppressWarnings(suppressMessages({
  library(data.table)
}))

## ---- Root locations ---------------------------------------------------------
# FS_ROOT   : root of this from-source pipeline's data tree (downloads + parsed)
# FS_OUT    : where final tables/figures/reports are written
# Override both via environment variables for a clean-machine run.
FS_ROOT  <- Sys.getenv("FS_ROOT",  unset = file.path(getwd(), "from_source"))
RAW      <- Sys.getenv("FS_RAW",   unset = file.path(FS_ROOT, "data", "raw"))
PARSED   <- Sys.getenv("FS_PARSED",unset = file.path(FS_ROOT, "data", "parsed"))
GDC_CACHE<- Sys.getenv("FS_GDC_CACHE", unset = file.path(FS_ROOT, "data", "gdc_cache"))

OUT      <- Sys.getenv("FS_OUT",   unset = file.path(FS_ROOT, "out"))
DIR_TAB  <- file.path(OUT, "tables")
DIR_FIG  <- file.path(OUT, "figures")
DIR_REP  <- file.path(OUT, "reports")
for (d in c(RAW, PARSED, GDC_CACHE, OUT, DIR_TAB, DIR_FIG, DIR_REP))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

## Path to the provenance manifest (appended to by every download stage).
SOURCE_MANIFEST <- file.path(DIR_REP, "source_manifest.csv")

## ---- Parsed object files (the CONTRACT the assembly chain consumes) ---------
# These names/shapes MUST match what the assembly stages (06..13) read. They are
# identical to the objects the original download.data.R produced.
#
# IMPORTANT: cptac.tmt.RData, cptac.matched.RData and drugs.RData each contain an
# object whose name COLLIDES with another file's object (both cptac.tmt and
# cptac.norm.qc.imp hold 'cptac.data'; both matched and proteome hold
# 'cptac.pheno'). The shared-IN loader (_load_inputs.R) must therefore NOT load
# them, or the pre-impute data.frame would clobber the imputed matrix. They are
# instead loaded explicitly, into their own environments, by the stages that use
# them (07_proteome_relative loads the TMT file; 13_feature_adapter loads matched
# + drugs). INPUTS below is exactly the set the original pipeline hoisted into IN.
INPUTS <- list(
  proteome    = file.path(PARSED, "cptac.norm.qc.imp.RData"), # cptac.data(matrix) + cptac.pheno
  cn          = file.path(PARSED, "cptac.cndata.RData"),      # cndata + cnpheno
  rna         = file.path(PARSED, "cptac.rnadata.RData"),     # rnadata, rnadata.norm + rnapheno
  meth_genes  = file.path(PARSED, "methdata.genes.RData"),    # methdata.genes + methpheno
  meth_tss    = file.path(PARSED, "methdata.tss.RData"),      # methdata.tss + methpheno
  snv         = file.path(PARSED, "cptac.snvdata.RData"),     # snvdata + snvpheno
  genemap     = file.path(PARSED, "genemap.RData"),           # genemap
  gtex        = file.path(PARSED, "gtex.tpms.RData"),         # gtex (gene x cancer-type)
  surface     = file.path(PARSED, "surface.genes.RData"),     # surface.genes
  secreted    = file.path(PARSED, "secreted.genes.RData")     # secreted.genes (+ .full)
)

# ASCAT tumour purity + ploidy (per case), pulled from GDC file-level fields on
# the allele-specific segment files. Loaded explicitly by 08_copynumber.R (not
# via IN), so it is optional: absent -> purity columns are NA, pipeline unchanged.
INPUTS_AUX_PURITY <- file.path(PARSED, "ascat_purity.RData")   # ascat_purity (caseid, tumor_purity, tumor_ploidy_ascat)

## Auxiliary parsed files loaded EXPLICITLY (not via IN) to avoid the name
## collisions noted above. The download stages still write these paths.
INPUTS_AUX <- list(
  proteome_tmt = file.path(PARSED, "cptac.tmt.RData"),        # cptac.data(pre-impute df)+pheno
  matched      = file.path(PARSED, "cptac.matched.RData"),    # cptac.matched.data/pheno
  drugs        = file.path(PARSED, "drugs.RData")             # drugs (TCSA S38)
)

## ---- Analysis parameters ----------------------------------------------------
# Cancer types carried through the proteome-core modelling table. HNSCC is
# downloaded but excluded from TUMOR_CODES because it has no GTEx tissue mapping
# and the original proteome-core scope excluded it; keep this list identical to
# the cached pipeline so downstream tables match.
TUMOR_CODES <- c("CCRCC", "GBM", "LSCC", "LUAD", "PDA", "UCEC")

# GTEx cancer-type -> tissue mapping (documented Lung collision: LUAD == LSCC).
GTEX_LUNG_COLLISION <- c("LUAD", "LSCC")

# Autosomes + sex chromosomes used for continuous ploidy.
PLOIDY_CHROMS <- c(as.character(1:22), "X", "Y")

# Minimum normal aliquots (per cancer type) to form a T-vs-N reference.
MIN_NORMAL_REF <- 3L

## ---- SOURCE IDENTIFIERS (verified live; see DOWNLOAD_PIPELINE.md) -----------

## PDC (Proteomic Data Commons) — CPTAC Discovery Proteome studies.
## GraphQL endpoint + the 7 Discovery Proteome pdc_study_ids (latest versions).
## Each maps to a TUMOR_CODE. HNSCC included for completeness.
PDC_GRAPHQL <- "https://pdc.cancer.gov/graphql"
PDC_STUDIES <- data.table(
  tumor_code   = c("UCEC",     "CCRCC",    "LUAD",     "GBM",      "HNSCC",    "LSCC",     "PDA"),
  pdc_study_id = c("PDC000125","PDC000127","PDC000153","PDC000204","PDC000221","PDC000234","PDC000270"),
  study_name   = c("CPTAC UCEC Discovery Study - Proteome",
                   "CPTAC CCRCC Discovery Study - Proteome",
                   "CPTAC LUAD Discovery Study - Proteome",
                   "CPTAC GBM Discovery Study - Proteome",
                   "CPTAC HNSCC Discovery Study - Proteome",
                   "CPTAC LSCC Discovery Study - Proteome",
                   "CPTAC PDA Discovery Study - Proteome")
)

## GDC (Genomic Data Commons) — CPTAC-3 program.
## RNA workflow modernized: original build used 'HTSeq - Counts' (RETIRED by
## GDC); current equivalent is 'STAR - Counts'. Documented drift in DOWNLOAD_PIPELINE.md.
GDC_PROGRAM      <- "CPTAC"
GDC_PROJECT      <- "CPTAC-3"
GDC_RNA_WORKFLOW <- Sys.getenv("FS_GDC_RNA_WORKFLOW", unset = "STAR - Counts")
# Optional GDC data-release pin (informational; recorded to manifest). Empty = latest.
GDC_RELEASE_PIN  <- Sys.getenv("FS_GDC_RELEASE", unset = "")

## GTEx v8 gene median TPM (public download).
## NOTE: GTEx migrated buckets. The legacy gtex_analysis_v8 bucket name contains
## underscores (breaks virtual-hosted HTTPS) and is deprecated; the current
## public path is the adult-gtex bucket, bulk-gex/v8/rna-seq/. Verified 200 OK,
## 56,200 genes x 54 tissues, all mapped tissues present.
GTEX_URL <- paste0("https://adult-gtex.storage.googleapis.com/",
                   "bulk-gex/v8/rna-seq/",
                   "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")
# GTEx tissue columns (as they appear in the .gct) mapped to cancer types.
GTEX_TISSUE_MAP <- c(CCRCC = "Kidney - Cortex",
                     GBM   = "Brain - Cortex",
                     LSCC  = "Lung",
                     LUAD  = "Lung",
                     PDA   = "Pancreas",
                     UCEC  = "Uterus")

## The Cancer Surfaceome Atlas (TCSA) — surface genes + druggability + drugs.
TCSA_S2  <- "http://fcgportal.org/TCSA/Download/Table%20S2.xlsx"   # surface genes
TCSA_S36 <- "http://fcgportal.org/TCSA/Download/Table%20S36.xlsx"  # druggability
TCSA_S38 <- "http://fcgportal.org/TCSA/Download/Table%20S38.xlsx"  # drugs

## Human Protein Atlas — predicted secreted proteins (tab-separated search export).
HPA_SECRETED_URL <- paste0("https://www.proteinatlas.org/api/search_download.php",
                           "?search=protein_class:Predicted%20secreted%20proteins",
                           "&format=tsv&columns=g,gs,eg,up,chr,chrp,pc,upbp",
                           "&compress=no")

## DescribePROT — per-protein disorder/structure features, human proteome (9606).
DESCRIBEPROT_URL <- "http://biomine.cs.vcu.edu/servers/DESCRIBEPROT/database_data_download/9606_value.csv"

## Human Protein Atlas bulk downloads (zipped TSVs) used by the off-target and
## single-cell analyses:
##   normal_tissue.tsv.zip     -> IHC protein levels across normal tissues/cells
##   rna_single_cell_type.tsv.zip -> normal single-cell RNA (nTPM) by cell type
## These are the canonical HPA "Downloadable data" assets, pinned to the v24
## release for reproducibility (the unversioned www path renamed normal_tissue
## to normal_ihc_data; the versioned host keeps the stable filename).
HPA_NORMAL_TISSUE_URL   <- Sys.getenv("FS_HPA_NORMAL_URL", unset =
  "https://v24.proteinatlas.org/download/tsv/normal_tissue.tsv.zip")
HPA_SINGLE_CELL_URL     <- Sys.getenv("FS_HPA_SC_URL", unset =
  "https://v24.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip")

## ADC Atlas — antibody-drug-conjugate target gene table (supplementary S3).
## Published static asset; recorded in manifest, path overridable.
ADC_ATLAS_URL <- Sys.getenv("FS_ADC_URL", unset =
  "https://static-content.springer.com/esm/art%3A10.1038%2Fs41417-023-00701-3/MediaObjects/41417_2023_701_MOESM2_ESM.xlsx")

## CORUM 5.x (Helmholtz Munich) — mammalian protein complexes.
## The old static /download/releases/current/*.zip paths are GONE (site is now a
## single-page app backed by a FastAPI service). Current download route:
##   GET {CORUM_API}/public/file/info                       -> list of file_ids
##   GET {CORUM_API}/public/file/download_current_file
##            ?file_id=human&file_format=txt                -> humanComplexes.txt
## NOTE: this host serves an INCOMPLETE TLS chain (missing intermediate CA), so
## httr/download.file may need ssl_verifypeer handling; see 05_corum_depmap.R.
## Verified: file_id=human,txt -> 7,812 rows, 28 columns (CORUM 5.x schema).
## The 5.x column names differ from the 2022 CORUM 4.x the analysis originally
## used; 05_corum_depmap.R maps them to the analysis-expected names.
CORUM_API      <- Sys.getenv("FS_CORUM_API", unset = "https://mips.helmholtz-muenchen.de/fastapi-corum")
CORUM_FILE_ID  <- "human"
CORUM_FORMAT   <- "txt"

## DepMap 23Q4 Public (figshare article 24667905). File IDs verified via
## figshare API. NOTE: 23Q4 uses NEW filenames vs the project's cached files;
## the rename map to the analysis-expected names lives in 05_corum_depmap.R.
DEPMAP_FIGSHARE_ARTICLE <- 24667905L
DEPMAP_FILES <- data.table(
  figshare_name = c("OmicsCNGene.csv",
                    "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                    "Model.csv"),
  figshare_url  = c("https://ndownloader.figshare.com/files/43346913",
                    "https://ndownloader.figshare.com/files/43347204",
                    "https://ndownloader.figshare.com/files/43746708"),
  analysis_name = c("Copy_Number_Public_23Q4.csv",
                    "Expression_Public_23Q4.csv",
                    "Model.csv")
)
# DepMap 23Q4 CRISPR essentiality files (same figshare article 24667905), written
# VERBATIM (no Entrez strip): the analysis reads the 'SYMBOL (EntrezID)' headers
# as-is (context-aware essentiality splits on ' (').
#   CRISPRGeneDependency.csv         -> Chronos dependency probability (models x genes), ~394 MB
#   AchillesCommonEssentialControls.csv -> common-essential gene control list
DEPMAP_CRISPR_FILES <- data.table(
  figshare_name = c("CRISPRGeneDependency.csv", "AchillesCommonEssentialControls.csv"),
  figshare_url  = c("https://ndownloader.figshare.com/files/43346574",
                    "https://ndownloader.figshare.com/files/43346361"),
  analysis_name = c("CRISPRGeneDependency.csv", "AchillesCommonEssentialControls.csv")
)
# DepMap proteomics (Gygi/CCLE) is NOT in the DepMap 23Q4 omics release; it comes
# from the CCLE proteomics quant table. Recorded here for the manifest.
DEPMAP_PROTEOMICS_URL <- Sys.getenv("FS_DEPMAP_PROT_URL", unset =
  "https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz")

## TCGA ATAC-seq atlas (Corces 2018, Science; hosted on the GDC ATAC-seq AWG
## publication page). Chromatin accessibility as the mediator of CN->mRNA
## transmission. The "All cancer type-specific count matrices in normalized
## counts [ZIP]" bundle contains 23 per-type <TYPE>_log2norm.txt matrices (hg38
## peaks). URL verified from https://gdc.cancer.gov/about-data/publications/ATACseq-AWG.
## Five of the six CPTAC types have matched tumour ATAC (PDA/pancreatic is absent).
ATAC_NORMCOUNTS_URL <- Sys.getenv("FS_ATAC_URL", unset =
  "https://api.gdc.cancer.gov/data/38b8f311-f3a4-4746-9829-b8e3edb9c157")
# ATAC cancer-type code -> our CPTAC tumor_code (LUSC=LSCC, KIRC=CCRCC; no PDA)
ATAC_TYPE_MAP <- c(LUAD = "LUAD", LUSC = "LSCC", KIRC = "CCRCC",
                   GBM = "GBM", UCEC = "UCEC")
ATAC_PROMOTER_WINDOW <- 2000L   # TSS +/- 2 kb promoter window for peak overlap

## ---- Small provenance helper ------------------------------------------------
# Append one provenance row (per downloaded file) to source_manifest.csv.
record_source <- function(stage, source_name, url, local_path,
                          note = "") {
  ok  <- file.exists(local_path)
  sz  <- if (ok) file.info(local_path)$size else NA_real_
  md5 <- if (ok) unname(tools::md5sum(local_path)) else NA_character_
  row <- data.table(
    stage       = stage,
    source_name = source_name,
    url         = url,
    local_file  = if (ok) normalizePath(local_path) else local_path,
    bytes       = sz,
    md5         = md5,
    access_date = as.character(Sys.Date()),
    note        = note
  )
  fwrite(row, SOURCE_MANIFEST,
         append = file.exists(SOURCE_MANIFEST))
  invisible(row)
}

## ---- Robust download helper -------------------------------------------------
# Skip-if-exists (resumable), mode 'wb' for binaries. Returns local path.
fetch <- function(url, dest, mode = "wb", overwrite = FALSE, quiet = TRUE, ...) {
  if (file.exists(dest) && !overwrite && file.info(dest)$size > 0) {
    message("[fetch] cached: ", basename(dest)); return(invisible(dest))
  }
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  message("[fetch] downloading: ", basename(dest))
  ok <- tryCatch({
    utils::download.file(url, destfile = dest, mode = mode, quiet = quiet, ...)
    TRUE
  }, error = function(e) { message("[fetch] FAILED: ", conditionMessage(e)); FALSE })
  if (!ok || !file.exists(dest)) stop("Download failed for ", url)
  invisible(dest)
}

norm_case <- function(x) gsub("\\.", "-", x)   # dot -> dash caseid canonicalization

message("[config] FS_ROOT = ", FS_ROOT)
message("[config] PARSED  = ", PARSED)
message("[config] OUT     = ", OUT)
