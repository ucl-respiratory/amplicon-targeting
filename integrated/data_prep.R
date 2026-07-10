#!/usr/bin/env Rscript
# =============================================================================
# integrated/data_prep.R  --  Drive the from-source data_download for the
# integrated paper, running every LOAD-BEARING + REACHABLE layer on real data
# and caching all outputs under data_download/from_source (in the Dropbox repo).
#
# It does NOT modify any file in data_download/. It sources data_download's own
# stage scripts and config, runs the layers whose dependencies are installable
# in this environment, and provides structurally-valid EMPTY stubs for the two
# layers whose R packages will not solve against R 4.5:
#     - methylation (needs ChAMP; minor regulatory proxy, ATAC is the headline)
#     - somatic SNV (needs maftools; NOT consumed by any downstream analysis)
# The stub layers carry NA / zero rows, never invented values, and the gap is
# recorded in reports/data_prep_layers.csv.
#
# Resumable: GDC downloads are cached; parsed .RData are skip-if-exists here.
# =============================================================================
t0 <- Sys.time()
DD <- "/Users/adam/Dropbox/Projects/CN-targeting/hackathon/data_download"
setwd(DD)                       # FS_ROOT defaults to <cwd>/from_source (Dropbox cache)
options(timeout = 3600)
HERE <- DD                      # stage scripts live directly in data_download/

## from-scratch download flags (as run_all.R sets them) -----------------------
Sys.setenv(FS_GDC_RUN = "1", FS_AUX_RUN = "1", FS_CORUM_DEPMAP_RUN = "1",
           FS_ATAC_RUN = "1", FS_TOPOLOGY_RUN = "1")

source(file.path(HERE, "00_config.R"))
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

layers <- list()
mark <- function(name, status, note="") {
  layers[[length(layers)+1L]] <<- data.frame(layer=name, status=status, note=note)
  log(sprintf("LAYER %-14s %s %s", name, status, note))
}

skip_if <- function(path) file.exists(path) && file.info(path)$size > 0

## ---- Stage 01: annotation (genemap) ----------------------------------------
if (skip_if(INPUTS$genemap)) { mark("annotation","cached") } else {
  source(file.path(HERE,"01_annotation.R")); mark("annotation","built") }

## ---- Stage 02: proteome (PDC) — already parsed & cached in from_source ------
if (skip_if(INPUTS$proteome)) { mark("proteome","cached") } else {
  source(file.path(HERE,"02_proteome_pdc.R")); mark("proteome","built") }

## ---- Stage 03: GDC genomics — CN, ASCAT purity, RNA (real); meth/snv stub ---
# Source only the function DEFINITIONS (everything before the auto-run block),
# then call the load-bearing builders ourselves.
s03 <- readLines(file.path(HERE,"03_gdc_genomics.R"))
run_blk <- grep("Run all four layers", s03)
if (!length(run_blk)) run_blk <- length(s03) + 1L
eval(parse(text = paste(s03[seq_len(run_blk[1]-1L)], collapse = "\n")))

if (skip_if(INPUTS$cn))  { mark("copy_number","cached") } else {
  build_cn(); mark("copy_number","built","GDC AscatNGS gene-level") }
if (skip_if(INPUTS_AUX_PURITY)) { mark("ascat_purity","cached") } else {
  tryCatch({ build_ascat_purity(); mark("ascat_purity","built") },
           error=function(e) mark("ascat_purity","skip", conditionMessage(e))) }
if (skip_if(INPUTS$rna)) { mark("rna","cached") } else {
  build_rna(); mark("rna","built","GDC STAR-Counts + DESeq2 norm") }

# ---- meth/snv EMPTY stubs (deps unsolvable against R 4.5) -------------------
# Build structurally-valid empty objects matching what stage 09 expects, so the
# assembly runs and emits meth_table (all-NA) + snv_table (empty). No values are
# invented; downstream methylation features resolve to NA and are reported as a
# documented gap.
if (!skip_if(INPUTS$snv)) {
  snvdata <- data.frame(gene=character(0), caseid=character(0),
                        classification=character(0), stringsAsFactors=FALSE)
  snvpheno <- data.frame(caseid=character(0), stringsAsFactors=FALSE)
  save(snvdata, snvpheno, file = INPUTS$snv)
  mark("snv","stub","maftools unsolvable (R4.5); unused downstream")
} else mark("snv","cached")

if (!skip_if(INPUTS$meth_genes)) {
  # empty gene x case matrix; stage 09 melts it to an empty meth_table
  methdata.genes <- matrix(numeric(0), nrow=0, ncol=0)
  methpheno <- data.frame(caseid=character(0), stringsAsFactors=FALSE)
  save(methdata.genes, methpheno, file = INPUTS$meth_genes)
  methdata.tss <- matrix(numeric(0), nrow=0, ncol=0)
  save(methdata.tss, methpheno, file = INPUTS$meth_tss)
  mark("methylation","stub","ChAMP unsolvable (R4.5); minor regulatory proxy")
} else mark("methylation","cached")

## ---- Stage 04: auxiliary (GTEx, TCSA surface, HPA, DescribePROT, ADC Atlas) -
tryCatch({ source(file.path(HERE,"04_auxiliary.R")); mark("auxiliary","built") },
         error=function(e) mark("auxiliary","partial", conditionMessage(e)))

## ---- Stage 05: CORUM + DepMap ----------------------------------------------
tryCatch({ source(file.path(HERE,"05_corum_depmap.R")); mark("corum_depmap","built") },
         error=function(e) mark("corum_depmap","partial", conditionMessage(e)))

## ---- Stage 14: TCGA ATAC-seq accessibility (headline chromatin gate) --------
tryCatch({ source(file.path(HERE,"14_atac_gdc.R")); mark("atac","built") },
         error=function(e) mark("atac","partial", conditionMessage(e)))

## ---- Stage 15: UniProt membrane topology (surface gate) --------------------
tryCatch({ source(file.path(HERE,"15_uniprot_topology.R")); mark("topology","built") },
         error=function(e) mark("topology","partial", conditionMessage(e)))

## ---- Assembly chain 06-13 --------------------------------------------------
## Populate the shared input environment IN once, from the correct location
## (stage scripts guard on `if (!exists("IN"))` and otherwise try to source a
## nonexistent relative path "from_source/_load_inputs.R").
source(file.path(HERE, "_load_inputs.R"))
for (st in c("06_harmonize.R","07_proteome_relative.R","08_copynumber.R",
             "09_rna_meth_snv.R","10_assemble.R","11_annotate.R",
             "12_protein_core.R","13_feature_adapter.R")) {
  log("assembly:", st)
  source(file.path(HERE, st))
}
mark("assembly","built","06-13 -> omic tables + feature matrix")

## ---- record layer status ---------------------------------------------------
lay <- do.call(rbind, layers)
outdir <- "/Users/adam/Dropbox/Projects/CN-targeting/hackathon/integrated/reports"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(lay, file.path(outdir, "data_prep_layers.csv"), row.names = FALSE)
log("DONE in", round(difftime(Sys.time(), t0, units="mins"),1), "min")
print(lay)
