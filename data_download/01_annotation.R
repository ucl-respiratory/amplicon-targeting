# =============================================================================
# 01_annotation.R  --  Gene location/annotation map from org.Hs.eg.db
# -----------------------------------------------------------------------------
# Builds `genemap`: one entry per Entrez gene with chromosome, start, cytoband,
# chromosome arm, and Ensembl id. This is the annotation backbone the assembly
# chain joins onto (see 12_annotate.R, which COLLAPSES genemap to one row per
# symbol before joining, so the multi-Ensembl rows produced here are expected).
#
# SOURCE: Bioconductor annotation package org.Hs.eg.db (installed from
#   Bioconductor; version recorded to source_manifest.csv). No web download --
#   the annotation ships with the package, which is the canonical, versioned,
#   reproducible source for NCBI gene-location mappings.
#
# EXACT REPRODUCTION of the original download.data.R genemap build:
#   1. org.Hs.egCHRLOC  -> gene_id, start_location, Chromosome
#   2. merge org.Hs.egSYMBOL   (gene_id -> symbol)
#   3. merge org.Hs.egMAP      (gene_id -> cytogenetic_location)
#   4. keep autosomes only (Chromosome %in% 1:22; drops X/Y + alt contigs)
#   5. drop duplicated symbols (keep first)   <-- BEFORE the Ensembl merge
#   6. arm  = <chr-number><p|q>  parsed from the cytoband
#      POS  = abs(start_location)
#   7. merge org.Hs.egENSEMBL  (gene_id -> ensembl_id; MULTIPLIES some symbols)
#
# Emits: PARSED/genemap.RData  (object `genemap`; columns gene_id,
#        start_location, Chromosome, symbol, cytogenetic_location, arm, POS,
#        ensembl_id). Matches the cached genemap.RData schema exactly.
# =============================================================================

suppressWarnings(suppressMessages({
  library(org.Hs.eg.db)
  library(stringr)
}))
if (!exists("PARSED")) source(file.path("from_source", "00_config.R"))

build_genemap <- function() {
  m  <- data.frame(org.Hs.egCHRLOC)                      # gene_id, start_location, Chromosome
  m2 <- data.frame(org.Hs.egSYMBOL)                      # gene_id, symbol
  m  <- merge(m, m2, by = "gene_id")
  m2 <- data.frame(org.Hs.egMAP)                         # gene_id, cytogenetic_location
  m  <- merge(m, m2, by = "gene_id")
  m  <- m[which(m$Chromosome %in% as.character(1:22)), ] # autosomes only (drops alt contigs)
  m  <- m[which(!duplicated(m$symbol)), ]                # one row per symbol (pre-Ensembl)
  m$arm <- vapply(m$cytogenetic_location, function(l) {
    paste0(unlist(strsplit(l, "p|q"))[1], stringr::str_extract(l, "p|q"))
  }, character(1))
  m$POS <- abs(m$start_location)
  m2 <- data.frame(org.Hs.egENSEMBL)                     # gene_id, ensembl_id
  m  <- merge(m, m2, by = "gene_id")                     # re-inflates multi-Ensembl symbols
  rownames(m) <- NULL
  m
}

genemap <- build_genemap()

## ---- Provenance -------------------------------------------------------------
orgver <- as.character(utils::packageVersion("org.Hs.eg.db"))
save(genemap, file = INPUTS$genemap)
record_source(
  stage       = "01_annotation",
  source_name = paste0("org.Hs.eg.db ", orgver),
  url         = "Bioconductor::org.Hs.eg.db (CHRLOC/SYMBOL/MAP/ENSEMBL)",
  local_path  = INPUTS$genemap,
  note        = sprintf("genemap %d rows x %d cols; %d unique symbols; autosomes only",
                        nrow(genemap), ncol(genemap), length(unique(genemap$symbol)))
)

message("[annotation] genemap: ", nrow(genemap), " rows, ",
        length(unique(genemap$symbol)), " unique symbols, org.Hs.eg.db ", orgver)
