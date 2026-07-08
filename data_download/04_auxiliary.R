# =============================================================================
# 04_auxiliary.R  --  GTEx, TCSA surfaceome, HPA secreted, DescribePROT, ADC Atlas
# -----------------------------------------------------------------------------
# Static reference downloads that annotate genes (they do not depend on the
# CPTAC cohort). Each is fetched from its canonical source and parsed to the
# object the assembly / adapter stages expect.
#
#   GTEx v8 gene median TPM   -> gtex (genes x 6 cancer-type columns)
#   TCSA Table S2 / S36 / S38 -> surface.genes (+ druggability), drugs
#   HPA predicted secreted    -> secreted.genes (+ .full)
#   DescribePROT 9606         -> describeprot feature table (used by 13_adapter)
#   ADC Atlas table S3        -> adc_atlas (used by 13_adapter)
#
# Emits:
#   PARSED/gtex.tpms.RData        gtex
#   PARSED/surface.genes.RData    surface.genes
#   PARSED/drugs.RData            drugs
#   PARSED/secreted.genes.RData   secreted.genes, secreted.genes.full
#   RAW/describeprot_9606_value.csv   (raw; parsed in 13_adapter)
#   RAW/adc_atlas_tableS3.xlsx        (raw; parsed in 13_adapter)
# =============================================================================

suppressWarnings(suppressMessages({
  library(data.table); library(readxl)
}))
if (!exists("PARSED")) source(file.path("from_source", "00_config.R"))

## ===========================================================================
## GTEx v8 gene median TPM  ->  gtex (genes x cancer types)
## ===========================================================================
# Original logic (download.data.R lines 794-815): read the .gct (skip 2 header
# lines), keep genes present in the proteome, mean-aggregate duplicate symbols,
# then select the 6 tissue columns and rename to cancer types. LUAD and LSCC
# BOTH map to "Lung" (documented identical-column collision).
build_gtex <- function(restrict_genes = NULL) {
  gz <- file.path(RAW, basename(GTEX_URL))
  fetch(GTEX_URL, gz)
  # .gct: line1 '#1.2', line2 'nrow ncol', line3 header (Name, Description, tissues...)
  gtex_raw <- fread(cmd = paste("gzip -dc", shQuote(gz)), skip = 2, header = TRUE,
                    sep = "\t", check.names = FALSE)
  desc_col <- "Description"
  if (!is.null(restrict_genes))
    gtex_raw <- gtex_raw[get(desc_col) %in% restrict_genes]
  genes <- gtex_raw[[desc_col]]
  tissue_cols <- setdiff(colnames(gtex_raw), c("Name", "Description"))
  mat <- as.matrix(gtex_raw[, ..tissue_cols])
  agg <- aggregate(mat, by = list(genes), FUN = mean)     # mean over dup symbols
  rownames(agg) <- agg$Group.1; agg$Group.1 <- NULL
  # Select tissue columns in TUMOR_CODES order (LUAD & LSCC both = Lung).
  sel <- GTEX_TISSUE_MAP[c("CCRCC","GBM","LSCC","LUAD","PDA","UCEC")]
  gtex <- agg[, sel]
  colnames(gtex) <- c("CCRCC","GBM","LSCC","LUAD","PDA","UCEC")
  save(gtex, file = INPUTS$gtex)
  record_source("04_auxiliary", "GTEx v8 gene median TPM", GTEX_URL, gz,
                note = sprintf("gtex %d genes x 6 cancer types (LUAD==LSCC==Lung)",
                               nrow(gtex)))
  message("[gtex] gtex: ", nrow(gtex), " genes x 6 cancer types; LUAD==LSCC: ",
          isTRUE(all.equal(gtex$LUAD, gtex$LSCC)))
  invisible(gtex)
}

## ===========================================================================
## TCSA surfaceome (Table S2 surface genes + S36 druggability + S38 drugs)
## ===========================================================================
# Original: read.xls(skip=1); first data row holds the real column names for a
# block of columns; drop the two header rows; merge druggability on
# ENSEMBL.Gene.ID. Reproduced with readxl.
build_surface <- function() {
  s2 <- file.path(RAW, "tcsa_surface_S2.xlsx")
  s36 <- file.path(RAW, "tcsa_druggability_S36.xlsx")
  s38 <- file.path(RAW, "tcsa_drugs_S38.xlsx")
  fetch(TCSA_S2, s2); fetch(TCSA_S36, s36); fetch(TCSA_S38, s38)

  sg <- as.data.frame(read_excel(s2, skip = 1))
  # Row 1 carries sub-headers for a block of columns (as in the original).
  sel <- which(as.character(unlist(sg[1, ])) != "" & !is.na(unlist(sg[1, ])))
  colnames(sg)[sel] <- as.character(unlist(sg[1, sel]))
  sg <- sg[-c(1, 2), ]
  drug <- as.data.frame(read_excel(s36, skip = 1))
  surface.genes <- merge(sg, drug, by = "ENSEMBL.Gene.ID", all.x = TRUE)
  drugs <- as.data.frame(read_excel(s38, skip = 1))

  save(surface.genes, file = INPUTS$surface)
  save(drugs,         file = INPUTS_AUX$drugs)
  record_source("04_auxiliary", "TCSA Table S2 surface genes", TCSA_S2, s2,
                note = sprintf("surface.genes %d rows x %d cols", nrow(surface.genes), ncol(surface.genes)))
  record_source("04_auxiliary", "TCSA Table S36 druggability", TCSA_S36, s36)
  record_source("04_auxiliary", "TCSA Table S38 drugs", TCSA_S38, s38,
                note = sprintf("drugs %d rows", nrow(drugs)))
  message("[surface] surface.genes: ", nrow(surface.genes), " rows; drugs: ", nrow(drugs), " rows")
}

## ===========================================================================
## HPA predicted-secreted proteins  ->  secreted.genes(+ .full)
## ===========================================================================
# Original read a manually-downloaded proteinatlas TSV; here we hit the HPA
# search_download API for the 'Predicted secreted proteins' protein_class.
build_secreted <- function() {
  tsv <- file.path(RAW, "hpa_predicted_secreted.tsv")
  fetch(HPA_SECRETED_URL, tsv, mode = "w")
  secreted.genes.full <- fread(tsv, sep = "\t", header = TRUE, data.table = FALSE)
  gene_col <- if ("Gene" %in% colnames(secreted.genes.full)) "Gene" else colnames(secreted.genes.full)[1]
  secreted.genes <- secreted.genes.full[[gene_col]]
  save(secreted.genes, secreted.genes.full, file = INPUTS$secreted)
  record_source("04_auxiliary", "HPA predicted secreted proteins", HPA_SECRETED_URL, tsv,
                note = sprintf("%d secreted genes", length(secreted.genes)))
  message("[secreted] secreted.genes: ", length(secreted.genes), " genes")
}

## ===========================================================================
## DescribePROT + ADC Atlas (raw downloads; parsed downstream by 13_adapter)
## ===========================================================================
build_describeprot <- function() {
  dp <- file.path(RAW, "describeprot_9606_value.csv")
  ok <- tryCatch({ fetch(DESCRIBEPROT_URL, dp, mode = "w"); TRUE },
                 error = function(e) { message("[describeprot] download failed: ",
                                                conditionMessage(e)); FALSE })
  if (ok) record_source("04_auxiliary", "DescribePROT human proteome (9606)",
                        DESCRIBEPROT_URL, dp)
}
build_adc <- function() {
  adc <- file.path(RAW, "adc_atlas_tableS3.xlsx")
  ok <- tryCatch({ fetch(ADC_ATLAS_URL, adc); TRUE },
                 error = function(e) { message("[adc] download failed: ",
                                                conditionMessage(e)); FALSE })
  if (ok) record_source("04_auxiliary", "ADC Atlas supplementary Table S3",
                        ADC_ATLAS_URL, adc)
}

## ---- Run (GTEx restricted to proteome genes when available) -----------------
if (identical(Sys.getenv("FS_AUX_RUN"), "1") || !interactive()) {
  restrict <- NULL
  if (file.exists(INPUTS$proteome)) {
    e <- new.env(); load(INPUTS$proteome, e); restrict <- rownames(e$cptac.data)
  }
  build_gtex(restrict_genes = restrict)
  build_surface(); build_secreted(); build_describeprot(); build_adc()
} else {
  message("[aux] functions defined. Set FS_AUX_RUN=1 to execute all downloads.")
}
