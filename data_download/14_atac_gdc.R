# =============================================================================
# 14_atac_gdc.R  --  TCGA tumour chromatin accessibility (Corces 2018 ATAC atlas)
# -----------------------------------------------------------------------------
# Downloads the pan-cancer TCGA ATAC-seq atlas cancer-type-specific normalized
# count matrices from the GDC ATAC-seq AWG publication page, and computes a
# per-gene promoter-accessibility value for each of the CPTAC cancer types that
# have matched tumour ATAC data. This is the chromatin layer CPTAC lacks; it is
# the strongest predictor of CN->mRNA dosage transmission (see the project memos
# chromatin_accessibility_memo.md / atac_transmission_integration.md).
#
# SOURCE
#   Corces et al. 2018, Science 362:eaav1898, "The chromatin accessibility
#   landscape of primary human cancers." 410 tumours, 23 TCGA types.
#   GDC page: https://gdc.cancer.gov/about-data/publications/ATACseq-AWG
#   "All cancer type-specific count matrices in normalized counts [ZIP]"
#   -> api.gdc.cancer.gov/data/38b8f311-... (23 files <TYPE>_log2norm.txt, hg38).
#
# CANCER-TYPE MATCHING (5 of our 6 CPTAC types; PDA has no TCGA ATAC):
#   LUAD->LUAD, LUSC->LSCC, KIRC->CCRCC, GBM->GBM, UCEC->UCEC.
#
# METHOD (reproduces the project's ATAC analysis):
#   * peak accessibility = mean log2norm insertion count across that type's
#     tumour samples (columns after seqnames/start/end/name/score).
#   * per gene, promoter accessibility = MAX peak accessibility among peaks
#     overlapping the TSS +/- ATAC_PROMOTER_WINDOW (2 kb) window.
#   * TSS comes from the pipeline's own genemap (org.Hs.eg.db CHRLOC is GRCh38:
#     EGFR TSS 55,019,016 vs hg38 reference 55,018,820 -- within the 2 kb window,
#     so no separate coordinate source is needed).
#
# Peaks are hg38, matching the GRCh38 gene coordinates used throughout the
# pipeline. ~6,000-6,200 genes get a promoter-accessibility value per type.
#
# Emits:
#   PARSED/atac_gene_promoter_accessibility.parquet   (gene x {LUAD,LSCC,CCRCC,GBM,UCEC})
#   RAW/atac/<TYPE>_log2norm.txt                       (extracted source matrices)
#   OUT/figures/atac_transmission.png                  (accessibility vs CN->mRNA, if
#                                                        omic_table_annotated.parquet exists)
#
# NOTE: joined at the gene/cancer-type level (ATAC and CPTAC are different
# cohorts of the same cancer types, not sample-matched). Bulk tumour ATAC =>
# purity/cell-mixture confounding. Promoter-only (distal enhancers not included).
# =============================================================================

suppressWarnings(suppressMessages({
  library(data.table); library(arrow)
}))
if (!exists("PARSED")) source(file.path("from_source", "00_config.R"))
if (!exists("genemap")) {
  if (file.exists(INPUTS$genemap)) load(INPUTS$genemap) else
    stop("genemap not found; run 01_annotation.R first")
}

ATAC_DIR <- file.path(RAW, "atac")
dir.create(ATAC_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- Download + extract the per-type normalized-count matrices --------------
atac_download <- function() {
  needed <- paste0(names(ATAC_TYPE_MAP), "_log2norm.txt")
  have   <- file.exists(file.path(ATAC_DIR, needed))
  if (all(have)) { message("[atac] all ", length(needed), " matrices cached"); return(invisible()) }
  zip <- file.path(ATAC_DIR, "tcga_atac_normcounts.zip")
  fetch(ATAC_NORMCOUNTS_URL, zip)                       # ~630 MB ZIP of 23 matrices
  # extract only the 5 matrices we need
  all_in_zip <- utils::unzip(zip, list = TRUE)$Name
  target <- intersect(needed, basename(all_in_zip))
  utils::unzip(zip, files = target, exdir = ATAC_DIR, junkpaths = TRUE)
  record_source("14_atac_gdc",
                "TCGA ATAC-seq atlas normalized counts (Corces 2018)",
                ATAC_NORMCOUNTS_URL, zip,
                note = paste("extracted:", paste(target, collapse = ", ")))
  message("[atac] extracted: ", paste(target, collapse = ", "))
}

## ---- Per-gene promoter accessibility for one cancer type --------------------
# Sweep-line overlap of peaks against TSS+/-window, MAX peak accessibility.
gene_promoter_accessibility <- function(path, prom) {
  df <- fread(path, sep = "\t", header = TRUE)
  samp <- setdiff(colnames(df), c("seqnames", "start", "end", "name", "score"))
  df[, peak_acc := rowMeans(.SD, na.rm = TRUE), .SDcols = samp]   # mean across tumour samples
  peaks <- df[, .(chrkey = seqnames, pstart = start, pend = end, peak_acc)]
  out <- new.env(parent = emptyenv())
  for (chrk in unique(prom$chrkey)) {
    pk <- peaks[chrkey == chrk]
    if (!nrow(pk)) next
    gp <- prom[chrkey == chrk]
    setorder(pk, pstart)
    ps <- pk$pstart; pe <- pk$pend; pa <- pk$peak_acc
    for (i in seq_len(nrow(gp))) {
      gs <- gp$p_start[i]; ge <- gp$p_end[i]
      lo <- findInterval(ge, ps)                 # peaks with pstart <= ge
      if (lo == 0) next
      cand <- which(pe[seq_len(lo)] >= gs)       # and pend >= gs  => overlap
      if (length(cand)) assign(gp$gene[i], max(pa[cand]), envir = out)
    }
  }
  g <- ls(out)
  setNames(vapply(g, function(k) get(k, out), numeric(1)), g)
}

build_atac <- function() {
  atac_download()
  # promoter windows from genemap TSS (autosomes; genemap is autosome-only)
  prom <- as.data.table(genemap)[!is.na(start_location) & !is.na(Chromosome)]
  prom[, chrkey := paste0("chr", Chromosome)]
  prom[, `:=`(gene = symbol,
              p_start = abs(start_location) - ATAC_PROMOTER_WINDOW,
              p_end   = abs(start_location) + ATAC_PROMOTER_WINDOW)]
  prom <- prom[, .(gene, chrkey, p_start, p_end)]

  acc <- NULL
  for (atac_ct in names(ATAC_TYPE_MAP)) {
    our_ct <- ATAC_TYPE_MAP[[atac_ct]]
    f <- file.path(ATAC_DIR, paste0(atac_ct, "_log2norm.txt"))
    if (!file.exists(f)) { message("[atac] missing ", f, "; skipping ", our_ct); next }
    s <- gene_promoter_accessibility(f, prom)
    dt <- data.table(gene = names(s), acc = as.numeric(s)); setnames(dt, "acc", our_ct)
    acc <- if (is.null(acc)) dt else merge(acc, dt, by = "gene", all = TRUE)
    message(sprintf("[atac] %s (%s): %d genes with promoter accessibility",
                    our_ct, atac_ct, length(s)))
  }
  atac_gene_promoter_accessibility <- acc
  outf <- file.path(PARSED, "atac_gene_promoter_accessibility.parquet")
  write_parquet(atac_gene_promoter_accessibility, outf)
  record_source("14_atac_gdc", "Per-gene promoter accessibility (derived)",
                ATAC_NORMCOUNTS_URL, outf,
                note = sprintf("%d genes x %d cancer types (TSS+/-%dbp, max peak)",
                               nrow(acc), ncol(acc) - 1L, ATAC_PROMOTER_WINDOW))
  message("[atac] atac_gene_promoter_accessibility.parquet: ",
          nrow(acc), " genes x ", ncol(acc) - 1L, " cancer types")

  ## Optional QC figure: accessibility vs CN->mRNA transmission per type.
  omic_f <- file.path(DIR_TAB, "omic_table_annotated.parquet")
  if (file.exists(omic_f)) {
    suppressWarnings(suppressMessages(library(ggplot2)))
    om <- as.data.table(read_parquet(omic_f))[!is.na(cn) & !is.na(rna_log2)]
    # per-gene, per-type CN->mRNA transmission (Spearman) then join accessibility
    trans <- om[, .(transmission = if (.N >= 20)
                      suppressWarnings(cor(cn, rna_log2, method = "spearman")) else NA_real_),
                by = .(gene, tumor_code)]
    long <- melt(acc, id.vars = "gene", variable.name = "tumor_code",
                 value.name = "promoter_acc", na.rm = TRUE)
    m <- merge(trans[!is.na(transmission)], long, by = c("gene", "tumor_code"))
    if (nrow(m) > 100) {
      lab <- m[, .(r = cor(promoter_acc, transmission, use = "complete.obs"),
                   n = .N), by = tumor_code]
      p <- ggplot(m, aes(promoter_acc, transmission)) +
        geom_point(alpha = 0.1, size = 0.5, colour = "#4C72B0") +
        geom_smooth(method = "lm", se = FALSE, colour = "#C44E52") +
        facet_wrap(~ tumor_code, scales = "free_x") +
        labs(title = "Promoter chromatin accessibility vs CN->mRNA transmission",
             subtitle = "TCGA tumour ATAC (Corces 2018) joined to CPTAC per-gene dosage transmission",
             x = "promoter accessibility (max peak, mean log2norm)",
             y = "CN->mRNA transmission (Spearman rho)") +
        theme_minimal(base_size = 11)
      ggsave(file.path(DIR_FIG, "atac_transmission.png"), p, width = 9, height = 6, dpi = 200)
      message("[atac] accessibility~transmission correlation per type:")
      for (i in seq_len(nrow(lab)))
        message(sprintf("        %-6s r=%+.2f (n=%d)", lab$tumor_code[i], lab$r[i], lab$n[i]))
    }
  }
  invisible(acc)
}

## ---- Run --------------------------------------------------------------------
if (identical(Sys.getenv("FS_ATAC_RUN"), "1") || !interactive()) {
  build_atac()
} else {
  message("[atac] functions defined. Set FS_ATAC_RUN=1 to download + build.")
}
