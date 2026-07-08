#!/usr/bin/env Rscript
# =============================================================================
# run_all.R  --  Master runner for the from-source CN-targeting pipeline
# -----------------------------------------------------------------------------
# Runs the complete, self-contained pipeline that downloads every input from its
# canonical public repository and assembles the same final analysis tables the
# project uses. NOTHING here depends on the legacy cached .RData files.
#
# USAGE
#   Rscript from_source/run_all.R                 # full from-scratch run
#   Rscript from_source/run_all.R --from-cache    # skip downloads, assemble only
#                                                 # (reads existing PARSED/*.RData)
#   Rscript from_source/run_all.R --stages 06:13  # run a subset of stages
#
# ENVIRONMENT (all optional; see 00_config.R for defaults)
#   FS_ROOT        pipeline root (default ./from_source_run)
#   FS_RAW/FS_PARSED/FS_GDC_CACHE/FS_OUT   override individual dirs
#   FS_OLD_PLOIDY_CSV   legacy ploidy CSV for the old-vs-new comparison figure
#   FS_DEPMAP_ARM_URL   canonical source for DepMap arm-level CNAs (optional)
#
# The download stages (01-05) are guarded by their own FS_*_RUN flags; this
# runner sets them when doing a from-scratch run. Every stage is resumable:
# stages skip work whose output already exists (skip-if-exists), so an
# interrupted run can be restarted without redoing completed downloads.
#
# REQUIREMENTS: R >= 4.3 with data.table, arrow, jsonlite, httr, tidyverse,
# Hmisc, readxl, limma, org.Hs.eg.db, GenomicDataCommons, DESeq2, maftools,
# ChAMP, R.utils. See DOWNLOAD_PIPELINE.md for the full environment recipe.
# =============================================================================

t_start <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
from_cache <- "--from-cache" %in% args
HERE <- tryCatch(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
          error = function(e) "from_source")
if (is.na(HERE) || !nzchar(HERE)) HERE <- "from_source"

## stage selection (default all) --------------------------------------------
# --stages accepts "06:13" either glued (--stages=06:13) or as the next token
# (--stages 06:13).
stage_val <- NA_character_
si <- which(grepl("^--stages", args))
if (length(si)) {
  a <- args[si[1]]
  if (grepl("[:=]", sub("^--stages", "", a))) {
    stage_val <- sub("^--stages[= ]*", "", a)
  } else if (si[1] < length(args)) {
    stage_val <- args[si[1] + 1L]
  }
}
all_stages <- c(
  "01_annotation.R", "02_proteome_pdc.R", "03_gdc_genomics.R",
  "04_auxiliary.R",  "05_corum_depmap.R", "14_atac_gdc.R",
  "15_uniprot_topology.R", "16_cellxgene.py",
  "06_harmonize.R",  "07_proteome_relative.R", "08_copynumber.R",
  "09_rna_meth_snv.R","10_assemble.R", "11_annotate.R",
  "12_protein_core.R","13_feature_adapter.R")
# download stages (01-05 + ATAC atlas + UniProt topology + CELLxGENE) — skipped by --from-cache
download_stages <- c(all_stages[1:5], "14_atac_gdc.R",
                     "15_uniprot_topology.R", "16_cellxgene.py")
sel <- all_stages
if (!is.na(stage_val) && nzchar(stage_val)) {
  rng <- strsplit(stage_val, ":")[[1]]
  lo <- as.integer(rng[1]); hi <- as.integer(rng[length(rng)])
  sel <- all_stages[sapply(all_stages, function(s) {
    n <- as.integer(substr(s, 1, 2)); n >= lo && n <= hi })]
}
if (from_cache) sel <- setdiff(sel, download_stages)

## enable download stages for a from-scratch run ----------------------------
if (!from_cache) {
  Sys.setenv(FS_GDC_RUN = "1", FS_AUX_RUN = "1", FS_CORUM_DEPMAP_RUN = "1",
             FS_ATAC_RUN = "1", FS_TOPOLOGY_RUN = "1")
}

message("========================================================")
message(" CN-targeting from-source pipeline")
message(" mode: ", if (from_cache) "FROM-CACHE (assemble only)" else "FROM-SCRATCH (download + assemble)")
message(" stages: ", paste(sub("\\.R$", "", sel), collapse = ", "))
message("========================================================")

source(file.path(HERE, "00_config.R"))

run_log <- data.frame(stage = character(), status = character(),
                      seconds = numeric(), stringsAsFactors = FALSE)
for (s in sel) {
  message("\n>>> ", s)
  ts <- Sys.time()
  ok <- tryCatch({
    if (grepl("\\.py$", s)) {
      # Python stage (CELLxGENE). Pass OUT dir + census proxy via env. A non-zero
      # exit is tolerated (stage self-skips when census is unreachable) but logged.
      py <- Sys.getenv("FS_PYTHON", unset = "python")
      Sys.setenv(CNT_DATA_OUT = DIR_TAB)
      rc <- system2(py, shQuote(file.path(HERE, s)))
      if (rc != 0) message("  (python stage exit ", rc, " — Fig-6 slice optional)")
      "ok"
    } else { source(file.path(HERE, s)); "ok" }
  }, error = function(e) { message("  !! ERROR: ", conditionMessage(e)); "error" })
  dt <- as.numeric(difftime(Sys.time(), ts, units = "secs"))
  run_log <- rbind(run_log, data.frame(stage = s, status = ok, seconds = round(dt, 1)))
  if (ok == "error") stop("Stage ", s, " failed; see message above.")
}

## summary ------------------------------------------------------------------
write.csv(run_log, file.path(DIR_REP, "run_log.csv"), row.names = FALSE)
message("\n========================================================")
message(" DONE in ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1),
        " min")
message(" outputs in: ", OUT)
message(" final tables: omic_table_annotated.parquet, omic_table_protein_core.parquet,")
message("               str_omic_table_rebuilt.csv")
message(" provenance:   ", SOURCE_MANIFEST)
message("========================================================")
