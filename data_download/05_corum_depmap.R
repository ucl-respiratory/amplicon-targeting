# =============================================================================
# 05_corum_depmap.R  --  CORUM protein complexes + DepMap 23Q4 cell-line data
# -----------------------------------------------------------------------------
# CORUM (feature/validation code keys complexes by complex_id and expands
# subunits_gene_name) and the DepMap 23Q4 Public cell-line matrices used by the
# empirical-validation script (empirical_selection_depmap.R).
#
# CORUM 5.x  (Helmholtz Munich FastAPI):
#   GET {CORUM_API}/public/file/download_current_file?file_id=human&file_format=txt
#   -> humanComplexes.txt. Columns complex_id + subunits_gene_name already match
#   the analysis' corum_data_cleaned.csv contract (semicolon-joined subunits).
#   NOTE: host serves an incomplete TLS chain; we retry with peer verification
#   relaxed ONLY for this host if the verified request fails.
#
# DepMap 23Q4 Public  (figshare article 24667905):
#   OmicsCNGene.csv                               -> Copy_Number_Public_23Q4.csv
#   OmicsExpressionProteinCodingGenesTPMLogp1.csv -> Expression_Public_23Q4.csv
#   Model.csv                                     -> Model.csv
#   The analysis-cached CN/Expression files use BARE gene-symbol column headers,
#   whereas 23Q4 uses 'SYMBOL (EntrezID)'. We strip ' (EntrezID)' from the gene
#   columns so the reformatted files match the cached layout exactly (verified:
#   identical symbol order FAM87B,LINC01128,... / TSPAN6,TNMD,DPM1,...).
#   Cell lines are rows; genes are columns (as the analysis reads via row.names=1).
#
# DepMap proteomics + arm-level CNAs are NOT in the 23Q4 omics figshare release:
#   Proteomics.csv     <- CCLE proteomics (Nusinow/Gygi) protein_quant table.
#   Arm-level_CNAs.csv <- derived arm-level calls (documented; provide URL if
#                         a canonical source is configured, else skip w/ warning).
#
# Emits (to RAW/depmap/, names matching the analysis contract):
#   RAW/depmap/Copy_Number_Public_23Q4.csv
#   RAW/depmap/Expression_Public_23Q4.csv
#   RAW/depmap/Model.csv
#   RAW/depmap/Proteomics.csv            (if CCLE proteomics URL reachable)
#   RAW/depmap/Arm-level_CNAs.csv        (if source configured)
#   PARSED/corum_humanComplexes.txt      (CORUM 5.x)
# =============================================================================

suppressWarnings(suppressMessages({
  library(data.table); library(httr)
}))
if (!exists("PARSED")) source(file.path("from_source", "00_config.R"))
DEPMAP_DIR <- file.path(RAW, "depmap")
dir.create(DEPMAP_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- CORUM download (handles the incomplete TLS chain) ----------------------
corum_fetch <- function() {
  dest <- file.path(PARSED, "corum_humanComplexes.txt")
  if (file.exists(dest) && file.info(dest)$size > 0) {
    message("[corum] cached: ", basename(dest)); return(invisible(dest))
  }
  url <- sprintf("%s/public/file/download_current_file?file_id=%s&file_format=%s",
                 CORUM_API, CORUM_FILE_ID, CORUM_FORMAT)
  # Try verified first; the Helmholtz host presents an incomplete chain, so fall
  # back to relaxed peer verification for THIS host only.
  ok <- tryCatch({
    httr::GET(url, httr::write_disk(dest, overwrite = TRUE), httr::timeout(300)); TRUE
  }, error = function(e) FALSE)
  if (!ok || !file.exists(dest) || file.info(dest)$size == 0) {
    message("[corum] verified fetch failed (incomplete server TLS chain); ",
            "retrying with relaxed verification for mips.helmholtz-muenchen.de")
    httr::GET(url, httr::config(ssl_verifypeer = 0L),
              httr::write_disk(dest, overwrite = TRUE), httr::timeout(300))
  }
  record_source("05_corum_depmap", "CORUM 5.x humanComplexes",
                url, dest, note = "CORUM 5.x schema (complex_id, subunits_gene_name)")
  message("[corum] humanComplexes rows: ", length(readLines(dest)) - 1L)
  invisible(dest)
}

## ---- DepMap: download a figshare file + strip Entrez from gene columns ------
strip_entrez_header <- function(in_csv, out_csv) {
  # First column is the cell-line key (blank header); gene columns are
  # 'SYMBOL (EntrezID)'. Rewrite header only; stream the body unchanged.
  con <- file(in_csv, "r"); hdr <- readLines(con, n = 1); close(con)
  fields <- strsplit(hdr, ",")[[1]]
  fields[-1] <- sub("\\s*\\([0-9]+\\)$", "", fields[-1])  # drop ' (EntrezID)'
  # write new header + append original body (skip old header)
  tmp <- paste0(out_csv, ".tmp")
  writeLines(paste(fields, collapse = ","), tmp)
  # append remaining lines efficiently
  file.append(tmp, in_csv_body <- {
    b <- paste0(in_csv, ".body")
    system2("tail", c("-n", "+2", shQuote(in_csv)), stdout = b); b
  })
  file.rename(tmp, out_csv); unlink(in_csv_body)
  invisible(out_csv)
}

depmap_fetch <- function() {
  for (i in seq_len(nrow(DEPMAP_FILES))) {
    fig_url <- DEPMAP_FILES$figshare_url[i]
    fig_nm  <- DEPMAP_FILES$figshare_name[i]
    ana_nm  <- DEPMAP_FILES$analysis_name[i]
    raw_dl  <- file.path(DEPMAP_DIR, paste0(".dl_", fig_nm))
    out     <- file.path(DEPMAP_DIR, ana_nm)
    if (file.exists(out) && file.info(out)$size > 0) { message("[depmap] cached: ", ana_nm); next }
    fetch(fig_url, raw_dl)
    if (grepl("^Omics(CNGene|Expression)", fig_nm)) {
      strip_entrez_header(raw_dl, out); unlink(raw_dl)          # reformat gene headers
    } else {
      file.rename(raw_dl, out)                                  # Model.csv verbatim
    }
    record_source("05_corum_depmap", paste0("DepMap 23Q4 ", fig_nm, " -> ", ana_nm),
                  fig_url, out, note = "figshare article 24667905 (DepMap 23Q4 Public)")
    message("[depmap] wrote ", ana_nm)
  }
  ## CCLE proteomics (not part of the DepMap omics release).
  prot <- file.path(DEPMAP_DIR, "Proteomics.csv")
  if (!(file.exists(prot) && file.info(prot)$size > 0)) {
    ok <- tryCatch({ fetch(DEPMAP_PROTEOMICS_URL, paste0(prot, ".gz")); TRUE },
                   error = function(e) { message("[depmap] proteomics fetch failed: ",
                                                  conditionMessage(e)); FALSE })
    if (ok) {
      R.utils::gunzip(paste0(prot, ".gz"), destname = prot, overwrite = TRUE, remove = TRUE)
      record_source("05_corum_depmap", "CCLE proteomics (protein_quant_current_normalized)",
                    DEPMAP_PROTEOMICS_URL, prot,
                    note = "CCLE/Gygi proteomics; NOT in DepMap 23Q4 omics figshare")
    }
  }
  ## Arm-level CNAs: only if a canonical source URL is configured.
  arm_url <- Sys.getenv("FS_DEPMAP_ARM_URL", unset = "")
  if (nzchar(arm_url)) {
    arm <- file.path(DEPMAP_DIR, "Arm-level_CNAs.csv")
    tryCatch({ fetch(arm_url, arm)
               record_source("05_corum_depmap", "DepMap arm-level CNAs", arm_url, arm) },
             error = function(e) message("[depmap] arm-level fetch failed: ", conditionMessage(e)))
  } else {
    message("[depmap] NOTE: Arm-level_CNAs.csv source not configured (FS_DEPMAP_ARM_URL); ",
            "skipping. It is derived arm-level CN and used only by empirical validation.")
  }
}

## ---- Run --------------------------------------------------------------------
if (identical(Sys.getenv("FS_CORUM_DEPMAP_RUN"), "1") || !interactive()) {
  corum_fetch(); depmap_fetch()
} else {
  message("[corum_depmap] functions defined. Set FS_CORUM_DEPMAP_RUN=1 to execute downloads.")
}
