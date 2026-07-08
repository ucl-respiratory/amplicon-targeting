# =============================================================================
# 02_proteome_pdc.R  --  CPTAC Discovery proteome from the PDC GraphQL API
# -----------------------------------------------------------------------------
# Replaces the original download.data.R proteome block (~lines 51-343), which
# required (a) two hand-exported PDC manifest CSVs, (b) per-study TMT .tsv files
# parsed for the 'Unshared.Log.Ratio' columns, and (c) SEVEN differently
# formatted clinical/specimen .xlsx files hand-curated per cancer type to map
# aliquot -> case and tumour/normal group.
#
# ALL of that is now retrieved programmatically from the PDC GraphQL API:
#   * quantDataMatrix(data_type:"unshared_log2_ratio")  -> gene x aliquot matrix
#       (the CDAP gene-level unshared log2 ratios; the exact quantity the
#        original parsed from the *.Unshared.Log.Ratio columns).
#   * biospecimenPerStudy()                             -> aliquot -> case,
#       sample_type (Primary Tumor / Solid Tissue Normal) per study.
#
# The downstream pre-processing (normalize -> QC -> impute -> matched T-N) is
# reproduced EXACTLY as the original (limma::normalizeBetweenArrays; drop genes
# NA>10% and samples NA>40%; per-gene mean impute via Hmisc::impute; matched
# tumour-minus-normal per case).
#
# Emits (object names/shapes matching the cached pipeline contract):
#   PARSED/cptac.tmt.RData          cptac.data(df genes x aliquots, PRE-impute),
#                                   cptac.pheno(df aliquots x pheno)
#   PARSED/cptac.norm.qc.imp.RData  cptac.data(matrix, normalized+QC+imputed),
#                                   cptac.pheno(df, QC-filtered)
#   PARSED/cptac.matched.RData      cptac.matched.data(df genes x cases, T-N),
#                                   cptac.matched.pheno(df)
#
# NOTE on pheno columns: the original cptac.pheno carried ~14 columns including
# clinical fields (age, BMI, tumour_site, histology, smoking/alcohol history)
# parsed from the per-study clinical xlsx. Those clinical fields are NOT used by
# any assembly stage (only case_id, tumor_code, AliquotID, Group, Code are read
# downstream) and do not appear in any final analysis table, so this rebuild
# emits the 5 identity/group columns the pipeline actually consumes (plus a
# placeholder type_of_analyzed_samples). If full clinical annotation is later
# wanted, PDC exposes it via clinicalPerStudy(pdc_study_id) / clinicalMetadata.
# =============================================================================

suppressWarnings(suppressMessages({
  library(data.table); library(jsonlite); library(httr)
  library(limma); library(Hmisc)
}))
if (!exists("PARSED")) source(file.path("from_source", "00_config.R"))

## ---- PDC GraphQL helper -----------------------------------------------------
pdc_gql <- function(query, tries = 3, pause = 3) {
  for (k in seq_len(tries)) {
    r <- tryCatch(
      httr::POST(PDC_GRAPHQL,
                 httr::content_type_json(),
                 body = jsonlite::toJSON(list(query = query), auto_unbox = TRUE),
                 httr::timeout(300)),
      error = function(e) NULL)
    if (!is.null(r) && httr::status_code(r) == 200) {
      j <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"),
                              simplifyVector = FALSE)
      if (is.null(j$errors)) return(j$data)
      msg <- paste(vapply(j$errors, function(e) e$message %||% "", ""), collapse = "; ")
      if (!grepl("not found", msg, ignore.case = TRUE) && k < tries) { Sys.sleep(pause); next }
      stop("PDC GraphQL error: ", msg)
    }
    Sys.sleep(pause)
  }
  stop("PDC GraphQL request failed after ", tries, " tries")
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- Per-study biospecimen (aliquot -> case, sample_type) -------------------
fetch_biospecimen <- function(pdc_study_id) {
  q <- sprintf('{ biospecimenPerStudy(pdc_study_id:"%s" acceptDUA:true){
                 aliquot_submitter_id case_submitter_id sample_type disease_type } }',
               pdc_study_id)
  bio <- pdc_gql(q)$biospecimenPerStudy
  rbindlist(lapply(bio, function(b) data.table(
    aliquot   = b$aliquot_submitter_id %||% NA_character_,
    case_id   = b$case_submitter_id    %||% NA_character_,
    sample_type = b$sample_type        %||% NA_character_,
    disease   = b$disease_type         %||% NA_character_)))
}

## ---- Per-study gene x aliquot unshared log2-ratio matrix --------------------
fetch_quant_matrix <- function(pdc_study_id) {
  q <- sprintf('{ quantDataMatrix(pdc_study_id:"%s" data_type:"unshared_log2_ratio" acceptDUA:true) }',
               pdc_study_id)
  m <- pdc_gql(q)$quantDataMatrix           # list of rows; row[[1]] is header
  hdr <- unlist(m[[1]])                     # c("Gene/Aliquot", "uuid:CPTxxxx", ...)
  aliquots <- sub("^[^:]*:", "", hdr[-1])   # strip "study_uuid:" -> aliquot_submitter_id
  body <- m[-1]
  genes <- vapply(body, function(r) r[[1]], character(1))
  vals  <- t(vapply(body, function(r) as.numeric(unlist(r)[-1]), numeric(length(aliquots))))
  df <- as.data.frame(vals, stringsAsFactors = FALSE)
  colnames(df) <- aliquots
  rownames(df) <- genes
  df
}

## ---- Assemble pooled matrix + pheno across the 7 Discovery studies ----------
message("[proteome] querying PDC for ", nrow(PDC_STUDIES), " Discovery Proteome studies")
mat_list  <- list(); bio_list <- list()
for (i in seq_len(nrow(PDC_STUDIES))) {
  sid <- PDC_STUDIES$pdc_study_id[i]; tc <- PDC_STUDIES$tumor_code[i]
  message("  [", tc, " ", sid, "] quantDataMatrix + biospecimen ...")
  mat_list[[sid]] <- fetch_quant_matrix(sid)
  b <- fetch_biospecimen(sid); b[, tumor_code := tc]; bio_list[[sid]] <- b
  record_source("02_proteome_pdc", paste0("PDC ", tc, " ", sid),
                paste0(PDC_GRAPHQL, " quantDataMatrix(unshared_log2_ratio)+biospecimenPerStudy"),
                INPUTS_AUX$proteome_tmt,
                note = sprintf("%d genes x %d aliquots; %d biospecimens",
                               nrow(mat_list[[sid]]), ncol(mat_list[[sid]]), nrow(b)))
}

## Full-outer-join the per-study matrices on gene (matches original merge all=T).
merge_on_gene <- function(a, b) {
  m <- merge(a, b, by = "row.names", all = TRUE)
  rownames(m) <- m$Row.names; m$Row.names <- NULL; m
}
cptac.data <- Reduce(merge_on_gene, mat_list)

## Pheno: one row per aliquot present in the matrix.
bio <- rbindlist(bio_list, fill = TRUE)
bio <- bio[!duplicated(aliquot)]
bio[, Group := fifelse(grepl("Tumor", sample_type, ignore.case = TRUE), "Tumor",
                fifelse(grepl("Normal", sample_type, ignore.case = TRUE), "Normal", NA_character_))]
cptac.pheno <- data.frame(
  case_id                  = bio$case_id[match(colnames(cptac.data), bio$aliquot)],
  tumor_code               = bio$tumor_code[match(colnames(cptac.data), bio$aliquot)],
  type_of_analyzed_samples = NA,
  AliquotID                = colnames(cptac.data),
  Group                    = bio$Group[match(colnames(cptac.data), bio$aliquot)],
  Code                     = bio$tumor_code[match(colnames(cptac.data), bio$aliquot)],
  stringsAsFactors = FALSE
)
# Keep only aliquots we could map to a case (drops PDC 'Not Reported' pool refs).
keep <- !is.na(cptac.pheno$case_id) & !is.na(cptac.pheno$Group)
cptac.data  <- cptac.data[, keep, drop = FALSE]
cptac.pheno <- cptac.pheno[keep, , drop = FALSE]
save(cptac.data, cptac.pheno, file = INPUTS_AUX$proteome_tmt)   # pre-impute TMT object
message("[proteome] pooled TMT matrix: ", nrow(cptac.data), " genes x ",
        ncol(cptac.data), " aliquots (Tumor=", sum(cptac.pheno$Group=="Tumor"),
        " Normal=", sum(cptac.pheno$Group=="Normal"), ")")

## ---- Normalize -> QC -> impute (exact reproduction) -------------------------
cptac.norm <- limma::normalizeBetweenArrays(as.matrix(cptac.data))
samp.qc  <- apply(cptac.norm, 2, function(x) sum(is.na(x))) / nrow(cptac.norm)
genes.qc <- apply(cptac.norm, 1, function(x) sum(is.na(x))) / ncol(cptac.norm)
cptac.qc <- cptac.norm[genes.qc < 0.1, samp.qc < 0.4]
cptac.pheno.qc <- cptac.pheno[match(colnames(cptac.qc), cptac.pheno$AliquotID), ]

cptac.qc.imp <- t(apply(cptac.qc, 1, function(x) as.numeric(Hmisc::impute(x))))  # per-gene mean impute
colnames(cptac.qc.imp) <- colnames(cptac.qc)
cptac.data  <- cptac.qc.imp
cptac.pheno <- cptac.pheno.qc
cptac.pheno$Group[grep("Tumor",  cptac.pheno$Group, ignore.case = TRUE)] <- "Tumor"
cptac.pheno$Group[grep("Normal", cptac.pheno$Group, ignore.case = TRUE)] <- "Normal"
save(cptac.data, cptac.pheno, file = INPUTS$proteome)       # normalized+QC+imputed
message("[proteome] post-QC+impute: ", nrow(cptac.data), " genes x ", ncol(cptac.data), " aliquots")

## ---- Matched tumour-vs-normal (per case) ------------------------------------
dup_cases <- cptac.pheno$case_id[duplicated(cptac.pheno$case_id)]
mp <- rbindlist(lapply(dup_cases, function(x) {
  tid <- cptac.pheno$AliquotID[cptac.pheno$case_id == x & cptac.pheno$Group == "Tumor"]
  nid <- cptac.pheno$AliquotID[cptac.pheno$case_id == x & cptac.pheno$Group == "Normal"]
  if (length(tid) < 1 || length(nid) < 1 || is.na(tid[1]) || is.na(nid[1])) return(NULL)
  data.table(case_id = x, tumor_id = tid[1], normal_id = nid[1])
}))
cptac.matched.pheno <- merge(as.data.frame(mp),
                             cptac.pheno[!duplicated(cptac.pheno$case_id), ],
                             by = "case_id", all.x = TRUE, all.y = FALSE)
cptac.matched.data <- do.call(cbind, lapply(seq_len(nrow(cptac.matched.pheno)), function(i) {
  y <- data.frame(cptac.data[, cptac.matched.pheno$tumor_id[i]] -
                  cptac.data[, cptac.matched.pheno$normal_id[i]])
  colnames(y) <- cptac.matched.pheno$case_id[i]; y
}))
save(cptac.matched.data, cptac.matched.pheno, file = INPUTS_AUX$matched)
message("[proteome] matched T-N pairs: ", nrow(cptac.matched.pheno), " cases")
