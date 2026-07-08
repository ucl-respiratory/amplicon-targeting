# =============================================================================
# 06_assemble.R  --  Assemble one tidy omic table on (gene, caseid)
# -----------------------------------------------------------------------------
# Full outer join of the five per-omic tidy tables. The join key is
# (gene, caseid). INVARIANT (fixes BUG A row inflation at the source): the
# result has EXACTLY one row per (gene, caseid) -- asserted, not deduped after.
#
# The universe of rows is the union of (gene, caseid) pairs observed in the
# proteome, CN, RNA and methylation layers (continuous omics). SNV is a sparse
# overlay: absence of an SNV row = not mutated, so SNV missing flags become
# is_mutated = FALSE / n_variants = 0 rather than NA.
#
# Emits:
#   out/tables/omic_table_clean.parquet
#   out/reports/assembly_report.txt
# =============================================================================

suppressWarnings(suppressMessages({ library(data.table); library(arrow) }))
if (!exists("DIR_TAB")) source(file.path("from_source", "00_config.R"))

rd <- function(f) { d <- read_parquet(file.path(DIR_TAB, f), as_data_frame = TRUE); setDT(d); d }
prot <- rd("protein_relative.parquet")
cn   <- rd("cn_table.parquet")
rna  <- rd("rna_table.parquet")
meth <- rd("meth_table.parquet")
snv  <- rd("snv_table.parquet")

report <- c("ASSEMBLY REPORT", "===============", "",
            sprintf("proteome : %9d rows  %6d genes  %4d caseids", nrow(prot), uniqueN(prot$gene), uniqueN(prot$caseid)),
            sprintf("cn       : %9d rows  %6d genes  %4d caseids", nrow(cn),   uniqueN(cn$gene),   uniqueN(cn$caseid)),
            sprintf("rna      : %9d rows  %6d genes  %4d caseids", nrow(rna),  uniqueN(rna$gene),  uniqueN(rna$caseid)),
            sprintf("meth     : %9d rows  %6d genes  %4d caseids", nrow(meth), uniqueN(meth$gene), uniqueN(meth$caseid)),
            sprintf("snv      : %9d rows  %6d genes  %4d caseids", nrow(snv),  uniqueN(snv$gene),  uniqueN(snv$caseid)),
            "")

## ---- Per-omic (gene, caseid) key uniqueness (guard before joining) ---------
for (nm in c("prot","cn","rna","meth","snv")) {
  d <- get(nm); stopifnot(!any(duplicated(d[, .(gene, caseid)])))
}

## ---- Universe = union of continuous-omic keys ------------------------------
keyset <- unique(rbindlist(list(
  prot[, .(gene, caseid)], cn[, .(gene, caseid)],
  rna[,  .(gene, caseid)], meth[, .(gene, caseid)]
)))
setkey(keyset, gene, caseid)
report <- c(report, sprintf("union key universe (prot|cn|rna|meth): %d rows", nrow(keyset)), "")

## ---- Sequential left joins onto the key universe ---------------------------
omic <- keyset
omic <- merge(omic, prot[, .(gene, caseid, tumor_code, prot_abs, prot_relative,
                             protein_measured, has_normal_ref, n_normal_ref)],
              by = c("gene","caseid"), all.x = TRUE)
omic <- merge(omic, cn[,   .(gene, caseid, cn, cn_adjusted)],           by = c("gene","caseid"), all.x = TRUE)
omic <- merge(omic, rna[,  .(gene, caseid, rna_norm, rna_log2)],        by = c("gene","caseid"), all.x = TRUE)
omic <- merge(omic, meth[, .(gene, caseid, meth_beta)],                 by = c("gene","caseid"), all.x = TRUE)
omic <- merge(omic, snv[,  .(gene, caseid, n_variants, is_nonsilent)],  by = c("gene","caseid"), all.x = TRUE)

## ---- SNV sparse-overlay semantics: missing = not mutated -------------------
omic[, is_mutated  := !is.na(n_variants)]
omic[is.na(n_variants),   n_variants := 0L]
omic[is.na(is_nonsilent), is_nonsilent := FALSE]

## ---- Per-omic measured flags -----------------------------------------------
omic[, has_protein := !is.na(prot_abs)]
omic[, has_cn      := !is.na(cn)]
omic[, has_rna     := !is.na(rna_norm)]
omic[, has_meth    := !is.na(meth_beta)]

## ---- Backfill tumor_code from master where proteome absent -----------------
# Map caseid -> tumor_code from ALL non-GTEx aliquots (not just tumors): a case
# whose only proteome aliquot is an adjacent normal still has a well-defined
# cancer type, and its CN/RNA/meth rows must inherit it. Filtering to Tumor
# aliquots left such normal-only cases (e.g. C3N-02379, C3N-02587) with NA
# tumor_code. (caseid, tumor_code) is consistent across a case's aliquots.
master <- fread(file.path(DIR_TAB, "master_samples.csv"))
tc_map <- unique(master[is_gtex_ref == FALSE, .(caseid, tumor_code)])
tc_map <- tc_map[!duplicated(caseid)]
omic <- merge(omic, tc_map[, .(caseid, tc_fill = tumor_code)], by = "caseid", all.x = TRUE)
omic[is.na(tumor_code), tumor_code := tc_fill][, tc_fill := NULL]

## ---- INVARIANT: one row per (gene, caseid) ---------------------------------
setkey(omic, gene, caseid)
dup_n <- sum(duplicated(omic[, .(gene, caseid)]))
stopifnot(dup_n == 0)

## ---- Column order ----------------------------------------------------------
setcolorder(omic, c("gene","caseid","tumor_code",
                    "prot_abs","prot_relative","protein_measured","has_normal_ref","n_normal_ref",
                    "cn","cn_adjusted","rna_norm","rna_log2","meth_beta",
                    "n_variants","is_mutated","is_nonsilent",
                    "has_protein","has_cn","has_rna","has_meth"))
write_parquet(omic, file.path(DIR_TAB, "omic_table_clean.parquet"))

## ---- Report ----------------------------------------------------------------
report <- c(report, "FINAL TABLE", "-----------",
  sprintf("rows: %d", nrow(omic)),
  sprintf("genes: %d   caseids: %d", uniqueN(omic$gene), uniqueN(omic$caseid)),
  sprintf("duplicate (gene,caseid) keys: %d", dup_n), "",
  "per-omic coverage (rows with data / total rows):",
  sprintf("  protein : %9d  (%.1f%%)", sum(omic$has_protein), 100*mean(omic$has_protein)),
  sprintf("  cn      : %9d  (%.1f%%)", sum(omic$has_cn),      100*mean(omic$has_cn)),
  sprintf("  rna     : %9d  (%.1f%%)", sum(omic$has_rna),     100*mean(omic$has_rna)),
  sprintf("  meth    : %9d  (%.1f%%)", sum(omic$has_meth),    100*mean(omic$has_meth)),
  sprintf("  mutated : %9d  (%.1f%%)", sum(omic$is_mutated),  100*mean(omic$is_mutated)),
  "",
  "rows with ALL FOUR continuous omics (prot+cn+rna+meth):",
  sprintf("  %d (%.1f%%)", sum(omic$has_protein & omic$has_cn & omic$has_rna & omic$has_meth),
          100*mean(omic$has_protein & omic$has_cn & omic$has_rna & omic$has_meth)))
writeLines(report, file.path(DIR_REP, "assembly_report.txt"))
cat(paste(report, collapse = "\n"), "\n")
