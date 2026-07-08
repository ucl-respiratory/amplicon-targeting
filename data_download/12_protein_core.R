# =============================================================================
# 08_protein_core.R  --  Proteome-core view of the annotated omic table
# -----------------------------------------------------------------------------
# The full omic_table_annotated.parquet is the UNION of all gene universes
# (29,313 genes) so that RNA-only / methylation-only genes are not silently
# dropped. For most downstream modelling the relevant scope is genes that HAVE
# protein data (the CN->protein prediction target). Rather than load the 438 MB
# union and filter to has_protein on every run, we materialize that filtered
# view once here.
#
# This view is the modern equivalent of the old str_omic_table.csv scope
# (proteome gene grid), but carries every fix (continuous ploidy, collapsed
# genemap, methylation replicate averaging, tumor_code backfill) and all
# annotations. Same columns, same (gene, caseid) grain, same invariants.
#
# Emits: out/tables/omic_table_protein_core.parquet
# =============================================================================

suppressWarnings(suppressMessages({ library(data.table); library(arrow) }))
if (!exists("DIR_TAB")) source(file.path("from_source", "00_config.R"))

omic <- read_parquet(file.path(DIR_TAB, "omic_table_annotated.parquet"), as_data_frame = TRUE)
setDT(omic)
n_full <- nrow(omic)

core <- omic[has_protein == TRUE]

## ---- Invariants: still one row per (gene, caseid); protein complete --------
stopifnot(sum(duplicated(core[, .(gene, caseid)])) == 0)
stopifnot(all(core$has_protein))
stopifnot(!any(is.na(core$prot_abs)))          # protein present by definition

setkey(core, gene, caseid)
write_parquet(core, file.path(DIR_TAB, "omic_table_protein_core.parquet"))

## ---- Diagnostics -----------------------------------------------------------
message("[core] full union: ", n_full, " rows -> protein core: ", nrow(core),
        " rows (", round(100*nrow(core)/n_full, 1), "%)")
message("[core] genes: ", uniqueN(core$gene), "  caseids: ", uniqueN(core$caseid))
nr <- sapply(core[, .(prot_abs, prot_relative, cn, rna_log2, meth_beta, gtex_tpm)],
             function(x) round(mean(is.na(x)), 3))
print(nr)