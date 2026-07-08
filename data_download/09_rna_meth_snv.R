# =============================================================================
# 05_rna_meth_snv.R  --  RNA, methylation, SNV to tidy (gene, caseid) long
# -----------------------------------------------------------------------------
# RNA  : IN$rnadata.norm (23961 x 661, DESeq2-normalized counts, NOT logged).
#        Emit rna_norm (as-is) + rna_log2 = log2(rna_norm + 1) for modelling.
# METH : IN$methdata.genes (25186 x 626, gene-level beta values 0-1).
# SNV  : IN$snvdata (long MAF, already filtered to coding-impact classes).
#        Collapse to per (gene, caseid): n_variants, is_mutated, is_nonsilent.
#
# All three are canonicalized to dash-form caseids and written as tidy parquet
# with a strict one-row-per-(gene, caseid) invariant so downstream joins are 1:1.
#
# Emits:
#   out/tables/rna_table.parquet   (gene, caseid, rna_norm, rna_log2)
#   out/tables/meth_table.parquet  (gene, caseid, meth_beta)
#   out/tables/snv_table.parquet   (gene, caseid, n_variants, is_mutated, is_nonsilent)
#   out/figures/rna_vs_cn_scatter.png
#   out/figures/mutation_frequency_top30.png
# =============================================================================

suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(arrow)
}))
if (!exists("IN")) source(file.path("from_source", "_load_inputs.R"))
norm_case <- function(x) gsub("\\.", "-", x)

melt_matrix <- function(mat, value_name) {
  m <- as.matrix(mat); colnames(m) <- norm_case(colnames(m))
  dt <- as.data.table(m, keep.rownames = "gene")
  long <- melt(dt, id.vars = "gene", variable.name = "caseid",
               value.name = value_name, variable.factor = FALSE)
  long
}

## ---- RNA --------------------------------------------------------------------
rna <- melt_matrix(IN$rnadata.norm, "rna_norm")
rna[, rna_log2 := log2(rna_norm + 1)]
setkey(rna, gene, caseid); stopifnot(!any(duplicated(rna[, .(gene, caseid)])))
write_parquet(rna, file.path(DIR_TAB, "rna_table.parquet"))

## ---- Methylation ------------------------------------------------------------
# The raw beta matrix has DUPLICATE column names: cases with multiple
# methylation aliquots (e.g. C3L-00103 measured 3x) arrive as identically
# named columns, which R's make.names disambiguated to C3L.00103, C3L.00103.1,
# C3L.00103.2. A naive dot->dash normalization turns these into PHANTOM
# caseids (C3L-00103-1) that never match the canonical case. Fix: map every
# column to its canonical caseid (strip the make.names replicate index) and
# average beta across replicate aliquots per (gene, caseid).
meth_canon <- function(x) {
  x <- sub("^(C3[LN]\\.[0-9]+)\\.[0-9]+$", "\\1", x)   # drop make.names .N index
  gsub("\\.", "-", x)                                    # dots -> dashes
}
mm <- as.matrix(IN$methdata.genes)
colnames(mm) <- meth_canon(colnames(mm))
n_rep_cols <- sum(duplicated(colnames(mm)) | duplicated(colnames(mm), fromLast = TRUE))
meth <- as.data.table(mm, keep.rownames = "gene")
meth <- melt(meth, id.vars = "gene", variable.name = "caseid",
             value.name = "meth_beta", variable.factor = FALSE)
# average across replicate aliquots (single-aliquot cases unaffected: mean of 1)
meth <- meth[, .(meth_beta = mean(meth_beta, na.rm = TRUE)), by = .(gene, caseid)]
meth[is.nan(meth_beta), meth_beta := NA_real_]
setkey(meth, gene, caseid); stopifnot(!any(duplicated(meth[, .(gene, caseid)])))
write_parquet(meth, file.path(DIR_TAB, "meth_table.parquet"))
message("[meth] collapsed ", n_rep_cols, " replicate columns -> ",
        uniqueN(meth$caseid), " canonical caseids")

## ---- SNV --------------------------------------------------------------------
sv <- as.data.table(IN$snvdata)
sv[, caseid := norm_case(caseid)]
silent_classes <- c("Silent", "Synonymous_Variant")   # none present, robust anyway
snv <- sv[, .(n_variants  = .N,
              is_nonsilent = any(!classification %in% silent_classes)),
          by = .(gene, caseid)]
snv[, is_mutated := TRUE]                              # any row = mutated
setkey(snv, gene, caseid); stopifnot(!any(duplicated(snv[, .(gene, caseid)])))
write_parquet(snv, file.path(DIR_TAB, "snv_table.parquet"))

## ---- Figure 1: CN->RNA dosage sanity check ---------------------------------
# The dosage effect is WITHIN-gene, ACROSS-sample (not across genes: cross-gene
# means are dominated by intrinsic per-gene expression). For each gene we take
# the Spearman correlation of CN vs RNA over caseids; the distribution should be
# centred positive.
cn <- read_parquet(file.path(DIR_TAB, "cn_table.parquet"), as_data_frame = TRUE)
setDT(cn)
mrg <- merge(rna[, .(gene, caseid, rna_log2)], cn[, .(gene, caseid, cn)],
             by = c("gene", "caseid"))
per_gene <- mrg[, {
  if (.N >= 20 && length(unique(cn)) > 1 && length(unique(rna_log2)) > 1)
    .(rho = cor(cn, rna_log2, method = "spearman")) else .(rho = NA_real_)
}, by = gene][!is.na(rho)]
med_rho  <- median(per_gene$rho)
frac_pos <- mean(per_gene$rho > 0)
p1 <- ggplot(per_gene, aes(rho)) +
  geom_histogram(binwidth = 0.025, fill = "#4C72B0") +
  geom_vline(xintercept = 0,       linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = med_rho, colour = "#C44E52", linewidth = 0.8) +
  labs(title = "Per-gene CN->RNA dosage correlation (across caseids)",
       subtitle = sprintf("median rho = %.2f, %.0f%% of genes positive \u2014 expected dosage effect",
                          med_rho, 100 * frac_pos),
       x = "per-gene Spearman rho (copy number vs RNA)", y = "genes") +
  theme_minimal(base_size = 11)
ggsave(file.path(DIR_FIG, "rna_vs_cn_scatter.png"), p1, width = 6.5, height = 4.4, dpi = 200)

## ---- Figure 2: top-30 most frequently mutated genes ------------------------
mut_freq <- snv[, .(n_cases = uniqueN(caseid)), by = gene][order(-n_cases)][1:30]
mut_freq[, gene := factor(gene, levels = rev(gene))]
p2 <- ggplot(mut_freq, aes(n_cases, gene)) +
  geom_col(fill = "#4C72B0") +
  labs(title = "Top 30 most frequently mutated genes (CPTAC)",
       subtitle = paste0("Across ", uniqueN(snv$caseid), " caseids with SNV data"),
       x = "number of mutated caseids", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank())
ggsave(file.path(DIR_FIG, "mutation_frequency_top30.png"), p2, width = 6.5, height = 6, dpi = 200)

## ---- Diagnostics -----------------------------------------------------------
message("[rna]  rows: ", nrow(rna),  "  genes: ", uniqueN(rna$gene),  "  caseids: ", uniqueN(rna$caseid))
message("[meth] rows: ", nrow(meth), "  genes: ", uniqueN(meth$gene), "  caseids: ", uniqueN(meth$caseid))
message("[snv]  rows: ", nrow(snv),  "  genes: ", uniqueN(snv$gene),  "  caseids: ", uniqueN(snv$caseid),
        "  (long MAF rows: ", nrow(sv), ")")
message("[rna_vs_cn] per-gene median rho = ", round(med_rho, 3),
        "  frac positive = ", round(frac_pos, 3), "  (genes = ", nrow(per_gene), ")")
message("[snv] top gene: ", as.character(mut_freq$gene[nrow(mut_freq)]),
        " mutated in ", max(mut_freq$n_cases), " caseids")
