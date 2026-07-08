# =============================================================================
# 04_copynumber.R  --  Copy number + CONTINUOUS ploidy + ploidy-adjusted CN
# -----------------------------------------------------------------------------
# Input : IN$cndata  (8646 genes x 657 caseids, ASCAT gene-level absolute CN,
#                     integer, median 2)
#         IN$genemap (gene symbol -> Chromosome)
#
# FIX B (ploidy truncation): the old pipeline fed per-chromosome float CN into
# Python np.bincount, which truncates toward zero (1.8 -> 1) and forced 129
# tumours to an implausible ploidy = 1. Here ploidy is CONTINUOUS:
#     ploidy(caseid) = median over chromosomes of
#                      [ per-(caseid, chromosome) median gene CN ]
# No rounding at any step.
#
# cn_adjusted = cn / ploidy   (relative to the sample's own genome baseline).
#
# NOTE (inherited upstream limitation, documented not re-derived): the gene-level
# CN in cndata was produced by ASCAT gene_level files; gene CN at segment
# breakpoints reflects the upstream gene<-segment assignment. Raw .seg files are
# not present, so segment-level re-aggregation is out of scope for this rebuild.
#
# Emits:
#   out/tables/cn_table.parquet        (tidy: gene, caseid, cn, cn_adjusted)
#   out/tables/ploidy_table.csv        (caseid, ploidy_continuous, n_chrom)
#   out/figures/ploidy_distribution.png (old integer vs new continuous)
#   out/figures/cn_adjusted_density.png
# =============================================================================

suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(arrow)
}))
if (!exists("IN")) source(file.path("from_source", "_load_inputs.R"))

norm_case <- function(x) gsub("\\.", "-", x)

cn  <- as.matrix(IN$cndata)
colnames(cn) <- norm_case(colnames(cn))
gm  <- as.data.table(IN$genemap)[, .(symbol, Chromosome)]
gm  <- unique(gm[!is.na(Chromosome)], by = "symbol")     # one chrom per symbol
chrom_of <- gm$Chromosome[match(rownames(cn), gm$symbol)]

## ---- Continuous ploidy per caseid ------------------------------------------
# Long form for chromosome-wise medians
cn_dt <- as.data.table(cn, keep.rownames = "gene")
cn_long <- melt(cn_dt, id.vars = "gene", variable.name = "caseid",
                value.name = "cn", variable.factor = FALSE)
# Precompute chromosome vector OUTSIDE the data.table j-expression: match() with
# rownames(cn) inside `:=` resolves in the wrong frame and yields all-NA.
cn_rn     <- rownames(cn)
chrom_vec <- chrom_of[match(cn_long$gene, cn_rn)]
cn_long[, chrom := chrom_vec]

# per (caseid, chromosome) median CN -> per caseid median-over-chromosomes
chrom_med <- cn_long[!is.na(chrom) & !is.na(cn),
                     .(chrom_med = median(cn)), by = .(caseid, chrom)]
ploidy <- chrom_med[, .(ploidy_continuous = median(chrom_med),
                        n_chrom = .N), by = caseid]
fwrite(ploidy, file.path(DIR_TAB, "ploidy_table.csv"))

## ---- Ploidy-adjusted CN (tidy) ---------------------------------------------
cn_long <- merge(cn_long[, .(gene, caseid, cn)],
                 ploidy[, .(caseid, ploidy_continuous)], by = "caseid", all.x = TRUE)
cn_long[, cn_adjusted := cn / ploidy_continuous]
setkey(cn_long, gene, caseid)
stopifnot(!any(duplicated(cn_long[, .(gene, caseid)])))
write_parquet(cn_long[, .(gene, caseid, cn, cn_adjusted)],
              file.path(DIR_TAB, "cn_table.parquet"))

## ---- Figure 1: continuous ploidy distribution (vs old build if available) --
# The old (integer, np.bincount) ploidy CSV exists only in the legacy tree. In a
# from-scratch run it is absent, so the comparison overlay is optional: if a
# path is provided via FS_OLD_PLOIDY_CSV (and exists) we overlay it, otherwise
# we plot the new continuous distribution alone.
old_ploidy_csv <- Sys.getenv("FS_OLD_PLOIDY_CSV", unset = "")
have_old <- nzchar(old_ploidy_csv) && file.exists(old_ploidy_csv)
if (have_old) {
  old   <- fread(old_ploidy_csv)
  old_p <- old[!is.na(ploidy), .(caseid, ploidy_old = as.numeric(ploidy))]
  plot_dt <- rbind(
    data.table(ploidy = ploidy$ploidy_continuous, method = "new (continuous)"),
    data.table(ploidy = old_p$ploidy_old,          method = "old (integer, np.bincount)")
  )
  sub <- sprintf("old: %d samples at ploidy=1; new: %d samples <= 1",
                 sum(old_p$ploidy_old == 1), sum(ploidy$ploidy_continuous <= 1))
} else {
  plot_dt <- data.table(ploidy = ploidy$ploidy_continuous, method = "new (continuous)")
  sub <- sprintf("continuous ploidy; %d samples <= 1 (no implausible ploidy=1 pile-up)",
                 sum(ploidy$ploidy_continuous <= 1))
}
p1 <- ggplot(plot_dt, aes(ploidy, fill = method)) +
  geom_histogram(binwidth = 0.1, position = "identity", alpha = 0.55) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#C44E52") +
  scale_fill_manual(values = c("new (continuous)" = "#4C72B0",
                               "old (integer, np.bincount)" = "#DD8452")) +
  labs(title = "Continuous ploidy removes the implausible ploidy = 1 pile-up",
       subtitle = sub, x = "tumour ploidy", y = "caseids", fill = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(DIR_FIG, "ploidy_distribution.png"), p1, width = 7, height = 4.4, dpi = 200)

## ---- Figure 2: cn_adjusted density -----------------------------------------
set.seed(1)
samp <- cn_long[!is.na(cn_adjusted)][sample(.N, min(2e5, .N))]
p2 <- ggplot(samp, aes(cn_adjusted)) +
  geom_density(fill = "#4C72B0", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  coord_cartesian(xlim = c(0, 4)) +
  labs(title = "Ploidy-adjusted copy number (cn / continuous ploidy)",
       subtitle = "Mode near 1 = balanced relative to each tumour's own genome baseline",
       x = "cn_adjusted", y = "density") +
  theme_minimal(base_size = 11)
ggsave(file.path(DIR_FIG, "cn_adjusted_density.png"), p2, width = 7, height = 4.2, dpi = 200)

## ---- Diagnostics -----------------------------------------------------------
message("[cn] caseids: ", nrow(ploidy),
        "  ploidy median: ", round(median(ploidy$ploidy_continuous), 3),
        "  range: ", paste(round(range(ploidy$ploidy_continuous), 2), collapse=".."))
message("[cn] samples <=1 (new): ", sum(ploidy$ploidy_continuous <= 1))
if (have_old) {
  cmp <- merge(ploidy[, .(caseid, ploidy_continuous)], old_p, by = "caseid")
  message("[cn] old integer ploidy=1: ", sum(old_p$ploidy_old == 1))
  message("[cn] old-vs-new ploidy correlation (matched ", nrow(cmp), "): ",
          round(cor(cmp$ploidy_continuous, cmp$ploidy_old), 3))
}
message("[cn] tidy CN rows: ", nrow(cn_long))
