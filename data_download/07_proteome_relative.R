# =============================================================================
# 03_proteome.R  --  Tumor-vs-normal relative protein expression + flags
# -----------------------------------------------------------------------------
# Input : IN$cptac.data  (7495 genes x 1051 aliquots, imputed TMT log-ratios)
#         cptac.tmt.RData::cptac.data (13645 x 1051, PRE-imputation, has NAs)
#         master_samples.csv (aliquot -> caseid, tumor_code, group)
#
# Method (explicit, no silent imputation):
#   * Per cancer type with >= MIN_NORMAL_REF adjacent-normal aliquots, form a
#     normal reference = per-gene MEAN over that type's normal aliquots.
#   * relative expression (log fold change) = tumor log-ratio - normal_ref.
#     TMT values are already log-ratios, so subtraction = log2-scale LFC.
#   * has_normal_ref = FALSE for cancer types with < MIN_NORMAL_REF normals
#     (GBM has zero) -> relative expression left as NA, NOT imputed.
#   * protein_measured = was this gene x tumor-aliquot actually quantified in
#     the PRE-imputation matrix (TRUE) or filled by upstream imputation (FALSE)?
#   * n_normal_ref = number of normals backing each type's reference.
#
# Emits:
#   out/tables/protein_relative.parquet  (tidy: gene, caseid, tumor_code,
#         prot_relative, prot_abs, protein_measured, has_normal_ref, n_normal_ref)
#   out/figures/protein_missingness.png
#   out/figures/protein_qc.png
# =============================================================================

suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(arrow)
}))
if (!exists("IN"))     source(file.path("from_source", "_load_inputs.R"))
master <- fread(file.path(DIR_TAB, "master_samples.csv"))

## ---- Matrices --------------------------------------------------------------
Mimp <- IN$cptac.data                                   # imputed log-ratio
tmt_env <- new.env()
load(INPUTS_AUX$proteome_tmt, envir = tmt_env)   # pre-imputation TMT (own env: avoids cptac.data collision)
Mraw <- as.matrix(tmt_env$cptac.data)                   # pre-imputation (NAs)

genes    <- rownames(Mimp)
raw_has  <- genes %in% rownames(Mraw)                   # 7471/7495 present in raw

## ---- Sample groups ---------------------------------------------------------
tum  <- master[group == "Tumor"  & is_gtex_ref == FALSE & aliquot %in% colnames(Mimp)]
norm <- master[group == "Normal" & is_gtex_ref == FALSE & aliquot %in% colnames(Mimp)]

norm_by_type <- split(norm$aliquot, norm$tumor_code)
n_normal_ref <- sapply(norm_by_type, length)
types_with_ref <- names(n_normal_ref)[n_normal_ref >= MIN_NORMAL_REF]

## ---- Normal reference (per-gene mean over normals of each type) ------------
norm_ref <- sapply(types_with_ref, function(tc) {
  al <- intersect(norm_by_type[[tc]], colnames(Mimp))
  rowMeans(Mimp[, al, drop = FALSE])
})  # genes x types_with_ref

## ---- Build tidy tumor table ------------------------------------------------
build_type <- function(tc) {
  al <- tum[tumor_code == tc, aliquot]
  case <- tum[tumor_code == tc, caseid]
  sub  <- Mimp[, al, drop = FALSE]                      # genes x tumors
  has_ref <- tc %in% types_with_ref
  rel <- if (has_ref) sub - norm_ref[, tc] else matrix(NA_real_, nrow(sub), ncol(sub))
  # measured flag from raw matrix (aligned to genes; genes absent in raw = FALSE)
  meas <- matrix(FALSE, nrow(sub), ncol(sub), dimnames = dimnames(sub))
  g_in <- genes[raw_has]
  meas[g_in, ] <- !is.na(Mraw[g_in, al, drop = FALSE])
  data.table(
    gene            = rep(genes, times = ncol(sub)),
    caseid          = rep(case,  each  = nrow(sub)),
    tumor_code      = tc,
    prot_abs        = as.vector(sub),
    prot_relative   = as.vector(rel),
    protein_measured= as.vector(meas),
    has_normal_ref  = has_ref,
    n_normal_ref    = if (has_ref) n_normal_ref[[tc]] else 0L
  )
}
prot <- rbindlist(lapply(TUMOR_CODES, build_type))
setkey(prot, gene, caseid)

# Invariant: one row per (gene, caseid)
stopifnot(!any(duplicated(prot[, .(gene, caseid)])))

write_parquet(prot, file.path(DIR_TAB, "protein_relative.parquet"))

## ---- Figure 1: missingness (fraction measured) per cancer type -------------
miss <- prot[, .(frac_measured = mean(protein_measured)), by = tumor_code]
miss[, tumor_code := factor(tumor_code, levels = TUMOR_CODES)]
p1 <- ggplot(miss, aes(tumor_code, frac_measured)) +
  geom_col(fill = "#4C72B0", width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", frac_measured)), vjust = -0.3, size = 3) +
  scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  labs(title = "Fraction of protein values actually measured (pre-imputation)",
       subtitle = "Upstream QC removed high-missingness genes; residual gaps are flagged not imputed",
       x = NULL, y = "fraction measured") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank())
ggsave(file.path(DIR_FIG, "protein_missingness.png"), p1, width = 7, height = 4.2, dpi = 200)

## ---- Figure 2: relative-expression density per type (has-ref types) --------
pq <- prot[has_normal_ref == TRUE & !is.na(prot_relative)]
p2 <- ggplot(pq, aes(prot_relative, colour = tumor_code)) +
  geom_density() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(title = "Tumor-vs-normal relative protein expression (log fold change)",
       subtitle = "GBM excluded (no adjacent-normal proteome reference)",
       x = "relative expression (tumor - normal mean, log-ratio scale)",
       y = "density", colour = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "right")
ggsave(file.path(DIR_FIG, "protein_qc.png"), p2, width = 7, height = 4.2, dpi = 200)

## ---- Diagnostics -----------------------------------------------------------
message("[proteome] tidy rows: ", nrow(prot),
        "  genes: ", length(genes), "  tumor caseids: ", length(unique(prot$caseid)))
message("[proteome] types with normal ref (>=", MIN_NORMAL_REF, "): ",
        paste(types_with_ref, collapse=","))
message("[proteome] NO normal ref: ",
        paste(setdiff(TUMOR_CODES, types_with_ref), collapse=","))
message("[proteome] overall fraction measured: ", round(mean(prot$protein_measured), 4))
print(miss)
