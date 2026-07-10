# =============================================================================
# 07_annotate.R  --  Gene annotation joins WITHOUT row inflation
# -----------------------------------------------------------------------------
# Adds, per gene: chromosome/arm/cytoband, ensembl_id (deterministic), surface
# and secreted boolean flags, and GTEx normal-tissue TPM for the row's cancer
# type.
#
# BUG A ROOT CAUSE FIX: genemap has 1,282 symbols with >1 row (multi-ensembl).
# The old pipeline joined genemap onto the omic table on `symbol`, multiplying
# those rows. Here genemap is COLLAPSED to one deterministic row per symbol
# (smallest ensembl_id, first arm/band) BEFORE joining, and the one-row-per-
# (gene,caseid) invariant is re-asserted after every join.
#
# GTEx: columns ARE cancer types (tissue mapping applied upstream). LUAD and
# LSCC both map to Lung and are therefore identical -- documented collision.
#
# Emits:
#   out/tables/omic_table_annotated.parquet
#   out/figures/annotation_coverage.png
# =============================================================================

suppressWarnings(suppressMessages({ library(data.table); library(ggplot2); library(arrow) }))
if (!exists("IN")) source(file.path("from_source", "_load_inputs.R"))

omic <- read_parquet(file.path(DIR_TAB, "omic_table_clean.parquet"), as_data_frame = TRUE)
setDT(omic)
n0 <- nrow(omic)

## ---- genemap collapsed to ONE row per symbol -------------------------------
gm <- as.data.table(IN$genemap)[, .(symbol, Chromosome, arm, cytogenetic_location, ensembl_id)]
setorder(gm, symbol, ensembl_id)                       # deterministic pick
gm1 <- gm[, .(chromosome = Chromosome[1], arm = arm[1],
              cytoband = cytogenetic_location[1], ensembl_id = ensembl_id[1]),
          by = symbol]
stopifnot(!any(duplicated(gm1$symbol)))

omic <- merge(omic, gm1, by.x = "gene", by.y = "symbol", all.x = TRUE)
stopifnot(nrow(omic) == n0)                            # no inflation

## ---- surface / secreted boolean flags --------------------------------------
## Optional TCSA/HPA surfaceome: if those auxiliary downloads failed the sets
## are absent (IN has no surface.genes/secreted.genes); flags resolve empty and
## the gap is recorded, rather than fabricating membership.
if (!is.null(IN$surface.genes)) {
  sg <- as.data.table(IN$surface.genes)
  surface_sym <- unique(na.omit(sg$HGNC.Symbol.x))
} else {
  warning("surface.genes absent; is_surface set FALSE for all genes")
  surface_sym <- character(0)
}
secreted_sym <- if (!is.null(IN$secreted.genes))
  unique(na.omit(IN$secreted.genes)) else character(0)
omic[, is_surface  := gene %in% surface_sym]
omic[, is_secreted := gene %in% secreted_sym]

## ---- GTEx normal-tissue TPM for the row's cancer type ----------------------
gt <- as.data.table(IN$gtex, keep.rownames = "gene")
gt_long <- melt(gt, id.vars = "gene", variable.name = "tumor_code",
                value.name = "gtex_tpm", variable.factor = FALSE)
omic <- merge(omic, gt_long, by = c("gene", "tumor_code"), all.x = TRUE)
stopifnot(nrow(omic) == n0)

## ---- Re-assert invariant + save --------------------------------------------
setkey(omic, gene, caseid)
stopifnot(sum(duplicated(omic[, .(gene, caseid)])) == 0)
setcolorder(omic, c("gene","caseid","tumor_code","chromosome","arm","cytoband","ensembl_id"))
write_parquet(omic, file.path(DIR_TAB, "omic_table_annotated.parquet"))

## ---- Figure: annotation coverage (fraction of genes carrying each) ---------
genes_all <- unique(omic$gene)
cov <- data.table(
  annotation = c("chromosome","ensembl_id","surface","secreted","gtex_tpm"),
  frac = c(
    mean(genes_all %in% gm1[!is.na(chromosome), symbol]),
    mean(genes_all %in% gm1[!is.na(ensembl_id), symbol]),
    mean(genes_all %in% surface_sym),
    mean(genes_all %in% secreted_sym),
    mean(genes_all %in% gt$gene)
  )
)
cov[, annotation := factor(annotation, levels = annotation[order(frac)])]
p <- ggplot(cov, aes(frac, annotation)) +
  geom_col(fill = "#4C72B0") +
  geom_text(aes(label = sprintf("%.1f%%", 100*frac)), hjust = -0.1, size = 3) +
  scale_x_continuous(limits = c(0, 1.08), expand = c(0,0)) +
  labs(title = "Annotation coverage across genes in the omic table",
       subtitle = sprintf("%d distinct genes total", length(genes_all)),
       x = "fraction of genes carrying annotation", y = NULL) +
  theme_minimal(base_size = 11) + theme(panel.grid.major.y = element_blank())
ggsave(file.path(DIR_FIG, "annotation_coverage.png"), p, width = 6.5, height = 3.6, dpi = 200)

## ---- Diagnostics -----------------------------------------------------------
message("[annotate] rows: ", nrow(omic), " (unchanged from ", n0, ": ", nrow(omic)==n0, ")")
message("[annotate] genes with chromosome: ", round(100*cov$frac[cov$annotation=="chromosome"],1), "%")
message("[annotate] surface genes flagged: ", sum(omic$is_surface[!duplicated(omic$gene)]),
        "  secreted: ", sum(omic$is_secreted[!duplicated(omic$gene)]))
message("[annotate] gtex LUAD==LSCC identical: ",
        isTRUE(all.equal(gt$LUAD, gt$LSCC)))
print(cov)
