# =============================================================================
# 02_harmonize_samples.R  --  Harmonize sample identifiers to canonical caseid
# -----------------------------------------------------------------------------
# Canonical sample key = CPTAC caseid in DASH form (C3L-xxxxx / C3N-xxxxx;
# legacy LUAD pilot ids 11LUxxx are kept as-is). Methylation and RNA column
# names arrive in DOT form (C3L.xxxxx) and are normalized to dashes here.
#
# Emits:
#   out/tables/master_samples.csv  -- one row per proteome aliquot w/ caseid,
#                                      tumor_code, group, GTEx/legacy flags
#   out/tables/sample_overlap.csv  -- per-caseid membership across the 5 omics
#   out/figures/sample_overlap.png -- barplot of caseids per omic + 5-way core
# =============================================================================

suppressWarnings(suppressMessages({ library(data.table); library(ggplot2) }))
if (!exists("IN")) source(file.path("from_source", "_load_inputs.R"))

norm_case <- function(x) gsub("\\.", "-", x)   # dot -> dash canonicalization

## ---- Master proteome sample table ------------------------------------------
ph <- as.data.table(IN$cptac.pheno)
master <- ph[, .(aliquot = AliquotID, caseid = case_id, tumor_code, group = Group)]
master[, is_gtex_ref  := grepl("^GTEX-", caseid)]
master[, is_legacy_id := grepl("^11LU", caseid)]
fwrite(master, file.path(DIR_TAB, "master_samples.csv"))

## ---- Caseid sets per omic (tumor samples) ----------------------------------
omic_sets <- list(
  proteome = unique(master[group == "Tumor" & !is_gtex_ref, caseid]),
  cn       = norm_case(colnames(IN$cndata)),
  rna      = norm_case(colnames(IN$rnadata.norm)),
  meth     = norm_case(colnames(IN$methdata.genes)),
  snv      = norm_case(unique(IN$snvpheno$caseid))
)

all_case <- sort(unique(unlist(omic_sets)))
memb <- data.table(caseid = all_case)
for (nm in names(omic_sets)) memb[[nm]] <- memb$caseid %in% omic_sets[[nm]]
memb[, n_omic := rowSums(.SD), .SDcols = names(omic_sets)]
fwrite(memb, file.path(DIR_TAB, "sample_overlap.csv"))

## ---- Overlap summary + figure ----------------------------------------------
counts <- data.table(
  set   = c(names(omic_sets), "all_5_omics", "proteome_and_cn"),
  n     = c(sapply(omic_sets, length),
            sum(memb$n_omic == 5L),
            sum(memb$proteome & memb$cn))
)
counts[, kind := ifelse(set %in% names(omic_sets), "per-omic", "intersection")]
counts[, set := factor(set, levels = set)]

p <- ggplot(counts, aes(set, n, fill = kind)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("per-omic" = "#4C72B0", "intersection" = "#C44E52")) +
  labs(title = "CPTAC tumor caseids available per omic layer",
       subtitle = "530 caseids carry all 5 omics; 655 carry the proteome+CN core pair",
       x = NULL, y = "distinct tumor caseids", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = "top")
ggsave(file.path(DIR_FIG, "sample_overlap.png"), p, width = 7, height = 4.5, dpi = 200)

## ---- Console diagnostics ----------------------------------------------------
message("[harmonize] proteome aliquots: ", nrow(master),
        " (Tumor=", sum(master$group=="Tumor"), " Normal=", sum(master$group=="Normal"), ")")
message("[harmonize] GTEx-ref rows: ", sum(master$is_gtex_ref),
        "  legacy 11LU ids: ", sum(master$is_legacy_id))
message("[harmonize] HNSCC present in proteome? ",
        "HNSCC" %in% master$tumor_code, " (tumor_codes: ",
        paste(sort(unique(na.omit(master$tumor_code))), collapse=","), ")")
message("[harmonize] caseids all-5-omics: ", sum(memb$n_omic == 5L),
        "  proteome+cn: ", sum(memb$proteome & memb$cn))
print(counts)
