# =============================================================================
# 09_james_adapter.R  --  Emit a str_omic_table.csv-compatible file for James'
#                         ML pipeline (xg_boost.py / cnn.py), plus the
#                         adc_distinct_genes.csv gene filter list.
# -----------------------------------------------------------------------------
# WHY: James' Python models read a wide flat CSV with SHORT column names
#      (prot, cn, rna, meth, gtex, ...), the ECDF target prot.rel.tissue, and
#      18 DESCRIBEPROT protein-structure columns joined by UniProt accession.
#      Our rebuilt ingest uses long, explicit names and does NOT carry the ECDF
#      relative-protein columns or the structure features. This stage bridges
#      the two so James' code runs UNMODIFIED (only his two hard-coded
#      /Users/jameslam/data/ paths need repointing).
#
# INPUTS
#   out/tables/omic_table_protein_core.parquet   (our proteome-core table)
#   IN$cptac.data (Mimp) + master_samples.csv     (for the ECDF prot.rel.* cols)
#   old/data/raw/9606_value.csv (skip=21)          (DESCRIBEPROT, key=UniProt ACC)
#   old/data/ADC_atlas_..._tableS3.csv             (ADC target gene list source)
#   org.Hs.eg.db (symbol -> Entrez -> UniProt map)
#
# METHOD NOTES
#   * prot.rel.tissue / prot.rel.match / prot.rel.all reproduce James'
#     cache_gene_lookup_table.R definitions EXACTLY, computed on the imputed
#     proteome matrix Mimp (genes x aliquots) + master pheno:
#       - prot.rel.tissue = fraction of TISSUE-MATCHED normal aliquots the
#         tumour's protein value exceeds  (0..1 ECDF rank).  <-- TARGET
#       - prot.rel.all    = same but against ALL normals (any tumour type).
#       - prot.rel.match  = tumour minus the (single/mean) matched-normal value
#         for the SAME case, NA if that case has no normal aliquot.
#   * DESCRIBEPROT join: James deduped the symbol->uniprot map on gene_id, which
#     kept a non-canonical isoform for ~82% of genes and gave features to only
#     1,315 genes. We instead match on ANY of a gene's UniProt accessions and
#     pick the lexicographically-smallest ACC present in DESCRIBEPROT -> 7,280
#     genes covered. Set JAMES_STRICT_UNIPROT=1 to reproduce his 1,315-gene
#     behaviour instead.
#   * cn.capped = min(cn, 5) (James' name); we also already have cn_adjusted.
#   * cytogenetic_location = our cytoband.
#
# EMITS
#   out/tables/str_omic_table_rebuilt.csv     (James-compatible, all his cols)
#   out/tables/adc_distinct_genes.csv         (one column, no header)
#   out/tables/describeprot_gene_features.parquet  (gene-level structure table)
#   out/reports/james_adapter_coverage.txt    (column + coverage audit)
# =============================================================================

suppressWarnings(suppressMessages({
  library(data.table); library(arrow)
}))
if (!exists("DIR_TAB")) source(file.path("from_source", "00_config.R"))
.ulib <- if (nzchar(Sys.getenv("CNT_WORKSPACE"))) file.path(CNT_WORKSPACE, "rlibs") else "./rlibs"
.libPaths(c(.ulib, "./rlibs", .libPaths()))
suppressWarnings(suppressMessages(library(org.Hs.eg.db)))
STRICT <- nzchar(Sys.getenv("JAMES_STRICT_UNIPROT"))

# DescribePROT + ADC Atlas: downloaded by 04_auxiliary.R into RAW/. Fall back to
# the legacy locations if a from-scratch download was skipped.
RAW_DP  <- { p <- file.path(RAW, "describeprot_9606_value.csv")
             if (file.exists(p)) p else file.path(CNT_ROOT, "old", "data", "raw", "9606_value.csv") }
ADC_CSV <- { p <- file.path(RAW, "adc_atlas_tableS3.xlsx")
             if (file.exists(p)) p else file.path(CNT_ROOT, "old", "data",
                     "ADC_atlas_41417_2023_701_MOESM2_ESM_tableS3.csv") }

## ===========================================================================
## PART 1 -- ECDF relative-protein columns (prot.rel.tissue/all/match)
## ===========================================================================
if (!exists("IN")) source(file.path("from_source", "_load_inputs.R"))
master <- fread(file.path(DIR_TAB, "master_samples.csv"))
Mimp   <- as.matrix(IN$cptac.data)                     # genes x aliquots (imputed)

# aliquot -> group/type/case (exclude GTEx pseudo-aliquots)
ph <- master[is_gtex_ref == FALSE & aliquot %in% colnames(Mimp)]
tum_ph  <- ph[group == "Tumor"]
norm_ph <- ph[group == "Normal"]
norm_all_al <- norm_ph$aliquot                                     # all normals
norm_al_by_type <- split(norm_ph$aliquot, norm_ph$tumor_code)      # per-type

genes_M <- rownames(Mimp)

# Pre-extract normal sub-matrices once (genes x normals)
N_all <- Mimp[, norm_all_al, drop = FALSE]

# For each tumour aliquot compute the two ECDF ranks across genes at once.
# rank_all[g]    = mean( tumour[g] > N_all[g, ] )
# rank_tissue[g] = mean( tumour[g] > N_tissue[g, ] )   (tissue = same tumor_code)
ecdf_rank <- function(tvec, Nmat) {
  # tvec: length G ; Nmat: G x K ; returns length G fraction of columns exceeded
  if (ncol(Nmat) == 0L) return(rep(NA_real_, length(tvec)))
  rowMeans(Nmat < tvec)          # tvec recycled down columns (column-major): OK, G rows
}

rel_rows <- vector("list", nrow(tum_ph))
for (i in seq_len(nrow(tum_ph))) {
  al   <- tum_ph$aliquot[i]; cs <- tum_ph$caseid[i]; tc <- tum_ph$tumor_code[i]
  tvec <- Mimp[, al]
  N_tis <- { a <- norm_al_by_type[[tc]]; if (is.null(a)) Mimp[, integer(0), drop=FALSE] else Mimp[, a, drop=FALSE] }
  # matched normal: same case, Normal group
  mn_al <- norm_ph[caseid == cs, aliquot]
  rel_match <- if (length(mn_al) == 0L) rep(NA_real_, length(tvec)) else
                 tvec - rowMeans(Mimp[, mn_al, drop = FALSE])
  rel_rows[[i]] <- data.table(
    gene            = genes_M,
    caseid          = cs,
    prot.rel.all    = ecdf_rank(tvec, N_all),
    prot.rel.tissue = ecdf_rank(tvec, N_tis),
    prot.rel.match  = rel_match
  )
}
rel <- rbindlist(rel_rows)
setkey(rel, gene, caseid)
message("[adapter] ECDF relative-protein rows: ", nrow(rel),
        " (", uniqueN(tum_ph$aliquot), " tumour aliquots x ", length(genes_M), " genes)")

## ===========================================================================
## PART 2 -- DESCRIBEPROT structure features joined to genes via UniProt
## ===========================================================================
dp <- fread(RAW_DP, skip = 21)                 # ACC + 20 structure/annot cols
# numeric-ize the structure columns (cols 4..end); ACC/entry/name stay char
struct_cols <- names(dp)[4:ncol(dp)]
dp[, (struct_cols) := lapply(.SD, as.numeric), .SDcols = struct_cols]

uni <- as.data.frame(org.Hs.egUNIPROT)         # gene_id, uniprot_id
sym <- as.data.frame(org.Hs.egSYMBOL)          # gene_id, symbol
smap <- merge(sym, uni, by = "gene_id")        # symbol <-> uniprot (all pairs)
setDT(smap)

if (STRICT) {
  smap1 <- smap[!duplicated(gene_id)]                     # James' exact recipe
  gene2acc <- unique(smap1[, .(symbol, uniprot_id)])
} else {
  # keep pairs whose ACC exists in DESCRIBEPROT, pick smallest ACC per symbol
  cand <- smap[uniprot_id %in% dp$ACC, .(symbol, uniprot_id)]
  gene2acc <- unique(cand)[order(symbol, uniprot_id)][, .(uniprot_id = uniprot_id[1]), by = symbol]
}
setnames(dp, "ACC", "uniprot_id")
setnames(gene2acc, "symbol", "gene")
gene_struct <- merge(gene2acc, dp, by = "uniprot_id", all.x = TRUE)
# keep gene + the 18 modelling columns James references (+ ProteinName for trace)
keep_struct <- c("gene","uniprot_id", struct_cols)
gene_struct <- gene_struct[, intersect(keep_struct, names(gene_struct)), with = FALSE]
setkey(gene_struct, gene)
n_feat_genes <- gene_struct[!is.na(get(struct_cols[1])), uniqueN(gene)]
message("[adapter] DESCRIBEPROT: genes with structure features = ", n_feat_genes,
        if (STRICT) "  (STRICT/James recipe)" else "  (match-any-uniprot)")

## ===========================================================================
## PART 3 -- ADC distinct-genes list (union of tableS3 symbols)
## ===========================================================================
# tableS3: row1 = title, row2 = cancer-type codes, rows 3+ = gene symbols in a
# grid (blanks pad short columns). James: unlist -> unique -> drop "" ; the
# cancer codes never collide with real gene symbols so they're harmless, but we
# drop them explicitly for cleanliness.
# The ADC Atlas table S3 arrives either as the legacy pre-extracted CSV or as
# the published .xlsx supplement (04_auxiliary.R downloads the .xlsx). Read
# whichever we have; both share the grid layout (title row, cancer-code row,
# then a gene-symbol grid).
if (grepl("\\.xlsx?$", ADC_CSV, ignore.case = TRUE)) {
  suppressWarnings(suppressMessages(library(readxl)))
  # sheet holding table S3 (fall back to first sheet); skip the title row
  sheets <- readxl::excel_sheets(ADC_CSV)
  s3 <- sheets[grepl("S3|target", sheets, ignore.case = TRUE)][1]
  if (is.na(s3)) s3 <- sheets[1]
  adc_raw <- as.data.table(readxl::read_excel(ADC_CSV, sheet = s3, skip = 1,
                                              col_names = FALSE))
} else {
  adc_raw <- fread(ADC_CSV, header = FALSE, skip = 1, fill = TRUE)  # skip title row
}
cancer_codes <- as.character(unlist(adc_raw[1, ]))               # header row of codes
adc_vals <- unique(as.character(unlist(adc_raw[-1, ])))          # all gene cells
adc_genes <- sort(setdiff(adc_vals[nzchar(adc_vals)], cancer_codes))
message("[adapter] ADC distinct genes: ", length(adc_genes))

## ===========================================================================
## PART 4 -- Assemble the James-compatible wide CSV
## ===========================================================================
core <- as.data.table(read_parquet(file.path(DIR_TAB, "omic_table_protein_core.parquet")))

# rename our columns -> James' names
out <- core[, .(
  gene, caseid, tumor_code,
  cn,
  cn_adjusted,
  cn.capped            = pmin(cn, 5),
  rna                  = rna_norm,
  meth                 = meth_beta,
  gtex                 = gtex_tpm,
  prot                 = prot_abs,
  cytogenetic_location = cytoband,
  chromosome, arm,
  is_surface, is_secreted
)]

# attach ECDF relative columns (by gene+caseid)
out <- merge(out, rel, by = c("gene","caseid"), all.x = TRUE)
# attach structure features (by gene)
out <- merge(out, gene_struct, by = "gene", all.x = TRUE)

setcolorder(out, c("gene","caseid","tumor_code",
                   "cn","cn_adjusted","cn.capped",
                   "rna","meth","gtex","prot",
                   "prot.rel.match","prot.rel.tissue","prot.rel.all",
                   "cytogenetic_location","chromosome","arm",
                   "is_surface","is_secreted","uniprot_id", struct_cols))

## ===========================================================================
## PART 5 -- Write outputs
## ===========================================================================
fwrite(out, file.path(DIR_TAB, "str_omic_table_rebuilt.csv"))
fwrite(data.table(adc_genes), file.path(DIR_TAB, "adc_distinct_genes.csv"),
       col.names = FALSE)
write_parquet(gene_struct, file.path(DIR_TAB, "describeprot_gene_features.parquet"))

## ---- Coverage / column audit ----------------------------------------------
target_present <- out[!is.na(prot.rel.tissue), .N]
adc_in_table   <- out[gene %in% adc_genes, uniqueN(gene)]
rep_lines <- c(
  "James adapter -- coverage audit",
  paste0("rows: ", format(nrow(out), big.mark=",")),
  paste0("genes: ", uniqueN(out$gene), "   caseids: ", uniqueN(out$caseid)),
  paste0("prot.rel.tissue non-NA (target rows): ", format(target_present, big.mark=",")),
  paste0("genes with DESCRIBEPROT features: ", n_feat_genes,
         if (STRICT) " (STRICT)" else " (match-any)"),
  paste0("ADC distinct genes: ", length(adc_genes),
         "  | present in table: ", adc_in_table),
  "",
  "columns emitted:",
  paste0("  ", paste(names(out), collapse=", "))
)
writeLines(rep_lines, file.path(DIR_REP, "james_adapter_coverage.txt"))
cat(paste(rep_lines, collapse="\n"), "\n")