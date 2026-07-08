# =============================================================================
# 03_gdc_genomics.R  --  Copy number, RNA, methylation, SNV from the GDC
# -----------------------------------------------------------------------------
# Downloads the four genomic layers for the CPTAC-3 cohort from the NCI Genomic
# Data Commons via the GenomicDataCommons Bioconductor client, and parses each
# to the object shapes the assembly chain consumes.
#
# MODERNIZATION vs the original download.data.R (documented drift):
#   * RNA workflow: original used 'HTSeq - Counts', which GDC RETIRED. The
#     current bulk RNA workflow is 'STAR - Counts' (GENCODE v36). We read the
#     'unstranded' column (the HTSeq-count analogue) from the augmented
#     star_gene_counts.tsv, skipping the 4 STAR summary rows (N_unmapped, etc.).
#     Absolute counts differ from the 2022 HTSeq build (different aligner,
#     gene model) -> DESeq2-normalized values will not match to the digit.
#   * SNV workflow: original 'masked_somatic_mutation'; current data_type is
#     'Masked Somatic Mutation' produced by 'Aliquot Ensemble Somatic Variant
#     Merging and Masking'. MAF columns are unchanged (Hugo_Symbol, ...).
#   * CN: 'Gene Level Copy Number' (AscatNGS). The gene-level file now carries a
#     header (gene_id,gene_name,chromosome,start,end,copy_number,min,max); we
#     read gene_name + copy_number and max-aggregate per gene, as the original.
#   * Methylation: 'Methylation Beta Value' (SeSAMe level-3 betas, EPIC array).
#     Two-column (probe, beta) files, gene mapping via ChAMP probe.features.epic.
#
# Sample->case mapping is taken from GDC file metadata (cases.samples
# submitter_id + sample_type). caseid = first 9 chars of the tumour submitter_id
# (C3L-xxxxx / C3N-xxxxx), matching the original substr(.,1,9) convention.
#
# Emits (object names/shapes matching the cached pipeline contract):
#   PARSED/cptac.cndata.RData    cndata (genes x caseid), cnpheno
#   PARSED/cptac.rnadata.RData   rnadata, rnadata.norm (DESeq2), rnapheno
#   PARSED/methdata.genes.RData  methdata.genes (gene-body mean beta), methpheno
#   PARSED/methdata.tss.RData    methdata.tss (TSS200/1500 mean beta), methpheno
#   PARSED/cptac.snvdata.RData   snvdata (long MAF subset), snvpheno
#
# Requires (heavy Bioconductor deps): GenomicDataCommons, DESeq2, maftools, ChAMP.
# =============================================================================

suppressWarnings(suppressMessages({
  library(data.table)
  library(GenomicDataCommons)
}))
if (!exists("PARSED")) source(file.path("from_source", "00_config.R"))
if (!exists("genemap")) {
  if (file.exists(INPUTS$genemap)) load(INPUTS$genemap) else
    stop("genemap not found; run 01_annotation.R first")
}
# Proteome cohort defines which cases we keep (matches original cptac.pheno gate).
if (!file.exists(INPUTS$proteome))
  stop("proteome parsed object missing; run 02_proteome_pdc.R first")
.penv <- new.env(); load(INPUTS$proteome, .penv)
PROT_CASES  <- unique(.penv$cptac.pheno$case_id)
PROT_GENES  <- rownames(.penv$cptac.data)

gdc_set_cache(GDC_CACHE)

## ---- helper: files() query -> metadata data.table with case + sample_type ---
# Returns one row per file: file_id, file_name, tumor_sid, normal_sid, caseid.
gdc_file_meta <- function(data_type, workflow_type = NULL) {
  q <- files() |>
    GenomicDataCommons::filter(cases.project.project_id == GDC_PROJECT) |>
    GenomicDataCommons::filter(data_type == !!data_type)
  if (!is.null(workflow_type))
    q <- q |> GenomicDataCommons::filter(analysis.workflow_type == !!workflow_type)
  x <- q |>
    GenomicDataCommons::expand(c("cases", "cases.samples")) |>
    GenomicDataCommons::results(size = 100000)
  n <- length(x$file_id)
  meta <- data.table(file_id = x$file_id, file_name = x$file_name,
                     tumor_sid = NA_character_, normal_sid = NA_character_)
  for (i in seq_len(n)) {
    samps <- tryCatch(x$cases[[i]]$samples[[1]], error = function(e) NULL)
    if (is.null(samps)) next
    st  <- samps$sample_type; sid <- samps$submitter_id
    ti  <- grep("Tumor",  st); ni <- grep("Normal", st)
    if (length(ti)) meta$tumor_sid[i]  <- sid[ti[1]]
    if (length(ni)) meta$normal_sid[i] <- sid[ni[1]]
  }
  meta[, caseid := substr(fifelse(is.na(tumor_sid), normal_sid, tumor_sid), 1, 9)]
  meta[]
}

find_cached <- function(file_name) {
  f <- list.files(GDC_CACHE, pattern = file_name, recursive = TRUE, full.names = TRUE)
  if (length(f) == 0) NA_character_ else f[1]
}

## ===========================================================================
## COPY NUMBER  (AscatNGS gene-level; copy_number, max-aggregate per gene)
## ===========================================================================
build_cn <- function() {
  message("[gdc/cn] querying Gene Level Copy Number (AscatNGS)")
  cnpheno <- gdc_file_meta("Gene Level Copy Number", "AscatNGS")
  cnpheno <- cnpheno[!is.na(tumor_sid)]                      # tumour CN only
  cnpheno <- cnpheno[!duplicated(caseid)]
  cnpheno <- cnpheno[caseid %in% PROT_CASES]
  manifest <- gdcdata(cnpheno$file_id, progress = FALSE)     # download to cache
  cols <- list()
  for (i in seq_len(nrow(cnpheno))) {
    f <- find_cached(cnpheno$file_name[i]); if (is.na(f)) next
    x <- fread(f, sep = "\t", header = TRUE)
    x <- x[gene_name %in% PROT_GENES, .(gene = gene_name, cn = copy_number)]
    x <- x[, .(cn = max(cn, na.rm = TRUE)), by = gene]        # max-aggregate
    v <- setNames(x$cn, x$gene); cols[[cnpheno$caseid[i]]] <- v
  }
  genes <- sort(unique(unlist(lapply(cols, names))))
  cndata <- as.data.frame(sapply(cols, function(v) v[genes]))
  rownames(cndata) <- genes
  colnames(cndata) <- names(cols)
  # standardize cnpheno columns to the original (file_id,filename,tumorid,normalid,caseid)
  cnpheno <- cnpheno[caseid %in% colnames(cndata),
                     .(file_id, filename = file_name, tumorid = tumor_sid,
                       normalid = normal_sid, caseid)]
  save(cndata, cnpheno, file = INPUTS$cn)
  record_source("03_gdc_genomics", "GDC CPTAC-3 Gene Level Copy Number (AscatNGS)",
                "https://api.gdc.cancer.gov (GenomicDataCommons)", INPUTS$cn,
                note = sprintf("cndata %d genes x %d caseids", nrow(cndata), ncol(cndata)))
  message("[gdc/cn] cndata: ", nrow(cndata), " genes x ", ncol(cndata), " caseids")
}

## ===========================================================================
## RNA  (STAR - Counts; 'unstranded' column; ENSEMBL->symbol sum; DESeq2 norm)
## ===========================================================================
parse_star_counts <- function(path) {
  # Skip the '# gene-model' comment line + read header; drop 4 STAR summary rows.
  dt <- fread(path, sep = "\t", skip = 1, header = TRUE)
  dt <- dt[!gene_id %in% c("N_unmapped","N_multimapping","N_noFeature","N_ambiguous")]
  dt[, .(gene_id, unstranded)]
}
build_rna <- function() {
  message("[gdc/rna] querying Gene Expression Quantification (STAR - Counts)")
  rp <- gdc_file_meta("Gene Expression Quantification", GDC_RNA_WORKFLOW)
  rp <- rp[grepl("augmented_star_gene_counts", file_name)]   # the counts file
  rp[, caseid := substr(fifelse(is.na(tumor_sid), normal_sid, tumor_sid), 1, 9)]
  rp <- rp[!is.na(tumor_sid)]                                # tumour RNA
  rp <- rp[!duplicated(caseid)][caseid %in% PROT_CASES]
  gdcdata(rp$file_id, progress = FALSE)
  cols <- list(); ens_ref <- NULL
  for (i in seq_len(nrow(rp))) {
    f <- find_cached(rp$file_name[i]); if (is.na(f)) next
    d <- parse_star_counts(f)
    if (is.null(ens_ref)) ens_ref <- d$gene_id
    cols[[rp$caseid[i]]] <- setNames(d$unstranded, d$gene_id)[ens_ref]
  }
  rnadata <- as.data.frame(do.call(cbind, cols)); rownames(rnadata) <- ens_ref
  # ENSEMBL (strip version) -> symbol via genemap; sum transcripts per symbol.
  ens <- sub("\\..*$", "", rownames(rnadata))
  sym <- genemap$symbol[match(ens, genemap$ensembl_id)]
  sel <- which(!is.na(sym))
  agg <- rowsum(rnadata[sel, , drop = FALSE], group = sym[sel])
  rnadata <- agg
  rnapheno <- rp[caseid %in% colnames(rnadata),
                 .(file_id, filename = file_name, tumorid = tumor_sid,
                   normalid = normal_sid, caseid)]
  # DESeq2 size-factor normalization (design ~1), as original.
  suppressWarnings(suppressMessages(library(DESeq2)))
  dds <- DESeqDataSetFromMatrix(round(as.matrix(rnadata)),
                                DataFrame(rnapheno[match(colnames(rnadata), caseid)]),
                                design = ~1)
  dds <- estimateSizeFactors(dds)
  rnadata.norm <- as.data.frame(counts(dds, normalized = TRUE))
  save(rnadata, rnadata.norm, rnapheno, file = INPUTS$rna)
  record_source("03_gdc_genomics",
                paste0("GDC CPTAC-3 Gene Expression (", GDC_RNA_WORKFLOW, ")"),
                "https://api.gdc.cancer.gov (GenomicDataCommons)", INPUTS$rna,
                note = sprintf("rnadata %d genes x %d caseids; MODERNIZED from HTSeq-Counts",
                               nrow(rnadata), ncol(rnadata)))
  message("[gdc/rna] rnadata: ", nrow(rnadata), " genes x ", ncol(rnadata), " caseids")
}

## ===========================================================================
## SNV  (Masked Somatic Mutation; MAF -> long subset via maftools)
## ===========================================================================
build_snv <- function() {
  message("[gdc/snv] querying Masked Somatic Mutation")
  sp <- gdc_file_meta("Masked Somatic Mutation",
                      "Aliquot Ensemble Somatic Variant Merging and Masking")
  sp <- sp[!is.na(tumor_sid)][!duplicated(caseid)][caseid %in% PROT_CASES]
  gdcdata(sp$file_id, progress = FALSE)
  suppressWarnings(suppressMessages(library(maftools)))
  rows <- list()
  for (i in seq_len(nrow(sp))) {
    f <- find_cached(sp$file_name[i]); if (is.na(f)) next
    res <- tryCatch({
      m <- read.maf(f, verbose = FALSE)
      data.table(filename = sp$file_name[i], gene = m@data$Hugo_Symbol,
                 chr = m@data$Chromosome, start = m@data$Start_Position,
                 end = m@data$End_Position, type = m@data$Variant_Type,
                 classification = m@data$Variant_Classification,
                 ref = m@data$Reference_Allele, alt1 = m@data$Tumor_Seq_Allele1,
                 alt2 = m@data$Tumor_Seq_Allele2, HGVSp = m@data$HGVSp_Short)
    }, error = function(e) NULL)
    if (!is.null(res)) rows[[length(rows)+1L]] <- res
  }
  snvdata <- rbindlist(rows)
  snvdata[, caseid := sp$caseid[match(filename, sp$file_name)]]
  snvpheno <- sp[caseid %in% snvdata$caseid,
                 .(file_id, filename = file_name, tumorid = tumor_sid,
                   normalid = normal_sid, caseid)]
  save(snvdata, snvpheno, file = INPUTS$snv)
  record_source("03_gdc_genomics", "GDC CPTAC-3 Masked Somatic Mutation",
                "https://api.gdc.cancer.gov (GenomicDataCommons)", INPUTS$snv,
                note = sprintf("snvdata %d variant rows x %d caseids",
                               nrow(snvdata), uniqueN(snvdata$caseid)))
  message("[gdc/snv] snvdata: ", nrow(snvdata), " variants, ",
          uniqueN(snvdata$caseid), " caseids")
}

## ===========================================================================
## METHYLATION  (SeSAMe beta; ChAMP impute; probe.features gene/TSS mean beta)
## ===========================================================================
build_meth <- function() {
  message("[gdc/meth] querying Methylation Beta Value (SeSAMe, EPIC)")
  mp <- gdc_file_meta("Methylation Beta Value")
  mp[, caseid := substr(fifelse(is.na(tumor_sid), normal_sid, tumor_sid), 1, 9)]
  mp <- mp[!is.na(tumor_sid)][!duplicated(caseid)][caseid %in% PROT_CASES]
  gdcdata(mp$file_id, progress = FALSE)
  cols <- list(); probe_ref <- NULL
  for (i in seq_len(nrow(mp))) {
    f <- find_cached(mp$file_name[i]); if (is.na(f)) next
    d <- fread(f, sep = "\t", header = FALSE)          # V1 probe, V2 beta
    if (is.null(probe_ref)) probe_ref <- d$V1
    cols[[mp$caseid[i]]] <- setNames(d$V2, d$V1)[probe_ref]
  }
  methmat <- as.matrix(as.data.frame(cols)); rownames(methmat) <- probe_ref
  # ChAMP impute + probe->gene mapping (EPIC), as the original.
  suppressWarnings(suppressMessages(library(ChAMP)))
  methpheno <- mp[caseid %in% colnames(methmat),
                  .(file_id, filename = file_name, tumorid = tumor_sid,
                    normalid = normal_sid, caseid)]
  imp <- champ.impute(methmat, as.data.frame(methpheno),
                      ProbeCutoff = 0.2, SampleCutoff = 0.3)
  methdata.imp <- as.data.frame(imp$beta)
  data("probe.features.epic", package = "ChAMP", envir = environment())
  pf <- get("probe.features")
  # gene-body mean beta
  g   <- as.character(pf$gene[match(rownames(methdata.imp), rownames(pf))])
  sel <- which(g != "" & !is.na(g))
  methdata.genes <- rowsum(methdata.imp[sel, , drop = FALSE], g[sel]) /
                    as.vector(table(g[sel])[unique(g[sel])])   # mean
  methdata.genes <- as.data.frame(
    aggregate(methdata.imp[sel, ], by = list(g[sel]), FUN = mean))
  rownames(methdata.genes) <- methdata.genes$Group.1; methdata.genes$Group.1 <- NULL
  save(methdata.genes, methpheno, file = INPUTS$meth_genes)
  # TSS (TSS200/TSS1500) mean beta
  probes <- rownames(pf)[pf$feature %in% c("TSS200","TSS1500") & pf$gene != ""]
  probes <- intersect(probes, rownames(methdata.imp))
  tss <- aggregate(methdata.imp[probes, ], by = list(pf[probes, "gene"]), FUN = mean)
  rownames(tss) <- tss$Group.1; tss$Group.1 <- NULL
  methdata.tss <- tss
  save(methdata.tss, methpheno, file = INPUTS$meth_tss)
  record_source("03_gdc_genomics", "GDC CPTAC-3 Methylation Beta Value (SeSAMe EPIC)",
                "https://api.gdc.cancer.gov (GenomicDataCommons)", INPUTS$meth_genes,
                note = sprintf("methdata.genes %d genes x %d caseids",
                               nrow(methdata.genes), ncol(methdata.genes)))
  message("[gdc/meth] methdata.genes: ", nrow(methdata.genes), " genes x ",
          ncol(methdata.genes), " caseids")
}

## ---- Run all four layers (guarded so a partial run is resumable) ------------
if (identical(Sys.getenv("FS_GDC_RUN"), "1") || !interactive()) {
  build_cn(); build_rna(); build_snv(); build_meth()
} else {
  message("[gdc] functions defined. Set FS_GDC_RUN=1 to execute all four layers.")
}
