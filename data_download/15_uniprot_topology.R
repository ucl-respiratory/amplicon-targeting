# =============================================================================
# 15_uniprot_topology.R  --  UniProt membrane-topology annotation (surfaceome)
# -----------------------------------------------------------------------------
# Derives, from UniProt REST, the per-gene membrane topology used by the antigen
# ACCESSIBILITY filter (paper Fig 1 / Table 1). For every gene in the TCSA
# surfaceome (stage 04 `surface.genes`), we pull:
#   - transmembrane segment count      (ft_transmem)
#   - extracellular topological domains (ft_topo_dom "Extracellular")
#   - signal peptide                   (ft_signal)
#   - GPI-anchor lipidation            (ft_lipid "GPI-anchor amidated ...")
#   - subcellular location text        (cc_subcellular_location)
# and classify:
#   surface_TM_ectodomain  = >=1 TM AND >=1 extracellular topological domain
#   surface_GPI            = GPI-anchor present
#   membrane_associated_cytoplasmic_face / intracellular_membrane / intracellular
#   surface_accessible     = surface_TM_ectodomain OR surface_GPI
# (An antibody needs an extracellular epitope: a TM protein whose loops are all
#  cytoplasmic, or a peripheral membrane protein on the cytoplasmic face, is NOT
#  accessible.) Reviewed (Swiss-Prot) human entries only; one accession per gene
#  (smallest/canonical). rest.uniprot.org is queried in batches of 100 genes.
#
# Emits: OUT/tables/surface_topology_uniprot.parquet  (all surfaceome genes)
#        + a manifest row. analysis/ subsets this to its nominated gene set.
# =============================================================================

suppressWarnings(suppressMessages({
  library(data.table); library(httr); library(jsonlite); library(arrow)
}))
if (!exists("PARSED")) source(file.path("from_source", "00_config.R"))

UNIPROT_FIELDS <- "accession,gene_primary,cc_subcellular_location,ft_transmem,ft_topo_dom,ft_signal,ft_lipid"
UNIPROT_REST   <- "https://rest.uniprot.org/uniprotkb/search"

uniprot_query_batch <- function(genes) {
  q <- sprintf("(%s) AND (organism_id:9606) AND (reviewed:true)",
               paste(sprintf("(gene:%s)", genes), collapse = " OR "))
  resp <- httr::GET(UNIPROT_REST,
                    query = list(query = q, fields = UNIPROT_FIELDS,
                                 format = "json", size = 500),
                    httr::timeout(120))
  httr::stop_for_status(resp)
  fromJSON(rawToChar(resp$content), simplifyVector = FALSE)$results
}

classify_entry <- function(res) {
  acc  <- res$primaryAccession
  gene <- tryCatch(res$genes[[1]]$geneName$value, error = function(e) NA_character_)
  feats <- res$features
  ftype <- vapply(feats, function(f) f$type %||% "", "")
  n_tm  <- sum(ftype == "Transmembrane")
  topo  <- feats[ftype == "Topological domain"]
  n_extra <- sum(vapply(topo, function(f)
    grepl("Extracellular", f$description %||% "", ignore.case = TRUE), logical(1)))
  has_signal <- any(ftype == "Signal")
  lipids <- feats[ftype == "Lipidation"]
  gpi <- any(vapply(lipids, function(f)
    grepl("GPI-anchor", f$description %||% "", ignore.case = TRUE), logical(1)))
  subc <- tryCatch(paste(vapply(res$comments[[1]]$subcellularLocations,
            function(s) s$location$value, ""), collapse = ";"),
            error = function(e) NA_character_)
  data.table(gene = gene, acc_uniprot = acc, n_tm = n_tm,
             n_extracellular_topo = n_extra, has_signal = has_signal,
             gpi = gpi, subcell = subc %||% NA_character_)
}
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

build_topology <- function() {
  out <- file.path(DIR_TAB, "surface_topology_uniprot.parquet")
  if (file.exists(out) && file.info(out)$size > 0) {
    message("[topology] cached: ", basename(out)); return(invisible(out))
  }
  stopifnot(file.exists(INPUTS$surface))
  e <- new.env(); load(INPUTS$surface, e)
  genes <- unique(stats::na.omit(e$surface.genes))
  message("[topology] querying UniProt for ", length(genes), " surfaceome genes")
  batches <- split(genes, ceiling(seq_along(genes) / 100))
  rows <- list()
  for (i in seq_along(batches)) {
    res <- tryCatch(uniprot_query_batch(batches[[i]]),
                    error = function(err) { message("  batch ", i, " failed: ",
                                                     conditionMessage(err)); list() })
    for (r in res) rows[[length(rows) + 1L]] <- classify_entry(r)
    message("  batch ", i, "/", length(batches), " -> ", length(rows), " entries")
    Sys.sleep(0.2)
  }
  dt <- rbindlist(rows, fill = TRUE)
  # one canonical accession per gene: shortest accession string (canonical form)
  dt <- dt[order(gene, nchar(acc_uniprot), acc_uniprot)][!duplicated(gene)]
  dt[, surface_class := fifelse(n_tm >= 1 & n_extracellular_topo >= 1, "surface_TM_ectodomain",
                        fifelse(gpi, "surface_GPI",
                        fifelse(grepl("membrane", subcell, ignore.case = TRUE) & n_tm >= 1,
                                "intracellular_membrane",
                        fifelse(grepl("membrane", subcell, ignore.case = TRUE),
                                "membrane_associated_cytoplasmic_face", "intracellular"))))]
  dt[, surface_accessible := surface_class %in% c("surface_TM_ectodomain", "surface_GPI")]
  arrow::write_parquet(dt, out)
  record_source("15_uniprot_topology", "UniProt REST membrane topology (surfaceome)",
                UNIPROT_REST, out,
                note = sprintf("%d genes; %d surface-accessible", nrow(dt), sum(dt$surface_accessible)))
  message("[topology] wrote ", nrow(dt), " genes; ",
          sum(dt$surface_accessible), " surface-accessible")
  invisible(out)
}

if (identical(Sys.getenv("FS_TOPOLOGY_RUN"), "1") || !interactive()) {
  build_topology()
} else {
  message("[topology] function defined. Set FS_TOPOLOGY_RUN=1 to execute.")
}
