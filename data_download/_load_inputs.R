# =============================================================================
# _load_inputs.R  --  Load all parsed inputs (from_source) into IN + manifest
# -----------------------------------------------------------------------------
# Loads every parsed .RData into a single environment `IN` and writes
# out/tables/input_manifest.csv recording, for each source file and object:
#   source file, file size, file mtime, md5(file), object name, class, dim.
# This is the provenance anchor for the whole pipeline.
# =============================================================================

suppressWarnings(suppressMessages({
  library(data.table)
  library(tools)   # md5sum
}))

if (!exists("INPUTS")) source(file.path("from_source", "00_config.R"))

IN <- new.env(parent = emptyenv())

manifest_rows <- list()

load_one <- function(tag, path) {
  if (!file.exists(path)) stop("Missing input for '", tag, "': ", path)
  fsz  <- file.info(path)$size
  fmt  <- as.character(file.info(path)$mtime)
  fmd5 <- unname(tools::md5sum(path))
  tmp  <- new.env(parent = emptyenv())
  loaded <- load(path, envir = tmp)          # names of objects loaded
  for (nm in loaded) {
    obj <- get(nm, envir = tmp)
    assign(nm, obj, envir = IN)              # hoist into shared IN env
    dim_str <- if (!is.null(dim(obj))) paste(dim(obj), collapse = "x")
               else paste0("len=", length(obj))
    manifest_rows[[length(manifest_rows) + 1L]] <<- data.table(
      tag        = tag,
      source     = basename(path),
      file_bytes = fsz,
      file_mtime = fmt,
      file_md5   = fmd5,
      object     = nm,
      class      = class(obj)[1],
      dim        = dim_str
    )
  }
  invisible(loaded)
}

for (tag in names(INPUTS)) load_one(tag, INPUTS[[tag]])

manifest <- rbindlist(manifest_rows)
fwrite(manifest, file.path(DIR_TAB, "input_manifest.csv"))

message("[load] objects loaded into IN: ", paste(ls(IN), collapse = ", "))
message("[load] manifest rows: ", nrow(manifest))
print(manifest[, .(tag, object, class, dim)])
