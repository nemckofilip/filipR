#' Read a minimap2 PAF File
#'
#' @description
#' Reads the 12 mandatory PAF columns produced by [fn_minimap2_map()] into a
#' `data.table`, with the column names spelled out so callers do not have to.
#' Trailing key:type:value tags are ignored. Returns a correctly typed empty
#' table for an empty file, so callers can skip a `nrow()` guard.
#'
#' @param path Character. Path to a `.paf` file.
#' @param select Character or `NULL`. Subset of columns to keep, e.g.
#'   `c("qname", "tname", "nmatch")`. `NULL` (default) keeps all 12.
#' @param dedup Logical. Drop duplicate query-target pairs, which minimap2
#'   emits when a query aligns to the same target in several blocks.
#'   Default: `TRUE`.
#'
#' @return A `data.table` with columns `qname`, `qlen`, `qstart`, `qend`,
#'   `strand`, `tname`, `tlen`, `tstart`, `tend`, `nmatch`, `alen`, `mapq`
#'   (or the subset given by `select`).
#'
#' @examples
#' \dontrun{
#' p <- fn_read_paf("db/paf/ppLS1_mapont.paf")
#'
#' # best target per query, with the runner-up as an ambiguity margin
#' p <- fn_read_paf("db/paf/ppLS1_mapont.paf",
#'                  select = c("qname", "tname", "nmatch"))
#' p[order(qname, -nmatch),
#'   .(hit = tname[1], margin = if (.N > 1) nmatch[2] / nmatch[1] else 0),
#'   by = qname]
#' }
#' @import data.table
#' @export
fn_read_paf <- function(path, select = NULL, dedup = TRUE) {

  # ---- Input validation ----
  if (length(path) != 1) stop("fn_read_paf reads one file at a time.")
  if (!file.exists(path)) stop("file not found: ", path)

  paf_cols <- c("qname", "qlen", "qstart", "qend", "strand",
                "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq")
  paf_type <- c("character", "integer", "integer", "integer", "character",
                "character", "integer", "integer", "integer", "integer",
                "integer", "integer")

  keep <- if (is.null(select)) paf_cols else select
  if (!all(keep %in% paf_cols))
    stop("select must be a subset of: ", paste(paf_cols, collapse = ", "))
  idx <- match(keep, paf_cols)

  # ---- Empty file: return the right shape rather than NULL ----
  if (file.size(path) == 0) {
    empty <- lapply(paf_type[idx], function(ty) vector(ty, 0L))
    return(data.table::setnames(data.table::as.data.table(empty), keep))
  }

  # fill = TRUE because the optional tag columns vary in number per row
  p <- data.table::fread(path, sep = "\t", header = FALSE, fill = TRUE,
                         select = idx, col.names = keep)

  if (dedup && all(c("qname", "tname") %in% keep))
    p <- unique(p, by = c("qname", "tname"))
  p[]
}
