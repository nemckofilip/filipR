#' Test sgRNA/Gene Enrichment Using MAGeCK test
#'
#' @description
#' Builds a shell command for `mageck test`, comparing treatment vs control
#' samples from a MAGeCK count table and producing gene- and sgRNA-level
#' enrichment/depletion statistics (RRA). Submit with [fn_submit()].
#'
#' @param count.table Path to a MAGeCK count table (`.count.txt`).
#' @param treatment Character vector of treatment sample labels (must match the
#'   count-table header).
#' @param control Character vector of control sample labels.
#' @param base.name Output prefix name (MAGeCK `-n`).
#' @param control.sgrna Optional path to a control sgRNA ID list for the null
#'   distribution (`--control-sgrna`). Default NULL.
#' @param paired Logical. Paired comparison (treatment/control matched by order,
#'   e.g. rep1 vs rep1)? Default FALSE.
#' @param norm.method Normalisation: "none", "median", "total", or "control".
#'   Default "control" (normalise by the control sgRNAs).
#' @param output.dir Output directory. Default "db/counts/".
#' @return data.table with columns: gene_summary, sgrna_summary, cmd, path.
#' @import data.table
#' @export
fn_mageck_test <- function(count.table,
                           treatment,
                           control,
                           base.name,
                           control.sgrna = NULL,
                           paired = FALSE,
                           norm.method = "control",
                           output.dir = "db/counts/") {

  # ---- Input validation ----
  if (length(count.table) != 1) stop("Provide exactly one count.table.")
  if (paired && length(treatment) != length(control))
    stop("Paired mode requires equal numbers of treatment and control samples.")
  if (!norm.method %in% c("none", "median", "total", "control"))
    stop("norm.method must be one of: none, median, total, control.")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  # ---- Output paths ----
  prefix      <- file.path(output.dir, base.name)
  gene_out    <- paste0(prefix, ".gene_summary.txt")
  sgrna_out   <- paste0(prefix, ".sgrna_summary.txt")

  # ---- Build command ----
  cmd <- paste(
    "mageck test",
    "-k", shQuote(count.table),
    "-t", shQuote(paste(treatment, collapse = ",")),
    "-c", shQuote(paste(control, collapse = ",")),
    "-n", shQuote(prefix),
    "--norm-method", norm.method
  )
  if (paired) cmd <- paste(cmd, "--paired")
  if (!is.null(control.sgrna))
    cmd <- paste(cmd, "--control-sgrna", shQuote(control.sgrna))

  # ---- Return data.table ----
  data.table::data.table(
    gene_summary = gene_out,
    sgrna_summary = sgrna_out,
    cmd = cmd,
    path = gene_out
  )
}
