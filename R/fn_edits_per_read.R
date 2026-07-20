#' Per-read edit-count distribution from an edit_calls table
#'
#' @description
#' Builds a shell command to summarise, from the output of [fn_edit_calls()],
#' how many on-target edits each read carries — i.e. for every read, the number
#' of its target bases that are edited — and returns the distribution across
#' reads. Downstream this gives the fraction of reads with at least N edits
#' (survival curve) or a per-read edit-count histogram. The count runs in
#' `inst/Rscript/edits_per_read_summary.R`. Submit with [fn_submit()].
#'
#' @param edit.calls Path to a `*_edit_calls.tsv.gz` from [fn_edit_calls()].
#' @param base.name Base name for the output file.
#' @param output.dir Output directory. Default `"db/alignment_stats/edits_per_read/"`.
#'
#' @return A one-row `data.table` with columns: `edit_calls`, `output`, `cmd`, `path`.
#'   The output RDS holds `data.table(n_edited, n_reads)` with the total read
#'   count in attribute `total_reads`.
#' @export
fn_edits_per_read <- function(edit.calls,
                              base.name,
                              output.dir = "db/alignment_stats/edits_per_read/") {

  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  script <- system.file("Rscript", "edits_per_read_summary.R", package = "filipR")
  out    <- file.path(output.dir, paste0(base.name, "_edits_per_read.rds"))

  cmd <- paste("Rscript", shQuote(script), shQuote(edit.calls), shQuote(out))

  data.table::data.table(edit_calls = edit.calls, output = out,
                         cmd = cmd, path = out)
}
