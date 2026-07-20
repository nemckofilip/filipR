#' Annotate per-observation edit calls with gene / mature-transcript position
#'
#' @description
#' Builds a shell command to map the output of [fn_edit_calls()] onto each
#' gene's mature-transcript coordinate, so downstream queries (e.g. taking reads
#' with at least N edits and asking where those edits fall along the gene) need
#' no exon arithmetic. The heavy join runs in `inst/Rscript/annotate_edit_calls.R`.
#' Submit with [fn_submit()].
#'
#' The `gene.model` RDS is a one-time exon-union transcript model (built by the
#' calling subscript); keeping it out of [fn_edit_calls()] means the transcript
#' definition can be re-cut without re-running the BAM walk.
#'
#' @param edit.calls Path to a `*_edit_calls.tsv.gz` from [fn_edit_calls()].
#' @param gene.model Path to the exon-union gene-model RDS.
#' @param base.name Base name for the output file.
#' @param output.dir Output directory. Default `"db/alignment_stats/edit_calls_annot/"`.
#'
#' @return A one-row `data.table` with columns: `edit_calls`, `output`, `cmd`, `path`.
#' @export
fn_annotate_edit_calls <- function(edit.calls,
                                   gene.model,
                                   base.name,
                                   output.dir = "db/alignment_stats/edit_calls_annot/") {

  if (!file.exists(gene.model))
    stop("gene.model does not exist: ", gene.model)
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  script <- system.file("Rscript", "annotate_edit_calls.R", package = "filipR")
  out    <- file.path(output.dir, paste0(base.name, "_edit_calls_annot.rds"))

  cmd <- paste(
    "Rscript", shQuote(script),
    shQuote(edit.calls),
    shQuote(gene.model),
    shQuote(out)
  )

  data.table::data.table(edit_calls = edit.calls, output = out,
                         cmd = cmd, path = out)
}
