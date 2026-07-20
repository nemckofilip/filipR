#' Per-read edit-count distribution split by 3'-proximity
#'
#' @description Builds a shell command that, from an annotated edit-calls RDS,
#'   tallies the per-read on-target edit-count distribution for all genic reads
#'   and for reads overlapping the 3'-proximal window (closest target base within
#'   `dist.thr` nt of the 3' end). Used to ask whether restricting to
#'   3'-proximal reads shifts the edits-per-read survival. Heavy lifting runs in
#'   `inst/Rscript/edits_per_read_3p.R`. Submit with [fn_submit()].
#' @param edit.calls.annot Character. Path to an annotated edit-calls RDS
#'   (output of [fn_annotate_edit_calls()]).
#' @param base.name Base name for output files.
#' @param output.dir Output directory. Default 'db/alignment_stats/edits_per_read_3p/'.
#' @param dist.thr Integer. 3'-proximal window (nt): a read is "proximal" if its
#'   closest target base to the 3' end is within this distance. Default 500.
#' @return data.table with columns: output, cmd, path. Output TSV columns:
#'   region ("all" | "prox"), n_edited, n_reads.
#' @import data.table
#' @export
fn_edits_per_read_3p <- function(edit.calls.annot,
                                 base.name,
                                 output.dir = "db/alignment_stats/edits_per_read_3p/",
                                 dist.thr   = 500) {

  if (length(edit.calls.annot) != 1)
    stop("fn_edits_per_read_3p processes one sample at a time.")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  script <- system.file("Rscript", "edits_per_read_3p.R", package = "filipR")
  out    <- file.path(output.dir, paste0(base.name, "_edits_per_read_3p.tsv"))

  cmd <- paste("Rscript", shQuote(script),
               shQuote(edit.calls.annot), shQuote(out), dist.thr)

  data.table::data.table(output = out, cmd = cmd, path = out)
}
