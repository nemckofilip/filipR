#' Metagene editing profile from annotated edit calls
#'
#' @description Builds a shell command that bins per-observation edit calls by
#'   position along the mature transcript and tallies edited / total target
#'   bases per bin, producing two positional editing profiles: editing rate vs
#'   absolute distance from the 3' end, and editing rate vs fractional 5'->3'
#'   position (metagene). Heavy lifting runs in `inst/Rscript/metagene_editing.R`.
#'   Submit with [fn_submit()].
#' @param edit.calls.annot Character. Path to an annotated edit-calls RDS
#'   (output of [fn_annotate_edit_calls()]).
#' @param base.name Base name for output files.
#' @param output.dir Output directory. Default 'db/alignment_stats/metagene_editing/'.
#' @param min.gene.len Integer. Minimum gene length (nt) to include. Default 1000.
#' @param bin.size Integer. Bin size (nt) for the 3'-distance profile. Default 100.
#' @param max.dist Integer. Maximum distance from 3' end (nt) covered by the
#'   distance profile. Default 20000.
#' @param n.frac.bins Integer. Number of bins for the fractional metagene. Default 20.
#' @return data.table with columns: output, cmd, path. Output TSV is split by
#'   gene-length band (len_band).
#' @import data.table
#' @export
fn_metagene_editing <- function(edit.calls.annot,
                                 base.name,
                                 output.dir   = "db/alignment_stats/metagene_editing/",
                                 min.gene.len = 1000,
                                 bin.size     = 100,
                                 max.dist     = 20000,
                                 n.frac.bins  = 20) {

  if (length(edit.calls.annot) != 1)
    stop("fn_metagene_editing processes one sample at a time.")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  script <- system.file("Rscript", "metagene_editing.R", package = "filipR")
  out    <- file.path(output.dir, paste0(base.name, "_metagene_editing.tsv"))

  cmd <- paste("Rscript", shQuote(script),
               shQuote(edit.calls.annot), shQuote(out),
               min.gene.len, bin.size, max.dist, n.frac.bins)

  data.table::data.table(output = out, cmd = cmd, path = out)
}
