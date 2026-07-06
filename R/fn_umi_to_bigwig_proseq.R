#' Convert PRO-seq UMI Counts to BigWig
#'
#' @description
#' Builds a shell command to convert PRO-seq UMI counts into CPM-normalized,
#' per-strand bigwig tracks (plus strand `.ps.bw`, minus strand `.ns.bw`). Wraps
#' `inst/Rscript/umi_to_bigwig_proseq.R`. Submit with [fn_submit()].
#'
#' @param umi.counts Path to the UMI counts file from [fn_umi_collapse_proseq()] (single file).
#' @param base.name Base name / output prefix. Default NULL (derived from the counts file name).
#' @param bw.output.dir Directory for the bigwig files. Default 'db/bw/'.
#' @return data.table with columns: umi_counts, ps_bw, ns_bw, cmd, path.
#' @import data.table
#' @export
fn_umi_to_bigwig_proseq <- function(umi.counts,
                                     base.name = NULL,
                                     bw.output.dir = "db/bw/") {

  # ---- Input validation ----
  if (length(umi.counts) != 1) stop("fn_umi_to_bigwig_proseq processes one file at a time.")
  if (is.null(base.name)) base.name <- sub("\\.txt$", "", basename(umi.counts))
  if (!dir.exists(bw.output.dir)) dir.create(bw.output.dir, recursive = TRUE)

  # ---- Output paths ----
  output_prefix <- file.path(bw.output.dir, base.name)
  ps_bw <- paste0(output_prefix, ".ps.bw")
  ns_bw <- paste0(output_prefix, ".ns.bw")

  # ---- Build command ----
  script <- system.file("Rscript", "umi_to_bigwig_proseq.R", package = "filipR")
  cmd <- paste("Rscript", shQuote(script),
               shQuote(umi.counts),
               shQuote(output_prefix))

  # ---- Return data.table ----
  data.table::data.table(
    umi_counts = umi.counts,
    ps_bw = ps_bw,
    ns_bw = ns_bw,
    cmd = cmd,
    path = ps_bw
  )
}
