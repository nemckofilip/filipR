#' Collapse UMIs from a PRO-seq BAM
#'
#' @description
#' Builds a shell command to collapse UMIs from an aligned PRO-seq BAM and count
#' unique molecules per single-nucleotide genomic position. Reads are resized to
#' their 3' end and (for PRO-seq) the strand is flipped, then UMIs within 1
#' mismatch are collapsed per position. Wraps `inst/Rscript/umi_collapse_proseq.R`.
#' Submit with [fn_submit()].
#'
#' @param bam Path to the input BAM file (single file).
#' @param base.name Base name for output files. Default NULL (derived from the BAM name).
#' @param counts.output.dir Directory for the UMI counts file. Default 'db/counts/'.
#' @param stats.output.dir Directory for the UMI stats file. Default 'db/counts_stats/'.
#' @param flip.strand Logical. Flip the read strand? TRUE for PRO-seq (default).
#' @param umi.length Integer. UMI length in nt. Default 8.
#' @return data.table with columns: bam, umi_counts, umi_stats, cmd, path.
#' @import data.table
#' @export
fn_umi_collapse_proseq <- function(bam,
                                    base.name = NULL,
                                    counts.output.dir = "db/counts/",
                                    stats.output.dir = "db/counts_stats/",
                                    flip.strand = TRUE,
                                    umi.length = 8) {

  # ---- Input validation ----
  if (length(bam) != 1) stop("fn_umi_collapse_proseq processes one BAM at a time.")
  if (is.null(base.name)) base.name <- sub("\\.bam$", "", basename(bam))
  if (!dir.exists(counts.output.dir)) dir.create(counts.output.dir, recursive = TRUE)
  if (!dir.exists(stats.output.dir)) dir.create(stats.output.dir, recursive = TRUE)

  # ---- Output paths ----
  counts_out <- file.path(counts.output.dir, paste0(base.name, "_UMI_counts.txt"))
  stats_out <- file.path(stats.output.dir, paste0(base.name, "_UMI_stats.txt"))

  # ---- Build command ----
  script <- system.file("Rscript", "umi_collapse_proseq.R", package = "filipR")
  cmd <- paste("Rscript", shQuote(script),
               shQuote(bam),
               shQuote(counts_out),
               shQuote(stats_out),
               flip.strand,
               umi.length)

  # ---- Return data.table ----
  data.table::data.table(
    bam = bam,
    umi_counts = counts_out,
    umi_stats = stats_out,
    cmd = cmd,
    path = counts_out
  )
}
