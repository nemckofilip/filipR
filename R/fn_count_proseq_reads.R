#' Count PRO-seq Reads Over Features
#'
#' @description
#' Builds a shell command to count PRO-seq UMI reads over annotated features
#' (strand-specific), optionally subtracting reads in blacklisted regions. Wraps
#' `inst/Rscript/count_proseq_reads.R`. Submit with [fn_submit()].
#'
#' @param umi.count.file Path to the UMI counts file from [fn_umi_collapse_proseq()] (single .txt).
#' @param annotation.file Path to a `.rds` annotation (a data.table with columns
#'   seqnames, start, end, strand, cluster.id). One row per feature (e.g. per gene).
#' @param blacklist.file Optional `.rds` of blacklisted regions to subtract. Default NULL.
#' @param feature Feature label used in the output name. Default NULL (derived from annotation file name).
#' @param output.prefix Output prefix. Default NULL (derived from the UMI counts file name).
#' @param count.tables.output.dir Directory for the counts table. Default 'db/count_tables/'.
#' @return data.table with columns: umi_count_file, count_table, cmd, path.
#' @import data.table
#' @export
fn_count_proseq_reads <- function(umi.count.file,
                                   annotation.file,
                                   blacklist.file = NULL,
                                   feature = NULL,
                                   output.prefix = NULL,
                                   count.tables.output.dir = "db/count_tables/") {

  # ---- Input validation ----
  if (length(umi.count.file) != 1) stop("A single umi.count.file should be provided.")
  if (!grepl("\\.txt$", umi.count.file)) stop("umi.count.file should be in .txt format.")
  if (length(annotation.file) != 1) stop("A single annotation.file should be provided.")
  if (!grepl("\\.rds$", annotation.file)) stop("annotation.file should be in .rds format.")
  if (!is.null(blacklist.file) && !grepl("\\.rds$", blacklist.file))
    stop("blacklist.file should be in .rds format.")
  if (is.null(output.prefix)) output.prefix <- sub("\\.txt$", "", basename(umi.count.file))
  if (is.null(feature)) feature <- sub("\\.rds$", "", basename(annotation.file))
  if (!dir.exists(count.tables.output.dir)) dir.create(count.tables.output.dir, recursive = TRUE)

  # ---- Output path ----
  count_table <- file.path(count.tables.output.dir,
                           paste0(output.prefix, "_", feature, "_counts.txt"))

  # ---- Build command ----
  script <- system.file("Rscript", "count_proseq_reads.R", package = "filipR")
  cmd <- paste("Rscript", shQuote(script),
               shQuote(umi.count.file),
               shQuote(annotation.file),
               shQuote(count_table))
  if (!is.null(blacklist.file)) cmd <- paste(cmd, shQuote(blacklist.file))

  # ---- Return data.table ----
  data.table::data.table(
    umi_count_file = umi.count.file,
    count_table = count_table,
    cmd = cmd,
    path = count_table
  )
}
