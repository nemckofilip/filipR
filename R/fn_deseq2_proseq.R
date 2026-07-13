#' DESeq2 Differential Expression for PRO-seq
#'
#' @description
#' Builds a shell command to run DESeq2 differential expression on PRO-seq count
#' tables. No spike-in: size factors come from library size (`libSize`, using the
#' `umi_counts` column of the stats files) or DESeq2 defaults. Outputs a DESeq2
#' object, a fold-change table, and MA plots. Wraps `inst/Rscript/deseq2_proseq_analysis.R`.
#' Submit with [fn_submit()].
#'
#' @param umi.count.tables Character vector of count tables (`.txt`) from [fn_count_proseq_reads()].
#' @param sample.names Character vector of sample names (one per count table).
#' @param conditions Character vector of conditions (one per count table).
#' @param ref.genome.stat.files Character vector of `_UMI_stats.txt` files (one per count table).
#' @param ctl.conditions Character vector of control conditions. Default: all unique conditions.
#' @param normalization "libSize" (default) or "default".
#' @param output.prefix Prefix for output files.
#' @param feature Feature label (e.g. "genes").
#' @param dds.output.dir Directory for the DESeq2 object. Default 'db/dds/'.
#' @param FC.tables.output.dir Directory for the fold-change table. Default 'db/FC_tables/'.
#' @param MAplots.output.dir Directory for the MA plots. Default 'pdf/MAplots/'.
#' @return data.table with columns: dds, FC_table, MA_plot, cmd, path.
#' @import data.table
#' @export
fn_deseq2_proseq <- function(umi.count.tables,
                             sample.names,
                             conditions,
                             ref.genome.stat.files,
                             ctl.conditions = unique(conditions),
                             normalization = "libSize",
                             output.prefix,
                             feature,
                             dds.output.dir = "db/dds/",
                             FC.tables.output.dir = "db/FC_tables/",
                             MAplots.output.dir = "pdf/MAplots/") {

  # ---- Input validation ----
  if (any(duplicated(umi.count.tables))) stop("Some umi.count.tables are duplicated.")
  if (!all(grepl("\\.txt$", umi.count.tables))) stop("umi.count.tables should all be .txt.")
  if (uniqueN(lengths(list(umi.count.tables, sample.names, conditions, ref.genome.stat.files))) != 1)
    stop("umi.count.tables, sample.names, conditions and ref.genome.stat.files should have the same length.")
  if (any(!ctl.conditions %in% conditions)) stop("All ctl.conditions should exist in conditions.")
  if (length(normalization) != 1 || !normalization %in% c("default", "libSize"))
    stop("normalization should be one of 'default' or 'libSize'.")
  if (length(feature) != 1) stop("feature should be unique.")
  for (d in c(dds.output.dir, FC.tables.output.dir, MAplots.output.dir))
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)

  # ---- Output paths ----
  base <- paste0(output.prefix, "_", feature, "_", normalization, "_norm")
  dds_file  <- file.path(dds.output.dir, paste0(base, "_DESeq2.dds"))
  FC_table  <- file.path(FC.tables.output.dir, paste0(base, "_DESeq2_FC.txt"))
  MA_plots  <- file.path(MAplots.output.dir, paste0(base, "_MAplots.pdf"))

  # ---- Build command ----
  script <- system.file("Rscript", "deseq2_proseq_analysis.R", package = "filipR")
  cmd <- paste(
    "Rscript", shQuote(script),
    shQuote(paste(umi.count.tables, collapse = ",")),
    shQuote(paste(ref.genome.stat.files, collapse = ",")),
    shQuote(paste(sample.names, collapse = ",")),
    shQuote(paste(conditions, collapse = ",")),
    shQuote(paste(ctl.conditions, collapse = ",")),
    shQuote(dds.output.dir),
    shQuote(FC.tables.output.dir),
    shQuote(MAplots.output.dir),
    shQuote(output.prefix),
    shQuote(feature),
    shQuote(normalization)
  )

  # ---- Return data.table ----
  data.table::data.table(
    dds = dds_file,
    FC_table = FC_table,
    MA_plot = MA_plots,
    cmd = cmd,
    path = dds_file
  )
}
