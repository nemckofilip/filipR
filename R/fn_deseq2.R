#' Differential Expression Analysis with DESeq2
#'
#' Generates a shell command to run DESeq2 via [fn_submit()]. Compares all
#' non-control conditions against `ctl.condition`. Optionally uses spike-in
#' counts (e.g. ERCC) to set size factors instead of DESeq2's default
#' library-size estimation.
#'
#' @param count.files Character vector of paths to count `.txt` files (output
#'   of [fn_featureCounts()]). One file per sample, columns: gene_id,
#'   gene_name, count.
#' @param sample.names Character vector of sample names, same length as
#'   `count.files`.
#' @param conditions Character vector of condition labels, same length as
#'   `count.files`.
#' @param ctl.conditions Character vector of control condition names. All
#'   pairwise comparisons of non-control vs each control are computed.
#' @param norm.counts Optional numeric vector of spike-in totals (e.g. ERCC
#'   read counts) for manual size factor normalization. Same length as
#'   `count.files`. Default `NULL` uses DESeq2's default estimation.
#' @param padj.cutoff Adjusted p-value cutoff for the `significant` flag
#'   column in the output table. Default `0.05`.
#' @param log2FC.cutoff log2FoldChange cutoff for the `significant` flag.
#'   Default `log2(1.5)`.
#' @param output.prefix Prefix for all output file names.
#' @param dds.output.folder Folder for the saved DESeqDataSet RDS. Default
#'   `"db/dds/"`.
#' @param FC.tables.output.folder Folder for the fold-change table TSV.
#'   Default `"db/FC_tables/"`.
#' @param MAplots.output.folder Folder for the MA plot PDF. Default
#'   `"pdf/MAplots/"`.
#'
#' @return A one-row `data.table` with columns `dds`, `FC.table`, `MA.plot`,
#'   `cmd`, and `path` (= `dds`, used by [fn_submit()] for skip-if-exists).
#' @export
fn_deseq2 <- function(count.files,
                      sample.names,
                      conditions,
                      ctl.conditions,
                      output.prefix,
                      norm.counts            = NULL,
                      padj.cutoff            = 0.05,
                      log2FC.cutoff          = log2(1.5),
                      dds.output.folder      = "db/dds/",
                      FC.tables.output.folder = "db/FC_tables/",
                      MAplots.output.folder  = "pdf/MAplots/") {

  if (length(count.files) != length(sample.names) ||
      length(count.files) != length(conditions))
    stop("count.files, sample.names, and conditions must all have the same length.")
  ctl.conditions <- unique(ctl.conditions)
  if (!all(ctl.conditions %in% conditions))
    stop("All ctl.conditions must be present in conditions.")
  if (!length(setdiff(unique(conditions), ctl.conditions)))
    stop("No non-control conditions found to compare against ctl.conditions.")
  if (!is.null(norm.counts) &&
      (length(norm.counts) != length(count.files) || !is.numeric(norm.counts)))
    stop("norm.counts must be a numeric vector of the same length as count.files.")
  if (any(grepl(",", c(count.files, sample.names, conditions, ctl.conditions))))
    stop("count.files, sample.names, conditions, and ctl.conditions must not contain commas.")

  script   <- system.file("Rscript", "deseq2_analysis.R", package = "filipR")
  dds_file <- file.path(dds.output.folder,       paste0(output.prefix, "_DESeq2.dds"))
  FC_file  <- file.path(FC.tables.output.folder, paste0(output.prefix, "_DESeq2_FC.txt"))
  MA_file  <- file.path(MAplots.output.folder,   paste0(output.prefix, "_MAplots.pdf"))
  norm_arg <- if (!is.null(norm.counts)) paste0(norm.counts, collapse = ",") else "NULL"

  cmd <- paste(
    "Rscript", shQuote(script),
    shQuote(paste0(count.files,  collapse = ",")),
    shQuote(paste0(sample.names, collapse = ",")),
    shQuote(paste0(conditions,    collapse = ",")),
    shQuote(paste0(ctl.conditions, collapse = ",")),
    padj.cutoff,
    log2FC.cutoff,
    shQuote(dds.output.folder),
    shQuote(FC.tables.output.folder),
    shQuote(MAplots.output.folder),
    shQuote(output.prefix),
    shQuote(norm_arg)
  )

  data.table::data.table(
    dds      = dds_file,
    FC.table = FC_file,
    MA.plot  = MA_file,
    cmd      = cmd,
    path     = FC_file
  )
}
