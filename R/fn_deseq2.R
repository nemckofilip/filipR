#' Differential Expression Analysis with DESeq2
#'
#' Generates a shell command to run DESeq2 via [fn_submit()]. Compares all
#' non-control conditions against `ctl.condition`. Optionally appends spike-in
#' count rows (e.g. ERCC) and uses them as DESeq2 `controlGenes`, so size
#' factors come from DESeq2's own median-of-ratios estimator restricted to the
#' spike-ins instead of the default all-gene estimation. Spike-in rows are
#' dropped from the reported results.
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
#' @param spikein.count.files Optional character vector of paths to per-spike-in
#'   count files (e.g. per-ERCC counts, same 3-column format as `count.files`),
#'   one per sample. When supplied these rows are appended to the count matrix
#'   and used as DESeq2 `controlGenes` for size-factor estimation. Same length
#'   as `count.files`. Default `NULL` uses DESeq2's default all-gene estimation.
#' @param keep.genes.file Optional path to a file with one gene_id per line.
#'   When supplied, gene rows are restricted to these IDs before testing (e.g.
#'   to drop pseudogenes), improving multiple-testing correction. Spike-in rows
#'   are always retained. Default `NULL` keeps all genes.
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
                      spikein.count.files    = NULL,
                      keep.genes.file        = NULL,
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
  if (!is.null(spikein.count.files) &&
      length(spikein.count.files) != length(count.files))
    stop("spikein.count.files must have the same length as count.files.")
  if (any(grepl(",", c(count.files, sample.names, conditions, ctl.conditions,
                       spikein.count.files))))
    stop("count.files, sample.names, conditions, ctl.conditions, and spikein.count.files must not contain commas.")
  if (!is.null(keep.genes.file) && !file.exists(keep.genes.file))
    stop("keep.genes.file does not exist: ", keep.genes.file)

  script   <- system.file("Rscript", "deseq2_analysis.R", package = "filipR")
  dds_file <- file.path(dds.output.folder,       paste0(output.prefix, "_DESeq2.dds"))
  FC_file  <- file.path(FC.tables.output.folder, paste0(output.prefix, "_DESeq2_FC.txt"))
  MA_file  <- file.path(MAplots.output.folder,   paste0(output.prefix, "_MAplots.pdf"))
  spikein_arg <- if (!is.null(spikein.count.files))
    paste0(spikein.count.files, collapse = ",") else "NULL"
  keep_arg <- if (!is.null(keep.genes.file)) keep.genes.file else "NULL"

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
    shQuote(spikein_arg),
    shQuote(keep_arg)
  )

  data.table::data.table(
    dds      = dds_file,
    FC.table = FC_file,
    MA.plot  = MA_file,
    cmd      = cmd,
    path     = FC_file
  )
}
