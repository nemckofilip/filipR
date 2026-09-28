#' Enrichment Between Two Fractions of a Pooled Screen
#'
#' @description Tests one fraction of a pooled screen against another with
#'   DESeq2 (e.g. supernatant vs lysate) and returns a fold-change table that
#'   carries the raw counts. Unlike [fn_deseq2()], which is built for RNA-seq,
#'   this runs **in-process** and is not submitted to the cluster - a pool of a
#'   few thousand constructs fits in under a second.
#'
#'   Size factors are DESeq2's default median-of-ratios over all constructs, so
#'   the null is "same share of the library". Spike-in normalisation would
#'   instead put the null at 100% of molecules transferred, which nothing
#'   reaches; absolute efficiency is a different question and not what this
#'   function answers.
#'
#' @param count.files Character vector of per-sample count files. Either the
#'   `construct` / `class` / `molecules` format written by `Count_constructs.R`
#'   (rows with `class != "construct"`, i.e. spike-ins, are dropped), or any
#'   table with a `construct` column and a `molecules` or `count` column. The
#'   counts are expected to be UMI-deduplicated molecules, but nothing in the
#'   function depends on that.
#' @param sample.names Character vector of sample names, same length as
#'   `count.files`.
#' @param conditions Character vector of condition labels, same length as
#'   `count.files`. Exactly two distinct values are allowed.
#' @param ctl.condition The condition to test against, i.e. the assay input
#'   (e.g. `"lysate"`). The reported log2FC is other-vs-this.
#' @param min.input Minimum counts required in **every** `ctl.condition`
#'   replicate. Applied to the input only, so the cut is independent of the
#'   contrast; cutting on the test fraction would select for the effect being
#'   measured. Default `20`, below which a log2 ratio has a standard error of
#'   1.7-2.8 and a construct can read as 60x enriched off a handful of counts.
#' @param padj.cutoff Adjusted p-value cutoff for the `hit` column.
#'   Default `0.001`.
#' @param log2FC.cutoff log2 fold-change cutoff for the `hit` column.
#'   Default `log2(2)`.
#' @param output.prefix Optional prefix for written output. When `NULL`
#'   (default) nothing is written and the table is only returned.
#' @param FC.tables.output.folder Folder for the fold-change table TSV. Only
#'   used when `output.prefix` is given. Default `"db/FC_tables/"`.
#' @param annotation Optional construct annotation, either a `data.table` or a
#'   path to an `.rds`/`.csv`, with a `construct` column and one or more
#'   columns to carry into the output (e.g. `category`, marking designs,
#'   negative and positive controls). Joined onto the result so a hit can be
#'   read against what kind of construct it is. Default `NULL`.
#' @param dds.output.folder Optional folder for the saved DESeqDataSet RDS.
#'   Default `NULL`, which does not save it.
#'
#' @return A `data.table`, one row per tested construct, sorted by decreasing
#'   `log2FC`, with columns: `construct`, any columns carried over from
#'   `annotation`, per-sample counts as `<sample>_raw` and `<sample>_cpm`,
#'   grouped by condition,
#'   `count.<condition>` (raw counts summed over that condition's replicates)
#'   and `total.<condition>` (that condition's library size), `log2FC`,
#'   `padj`, and `hit` (`enriched` / `depleted` / `unaffected`). Spike-ins are
#'   never included - they are dropped on load so they cannot reach the size
#'   factors or the tests. `baseMean`, `lfcSE` and the raw `pvalue` are
#'   deliberately omitted: the counts say what `baseMean` stood in for, and
#'   `padj` is the one that should be read.
#' @import data.table
#' @export
#' @examples
#' \dontrun{
#' fc <- fn_deseq2_enrichment(
#'   count.files   = meta$counts,
#'   sample.names  = meta$id,
#'   conditions    = meta$fraction,
#'   ctl.condition = "lysate",
#'   min.input     = 20,
#'   output.prefix = "MSAR_Peter",
#'   FC.tables.output.folder = "db/FC_tables/export_screen/")
#' }
fn_deseq2_enrichment <- function(count.files,
                                 sample.names,
                                 conditions,
                                 ctl.condition,
                                 min.input               = 20,
                                 padj.cutoff             = 0.001,
                                 log2FC.cutoff           = log2(2),
                                 annotation              = NULL,
                                 output.prefix           = NULL,
                                 FC.tables.output.folder = "db/FC_tables/",
                                 dds.output.folder       = NULL) {

  # ---- Validation ----
  if (length(count.files) != length(sample.names) ||
      length(count.files) != length(conditions))
    stop("count.files, sample.names and conditions must have the same length.")
  if (anyDuplicated(sample.names))
    stop("Duplicated entries in sample.names.")
  if (!all(file.exists(count.files)))
    stop("Missing count files: ",
         paste(count.files[!file.exists(count.files)], collapse = ", "))
  if (length(unique(conditions)) != 2)
    stop("Exactly two distinct conditions are required, got: ",
         paste(unique(conditions), collapse = ", "))
  if (!ctl.condition %in% conditions)
    stop("ctl.condition not found in conditions: ", ctl.condition)
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("DESeq2 is required by fn_deseq2_enrichment().")

  test.condition <- setdiff(unique(conditions), ctl.condition)

  # ---- Load counts ----
  # Spikes are dropped here so they can never enter the size factors
  dat <- data.table::rbindlist(lapply(seq_along(count.files), function(i) {
    d <- data.table::fread(count.files[i])
    cnt.col <- intersect(c("molecules", "count"), names(d))[1]
    if (is.na(cnt.col))
      stop("No 'molecules' or 'count' column in ", count.files[i])
    data.table::data.table(
      sample = sample.names[i], construct = d$construct,
      class = if ("class" %in% names(d)) d$class else "construct",
      count = d[[cnt.col]])
  }))

  cls <- unique(dat[, .(construct, class)])
  DF  <- data.table::dcast(dat, construct ~ factor(sample, sample.names),
                           value.var = "count")
  full <- as.matrix(DF[, !"construct"])
  rownames(full) <- DF$construct
  storage.mode(full) <- "integer"

  # Spikes are held back from the model entirely - they would otherwise steer
  # the size factors and take up multiple-testing budget for a guaranteed null
  spike.rows <- rownames(full) %in% cls[class == "spike", construct]
  mat <- full[!spike.rows, , drop = FALSE]

  # ---- Input floor ----
  # Required in EVERY control replicate, so one deep replicate cannot carry a
  # construct through on its own
  ctl.cols  <- sample.names[conditions == ctl.condition]
  test.cols <- sample.names[conditions == test.condition]
  keep <- rowSums(mat[, ctl.cols, drop = FALSE] >= min.input) == length(ctl.cols)
  message(sprintf("Input floor >= %g in every %s replicate: %d of %d kept",
                  min.input, ctl.condition, sum(keep), nrow(mat)))
  mat <- mat[keep, , drop = FALSE]

  # ---- DESeq2 ----
  coldata <- data.frame(
    condition = factor(conditions, c(ctl.condition, test.condition)),
    row.names = sample.names)
  dds <- DESeq2::DESeqDataSetFromMatrix(mat, coldata, ~ condition)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  if (!is.null(dds.output.folder) && !is.null(output.prefix)) {
    if (!dir.exists(dds.output.folder))
      dir.create(dds.output.folder, recursive = TRUE)
    saveRDS(dds, file.path(dds.output.folder,
                           paste0(output.prefix, "_DESeq2.dds")))
  }

  # independentFiltering off so every construct that passed the floor is
  # reported rather than silently returning padj = NA
  res <- DESeq2::results(dds,
                         contrast = c("condition", test.condition, ctl.condition),
                         independentFiltering = FALSE, cooksCutoff = FALSE)

  # ---- Fold-change table ----
  # Raw counts travel with the result: a log2FC of 6 built on 2 counts and one
  # built on 2,000 are not the same claim, and baseMean cannot tell them apart
  FC <- data.table::data.table(construct = rownames(res))

  # Carry the annotation so a hit can be read against what it is - a design, a
  # negative control or a positive control
  if (!is.null(annotation)) {
    ann <- if (is.character(annotation)) {
      if (grepl("\\.rds$", annotation, ignore.case = TRUE))
        as.data.table(readRDS(annotation)) else data.table::fread(annotation)
    } else as.data.table(annotation)
    if (!"construct" %in% names(ann))
      stop("annotation must have a 'construct' column.")
    FC <- ann[FC, on = "construct"]
    data.table::setcolorder(FC, "construct")
  }

  # Per-sample counts raw and as CPM, grouped by condition rather than
  # interleaved. Raw first, because it is the one that says what a fold change
  # actually rests on. Spike normalisation is deliberately not offered here -
  # it answers a different question (absolute efficiency) and needs the
  # per-fraction spike masses to be meaningful.
  ord <- c(ctl.cols, test.cols)
  for (s in ord) FC[, (paste0(s, "_raw")) := mat[, s]]

  cpm <- sweep(mat, 2, colSums(mat), "/") * 1e6
  for (s in ord) FC[, (paste0(s, "_cpm")) := round(cpm[, s], 2)]



  FC[, (paste0("count.", ctl.condition))  := rowSums(mat[, ctl.cols,  drop = FALSE])]
  FC[, (paste0("count.", test.condition)) := rowSums(mat[, test.cols, drop = FALSE])]
  FC[, (paste0("total.", ctl.condition))  := sum(mat[, ctl.cols])]
  FC[, (paste0("total.", test.condition)) := sum(mat[, test.cols])]

  FC[, log2FC := res$log2FoldChange]
  FC[, padj   := res$padj]
  FC[, hit := data.table::fcase(
    padj < padj.cutoff & log2FC >   log2FC.cutoff,  "enriched",
    padj < padj.cutoff & log2FC < (-log2FC.cutoff), "depleted",
    default = "unaffected")]
  data.table::setorder(FC, -log2FC)

  if (!is.null(output.prefix)) {
    if (!dir.exists(FC.tables.output.folder))
      dir.create(FC.tables.output.folder, recursive = TRUE)
    data.table::fwrite(
      FC, file.path(FC.tables.output.folder,
                    paste0(output.prefix, "_enrichment.txt")),
      sep = "\t", na = "NA")
  }

  FC[]
}
