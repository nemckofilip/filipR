#' Count RNA-Seq Reads Using featureCounts (Subread)
#'
#' @description
#' Counts reads from a BAM file using the featureCounts command-line tool.
#' Outputs a clean counts file (gene_id, gene_name, count) and a statistics file.
#' Supports single-end and paired-end data, stranded libraries.
#'
#' @param bam Character string. Path to the input BAM file. Only one file allowed.
#' @param base.name Base name for output files.
#' @param annotation.file Path to a GTF annotation file.
#' @param paired Logical. TRUE for paired-end, FALSE for single-end. Default FALSE.
#' @param strandSpecific Strand specificity: 0 = unstranded, 1 = forward, 2 = reverse. Default 1.
#' @param allowMultiOverlap Logical. Assign reads to all overlapping features? Default FALSE.
#' @param countMultiMapping Logical. Count multi-mapping reads? Default TRUE.
#' @param output.dir Directory for counts and stats files. Default "db/counts/".
#' @param cores Number of threads. Default 4.
#'
#' @return A data.table with columns: bam_in, counts, stats, cmd, job.name, path.
#' @export
fn_featureCounts <- function(bam,
                              base.name,
                              annotation.file,
                              paired = FALSE,
                              strandSpecific = 1,
                              allowMultiOverlap = FALSE,
                              countMultiMapping = TRUE,
                              output.dir = "db/counts/",
                              cores = 4) {
  
  # ---- Input validation ----
  if (length(bam) != 1) stop("fn_featureCounts processes one BAM file at a time.")
  if (!strandSpecific %in% c(0, 1, 2)) stop("strandSpecific must be 0, 1, or 2.")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
  
  # ---- Output paths ----
  counts_out <- file.path(output.dir, paste0(base.name, "_counts.txt"))
  counts_full <- paste0(counts_out, ".full")
  stats_out <- paste0(counts_full, ".summary")
  
  # ---- Build command ----
  cmd <- paste0(
    "featureCounts",
    if (paired) " -p --countReadPairs" else "",
    " -s ", strandSpecific,
    " -T ", cores,
    " -t gene",
    " -g gene_id",
    if (countMultiMapping) " -M" else "",
    " --extraAttributes gene_name",
    if (allowMultiOverlap) " -O" else "",
    " -a ", annotation.file,
    " -o ", counts_full,
    " ", bam,
    " && echo -e 'gene_id\\tgene_name\\tcount' > ", counts_out,
    " && tail -n +3 ", counts_full, " | cut -f1,7,8 | sed 's/\\.[0-9]*\\t/\\t/' >> ", counts_out,
    " && rm ", counts_full
  )
  
  # ---- Return data.table ----
  data.table::data.table(
    bam_in = bam,
    counts = counts_out,
    stats = stats_out,
    cmd = cmd,
    job.name = "featureCounts",
    path = counts_out
  )
}
