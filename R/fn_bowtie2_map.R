#' Map Reads with Bowtie2
#'
#' @name fn_bowtie2_map
#' @description
#' Maps reads to a reference index using Bowtie2.
#' Converts mapped reads to a sorted BAM file.
#' Preserves UMIs in read headers.
#' Optionally filters by MAPQ, sets max insert size, and writes unmapped reads.
#'
#' @param fq1 Character vector of R1 FASTQ files.
#' @param fq2 Character vector of R2 FASTQ files. Can be NULL for single-end mapping.
#' @param index Path to Bowtie2 index prefix.
#' @param base.name Base name for output files (BAM, stats, unmapped). Overrides fq filenames.
#' @param output.dir Directory for BAM and unmapped FASTQ. Default: 'db/bowtie2/'.
#' @param suffix Suffix to append to output file names (e.g., '_spikein', '_genome'). Default: ''.
#' @param save.unmapped Logical. If TRUE, writes unmapped reads. Default TRUE.
#' @param global Logical. If TRUE (default), performs end-to-end alignment; if FALSE, local alignment.
#' @param cores Total number of threads to use. Default 8.
#' @param mapq Optional minimum MAPQ to filter alignments. Default NULL (no filtering).
#' @param max.ins Optional maximum insert size for paired-end reads. Default NULL (Bowtie2 default).
#' @param multimapper.mode How to handle multi-mappers: "best" (default, report best/random), "all" (report all up to -k), or "unique" (only unique mappers).
#' @param alignment.stats.output.dir Optional directory for Bowtie2 alignment stats. Defaults to output.dir.
#'
#' @return A data.table with columns: fq1_in, fq2_in, bam, stats, fq1_unmapped, fq2_unmapped, cmd, job.name.
#' @export
fn_bowtie2_map <- function(fq1, 
                           fq2 = NULL, 
                           index, 
                           base.name,
                           output.dir = "db/bowtie2/", 
                           suffix = "",
                           save.unmapped = TRUE,
                           global = TRUE,
                           cores = 8,
                           mapq = NULL,
                           max.ins = NULL,
                           multimapper.mode = "best",
                           alignment.stats.output.dir = output.dir) {

  # Checks
  if(!is.null(fq2) && length(fq1) != length(fq2)) stop("fq1 and fq2 lengths differ")
  if(!multimapper.mode %in% c("best", "all", "unique")) stop("multimapper.mode must be 'best', 'all', or 'unique'")
  if(!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
  if(!dir.exists(alignment.stats.output.dir)) dir.create(alignment.stats.output.dir, recursive = TRUE)

  # Thread Balancing
  sort_cores <- 2
  bt2_cores <- max(1, cores - sort_cores)

  # Output paths based on base.name and suffix
  bam_out <- file.path(output.dir, paste0(base.name, suffix, ".bam"))
  stats_out <- file.path(alignment.stats.output.dir, paste0(base.name, suffix, ".stats"))

  if(!is.null(fq2)) {
    unmapped_base <- file.path(output.dir, paste0(base.name, suffix, "_unmapped_%.fq.gz"))
    fq1_unmapped <- sub("%", "1", unmapped_base)
    fq2_unmapped <- sub("%", "2", unmapped_base)
  } else {
    fq1_unmapped <- file.path(output.dir, paste0(base.name, suffix, "_unmapped.fq.gz"))
    fq2_unmapped <- NA_character_
  }

  # Build Bowtie2 command
  if(is.null(fq2)) {
    cmd_align <- paste(
      "bowtie2",
      "-p", bt2_cores,
      "-x", shQuote(index),
      "-U", shQuote(fq1),
      if(global) "--very-sensitive" else "--very-sensitive-local",
      "--no-unal"
    )
    if(!is.null(mapq)) warning("MAPQ filtering for single-end will be applied after BAM sorting")
    if(save.unmapped) cmd_align <- paste(cmd_align, "--un-gz", shQuote(fq1_unmapped))
  } else {
    cmd_align <- paste(
      "bowtie2",
      "-p", bt2_cores,
      "-x", shQuote(index),
      "-1", shQuote(fq1),
      "-2", shQuote(fq2),
      if(global) "--very-sensitive" else "--very-sensitive-local",
      "--no-unal"
    )
    if(!is.null(max.ins)) cmd_align <- paste(cmd_align, "-X", max.ins)
    if(save.unmapped) cmd_align <- paste(cmd_align, "--un-conc-gz", shQuote(unmapped_base))
  }

  # Handle multi-mapper mode
  if(multimapper.mode == "best") {
    cmd_align <- paste(cmd_align, "-k 1")
  } else if(multimapper.mode == "unique") {
    # Filter for MAPQ > 0 later (unique mappers have MAPQ > 0 in Bowtie2)
    if(is.null(mapq) || mapq < 1) mapq <- 1
  }
  # "all" mode: default Bowtie2 behavior, no changes needed

  # Pipeline: Align -> (Filter) -> Sort
  cmd <- paste(cmd_align, "2>", shQuote(stats_out))
  if(!is.null(mapq)) {
    cmd <- paste(cmd, "| samtools view -b -q", mapq)
  }
  cmd <- paste(cmd, "| samtools sort -@", sort_cores, "-o", shQuote(bam_out), "-")

  # Return metadata
  data.table::data.table(
    fq1_in = fq1,
    fq2_in = if(is.null(fq2)) NA_character_ else fq2,
    bam = bam_out,
    stats = stats_out,
    fq1_unmapped = if(save.unmapped) fq1_unmapped else NA_character_,
    fq2_unmapped = if(save.unmapped) fq2_unmapped else NA_character_,
    cmd = cmd,
    job.name = "bowtie2",
    path = bam_out
  )
}