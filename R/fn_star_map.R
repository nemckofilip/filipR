#' Map Reads with STAR (Single-End or Paired-End)
#'
#' @description
#' Maps reads to a reference genome using STAR.
#' Converts mapped reads to a sorted BAM file.
#' Preserves UMIs in read headers.
#'
#' @param fq1 Character string. Path to R1 FASTQ file.
#' @param fq2 Character string. Path to R2 FASTQ file. Default NULL (single-end).
#' @param index Path to STAR genome index directory.
#' @param base.name Base name for output files.
#' @param output.dir Directory for BAM and unmapped FASTQ. Default: 'db/bam/'.
#' @param save.unmapped Logical. If TRUE, writes unmapped reads. Default TRUE.
#' @param cores Total number of threads to use. Default 8.
#' @param mapq Optional minimum MAPQ to filter alignments. Default NULL.
#' @param outFilterMultimapNmax Maximum number of loci a read can map to. Default 10.
#' @param outFilterMismatchNmax Maximum number of mismatches per read. Default 6.
#' @param multimapper.mode How to handle multi-mappers: "best" (default, random assignment, output 1), "all" (output all alignments), or "unique" (only unique mappers).
#' @param alignment.stats.output.dir Optional directory for STAR log files. Default: 'db/alignment_stats/'.
#'
#' @return A data.table with columns: fq1_in, fq2_in, bam, stats, fq1_unmapped, fq2_unmapped, cmd, job.name.
#' @export
fn_star_map <- function(fq1, 
                        fq2 = NULL,
                        index, 
                        base.name,
                        output.dir = "db/bam/", 
                        save.unmapped = TRUE,
                        cores = 8,
                        mapq = NULL,
                        outFilterMultimapNmax = 10,
                        outFilterMismatchNmax = 6,
                        multimapper.mode = "best",
                        alignment.stats.output.dir = "db/alignment_stats/") {

  # ---- Input validation ----
  if (length(fq1) > 1) stop("fn_star_map processes one sample at a time.")
  if (!is.null(fq2) && length(fq2) > 1) stop("fn_star_map processes one sample at a time.")
  if (!multimapper.mode %in% c("best", "all", "unique")) stop("multimapper.mode must be 'best', 'all', or 'unique'")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
  if (!dir.exists(alignment.stats.output.dir)) dir.create(alignment.stats.output.dir, recursive = TRUE)

  is_paired <- !is.null(fq2)

  # ---- Thread & memory balancing ----
  sort_cores <- 2
  star_cores <- max(1, cores - sort_cores)
  sort_mem <- "2G"

  # ---- Output paths ----
  star_prefix <- file.path(output.dir, paste0(base.name, "_TMP_"))
  star_unsorted_bam <- paste0(star_prefix, "Aligned.out.bam")
  bam_out <- file.path(output.dir, paste0(base.name, "_genome.bam"))
  stats_out <- file.path(alignment.stats.output.dir, paste0(base.name, "_genome.stats"))

  # ---- Unmapped read paths ----
  if (is_paired) {
    unmapped_raw_1 <- paste0(star_prefix, "Unmapped.out.mate1")
    unmapped_raw_2 <- paste0(star_prefix, "Unmapped.out.mate2")
    fq1_unmapped <- file.path(output.dir, paste0(base.name, "_genome_unmapped_1.fq.gz"))
    fq2_unmapped <- file.path(output.dir, paste0(base.name, "_genome_unmapped_2.fq.gz"))
  } else {
    unmapped_raw <- paste0(star_prefix, "Unmapped.out.mate1")
    fq1_unmapped <- file.path(output.dir, paste0(base.name, "_genome_unmapped.fq.gz"))
    fq2_unmapped <- NA_character_
  }

  # ---- Handle input compression ----
  read_cmd_flag <- if (grepl("\\.gz$", fq1)) "--readFilesCommand zcat" else ""

  # ---- Handle multi-mapper mode ----
  if (multimapper.mode == "best") {
    # Random assignment, output only 1 alignment per read
    multimapper_flags <- paste(
      "--outMultimapperOrder Random",
      "--outSAMmultNmax 1"
    )
  } else if (multimapper.mode == "unique") {
    # Only unique mappers
    outFilterMultimapNmax <- 1
    multimapper_flags <- ""
  } else {
    # "all" mode: output all alignments (default STAR behavior)
    multimapper_flags <- ""
  }

  # ---- Build STAR command ----
  cmd_align <- paste(
    "STAR",
    "--runThreadN", star_cores,
    "--genomeDir", shQuote(index),
    "--readFilesIn", shQuote(fq1), if (is_paired) shQuote(fq2) else "",
    read_cmd_flag,
    "--outFileNamePrefix", shQuote(star_prefix),
    "--outSAMtype BAM Unsorted",
    "--outSAMattributes NH HI AS nM NM MD",
    "--outFilterMultimapNmax", outFilterMultimapNmax,
    "--outFilterMismatchNmax", outFilterMismatchNmax,
    "--outSAMstrandField intronMotif",
    "--alignEndsType Local",
    multimapper_flags
  )

  if (save.unmapped) cmd_align <- paste(cmd_align, "--outReadsUnmapped Fastx")

  # ---- Build pipeline ----
  cmd <- paste(cmd_align, "&&")

  # Sort BAM
  sort_cmd <- paste("samtools sort -@", sort_cores, "-m", sort_mem, "-o", shQuote(bam_out))
  if (!is.null(mapq)) {
    cmd <- paste(cmd, "samtools view -b -q", mapq, shQuote(star_unsorted_bam), "|", sort_cmd, "-")
  } else {
    cmd <- paste(cmd, sort_cmd, shQuote(star_unsorted_bam))
  }

  # Cleanup unsorted BAM
  cmd <- paste(cmd, "&&", "rm", shQuote(star_unsorted_bam))

  # Move log file
  log_default <- paste0(star_prefix, "Log.final.out")
  cmd <- paste(cmd, "&&", "mv", shQuote(log_default), shQuote(stats_out))

  # Cleanup STAR temp files
  cmd <- paste(cmd, "&&", "rm", 
               paste0(star_prefix, "Log.out"), 
               paste0(star_prefix, "Log.progress.out"), 
               paste0(star_prefix, "SJ.out.tab"))

  # Gzip unmapped reads
  if (save.unmapped) {
    if (is_paired) {
      cmd <- paste(cmd, "&& gzip -c", shQuote(unmapped_raw_1), ">", shQuote(fq1_unmapped))
      cmd <- paste(cmd, "&& gzip -c", shQuote(unmapped_raw_2), ">", shQuote(fq2_unmapped))
      cmd <- paste(cmd, "&& rm", shQuote(unmapped_raw_1), shQuote(unmapped_raw_2))
    } else {
      cmd <- paste(cmd, "&& gzip -c", shQuote(unmapped_raw), ">", shQuote(fq1_unmapped))
      cmd <- paste(cmd, "&& rm", shQuote(unmapped_raw))
    }
  }

  # ---- Return data.table ----
  data.table::data.table(
    fq1_in = fq1,
    fq2_in = if (is_paired) fq2 else NA_character_,
    bam = bam_out,
    stats = stats_out,
    fq1_unmapped = if (save.unmapped) fq1_unmapped else NA_character_,
    fq2_unmapped = if (save.unmapped) fq2_unmapped else NA_character_,
    cmd = cmd,
    job.name = "STAR",
    path = bam_out
  )
}