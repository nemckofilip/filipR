#' Map Reads with HISAT-3N (Conversion-Aware Splice-Aware Alignment)
#'
#' @description
#' Maps reads to a reference genome using HISAT-3N, a conversion-aware
#' extension of HISAT2. Designed for RNA editing experiments where the edited
#' base would cause excessive mismatches in a standard aligner (e.g. A-to-I by
#' TadA variants, C-to-U by APOBEC deaminases). Edited bases in the read are
#' treated as matches rather than mismatches, recovering reads that would
#' otherwise be unmapped. Produces a sorted, indexed BAM file. HISAT-3N adds
#' the \code{Yf:Z:} tag to each alignment recording the number of detected
#' conversions.
#'
#' @param fq1 Character. Path to R1 FASTQ file (gzip supported).
#' @param fq2 Character or \code{NULL}. Path to R2 FASTQ for paired-end;
#'   \code{NULL} for single-end.
#' @param index Character. Path prefix of the HISAT-3N index built with the
#'   same \code{base.change}.
#' @param base.name Character. Sample base name used for output file naming.
#' @param base.change Character. Conversion passed to \code{--base-change},
#'   e.g. \code{"A,G"} for TadA/A-to-I samples or \code{"C,T"} for
#'   APOBEC/C-to-U samples.
#' @param output.dir Character. Directory for output BAM and unmapped FASTQ.
#'   Default: \code{"db/bam/"}.
#' @param save.unmapped Logical. If \code{TRUE}, writes unmapped reads as
#'   gzipped FASTQ. Default: \code{TRUE}.
#' @param cores Integer. Total threads for HISAT-3N and samtools sort.
#'   Default: \code{8}.
#' @param mapq Integer or \code{NULL}. Minimum MAPQ for post-alignment
#'   samtools view filter. Default: \code{NULL}.
#' @param multimapper.mode Character. How to handle multi-mappers:
#'   \code{"best"} (one alignment per read, \code{-k 1}),
#'   \code{"all"} (up to HISAT-3N default of 5 alignments per read),
#'   \code{"unique"} (only uniquely mapped reads; forces \code{mapq = 60}).
#'   Default: \code{"best"}.
#' @param rna.strandness Character or \code{NULL}. Passed to
#'   \code{--rna-strandness}: \code{"RF"}, \code{"FR"}, or \code{NULL} for
#'   unstranded. Default: \code{NULL}.
#' @param alignment.stats.output.dir Character. Directory for alignment summary
#'   log. Default: \code{"db/alignment_stats/"}.
#'
#' @return A \code{data.table} with columns: \code{fq1_in}, \code{fq2_in},
#'   \code{bam}, \code{stats}, \code{fq1_unmapped}, \code{fq2_unmapped},
#'   \code{cmd}, \code{path}.
#' @export
fn_hisat3n_map <- function(fq1,
                            fq2                        = NULL,
                            index,
                            base.name,
                            base.change,
                            output.dir                 = "db/bam/",
                            save.unmapped              = TRUE,
                            cores                      = 8,
                            mapq                       = NULL,
                            multimapper.mode           = "best",
                            rna.strandness             = NULL,
                            alignment.stats.output.dir = "db/alignment_stats/") {

  # ---- Input validation ----
  if (length(fq1) != 1) stop("fn_hisat3n_map processes one sample at a time.")
  if (!is.null(fq2) && length(fq2) != 1) stop("fn_hisat3n_map processes one sample at a time.")
  if (!grepl("^[A-Z],[A-Z]$", base.change))
    stop("base.change must be two comma-separated bases, e.g. 'A,G'.")
  if (!multimapper.mode %in% c("best", "all", "unique"))
    stop("multimapper.mode must be 'best', 'all', or 'unique'.")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
  if (!dir.exists(alignment.stats.output.dir))
    dir.create(alignment.stats.output.dir, recursive = TRUE)

  is_paired <- !is.null(fq2)
  bc_tag    <- gsub(",", "", base.change)        # "AG" or "CT" for file naming

  # ---- Output paths ----
  bam_out   <- file.path(output.dir,
                         paste0(base.name, "_hisat3n_", bc_tag, ".bam"))
  stats_out <- file.path(alignment.stats.output.dir,
                         paste0(base.name, "_hisat3n_", bc_tag, ".stats"))

  # ---- Unmapped read paths ----
  fq1_unmapped <- fq2_unmapped <- NA_character_
  unmapped_flag <- ""
  if (save.unmapped) {
    if (is_paired) {
      fq1_unmapped    <- file.path(output.dir,
        paste0(base.name, "_hisat3n_", bc_tag, "_unmapped_1.fq.gz"))
      fq2_unmapped    <- file.path(output.dir,
        paste0(base.name, "_hisat3n_", bc_tag, "_unmapped_2.fq.gz"))
      unmapped_pattern <- file.path(output.dir,
        paste0(base.name, "_hisat3n_", bc_tag, "_unmapped_%.fq.gz"))
      unmapped_flag   <- paste("--un-conc-gz", shQuote(unmapped_pattern))
    } else {
      fq1_unmapped  <- file.path(output.dir,
        paste0(base.name, "_hisat3n_", bc_tag, "_unmapped.fq.gz"))
      unmapped_flag <- paste("--un-gz", shQuote(fq1_unmapped))
    }
  }

  # ---- Thread balancing ----
  sort_cores  <- 2L
  hisat_cores <- max(1L, cores - sort_cores)

  # ---- Multi-mapper handling ----
  multimapper_flag <- switch(multimapper.mode,
    best   = "-k 1",
    all    = "",
    unique = { mapq <- 60L; "-k 1" }
  )

  # ---- Read input flags ----
  read_flags <- if (is_paired) {
    paste("-1", shQuote(fq1), "-2", shQuote(fq2))
  } else {
    paste("-U", shQuote(fq1))
  }
  pe_flags <- if (is_paired) "--no-mixed --no-discordant" else ""

  # ---- Strandedness ----
  strand_flag <- if (!is.null(rna.strandness))
    paste("--rna-strandness", rna.strandness) else ""

  # ---- Build HISAT-3N command ----
  cmd_align <- paste(
    "hisat-3n",
    "--index",        shQuote(index),
    "--base-change",  base.change,
    read_flags,
    "-p",             hisat_cores,
    pe_flags,
    multimapper_flag,
    strand_flag,
    unmapped_flag,
    "--summary-file", shQuote(stats_out)
  )

  # Pipe into samtools sort
  sort_cmd <- paste(
    "samtools sort",
    "-@", sort_cores,
    "-m 2G",
    "-T", shQuote(paste0(bam_out, "_tmp")),
    "-o", shQuote(bam_out),
    "-"
  )

  if (!is.null(mapq)) {
    cmd <- paste(cmd_align, "| samtools view -b -q", mapq, "|", sort_cmd)
  } else {
    cmd <- paste(cmd_align, "|", sort_cmd)
  }

  cmd <- paste(cmd, "&&", "samtools index", shQuote(bam_out))

  # ---- Return ----
  data.table::data.table(
    fq1_in       = fq1,
    fq2_in       = if (is_paired) fq2 else NA_character_,
    bam          = bam_out,
    stats        = stats_out,
    fq1_unmapped = if (save.unmapped) fq1_unmapped else NA_character_,
    fq2_unmapped = if (save.unmapped && is_paired) fq2_unmapped else NA_character_,
    cmd          = cmd,
    path         = bam_out
  )
}
