#' Extract UMI from Read 2 and Move to Header
#'
#' @description
#' Uses a custom Perl script to extract the first N bases of Read 2 (UMI) 
#' and append it to the read ID of both Read 1 and Read 2.
#' Read 2 is then trimmed by N bases.
#'
#' @param fq1 Character string. Path to R1 FASTQ file.
#' @param fq2 Character string. Path to R2 FASTQ file.
#' @param base.name Base name for output files.
#' @param output.dir Directory for output FASTQ files. Default "db/fq_trimmed_umi/".
#' @param umi.len Length of the UMI to extract. Default 8.
#' @param script.path Path to the 'extract_umi_read2.pl' script. 
#' Defaults to the script bundled with this package.
#' @param cores Number of threads for compression (pigz). Default 4.
#'
#' @return A data.table with columns: fq1_in, fq2_in, fq1_out, fq2_out, cmd.
#' @export
fn_extract_umi_read2 <- function(fq1, 
                                 fq2, 
                                 base.name,
                                 output.dir = "db/fq_trimmed_umi/", 
                                 umi.len = 8,
                                 script.path = system.file("perl", "extract_umi_read2.pl", package = "filipR"),
                                 cores = 4) {
  
  # ---- Input validation ----
  if (length(fq1) != 1) stop("fn_extract_umi_read2 processes one sample at a time.")
  if (length(fq2) != 1) stop("fn_extract_umi_read2 processes one sample at a time.")
  if (script.path == "") stop("Perl script 'extract_umi_read2.pl' not found. Please ensure it is in inst/perl/ or provide a path.")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  # ---- Output paths ----
  fq1_out <- file.path(output.dir, paste0(base.name, "_umi_R1.fq.gz"))
  fq2_out <- file.path(output.dir, paste0(base.name, "_umi_R2.fq.gz"))

  # ---- Build command ----
  cmd <- paste("perl", script.path, 
               "--r1", fq1, 
               "--r2", fq2, 
               "--out1", fq1_out, 
               "--out2", fq2_out,
               "--umi-len", umi.len,
               "--threads", cores)

  # ---- Return data.table ----
  data.table::data.table(
    fq1_in = fq1,
    fq2_in = fq2,
    fq1_out = fq1_out,
    fq2_out = fq2_out,
    cmd = cmd,
    path = fq1_out
  )
}