#' Extract UMI from Read 2 and Move to Header
#'
#' @description
#' Uses a custom Perl script to extract the first N bases of Read 2 (UMI) 
#' and append it to the read ID of both Read 1 and Read 2.
#' Read 2 is then trimmed by N bases.
#'
#' @param fq1 Character vector of R1 FASTQ files.
#' @param fq2 Character vector of R2 FASTQ files.
#' @param output.dir Directory for output FASTQ files. Default "db/fq_umi/".
#' @param umi.len Length of the UMI to extract. Default 8.
#' @param script.path Path to the 'extract_umi_read2.pl' script. 
#' Defaults to the script bundled with this package.
#' @param cores Number of threads for compression (pigz). Default 4.
#'
#' @return A data.table with the commands and new file paths.
#' @export
fn_extract_umi_read2 <- function(fq1, 
                           fq2, 
                           output.dir = "db/fq_trimmed_umi/", 
                           umi.len = 8,
                           script.path = system.file("perl", "extract_umi_read2.pl", package = "filipR"),
                           cores = 4) {
  
  # Checks
  if(script.path == "") stop("Perl script 'extract_umi_read2.pl' not found. Please ensure it is in inst/perl/ or provide a path.")
  if(!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
  if(length(fq1) != length(fq2)) stop("fq1 and fq2 lengths differ")

  # Generate Paths
  fq1_out <- file.path(output.dir, paste0(sub("\\.(fq|fastq)(\\.gz)?$", "", basename(fq1)), "_umi.fq.gz"))
  fq2_out <- file.path(output.dir, paste0(sub("\\.(fq|fastq)(\\.gz)?$", "", basename(fq2)), "_umi.fq.gz"))

  # Generate Commands
  cmd <- paste("perl", script.path, 
               "--r1", fq1, 
               "--r2", fq2, 
               "--out1", fq1_out, 
               "--out2", fq2_out,
               "--umi-len", umi.len,
               "--threads", cores)

  dt <- data.table(
    file.type = rep(c("fq1.umi", "fq2.umi"), each = length(fq1)),
    path = c(fq1_out, fq2_out),
    cmd = rep(cmd, 2), 
    job.name = "umi_extract"
  )
  
  return(dt)
}

