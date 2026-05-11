#!/usr/bin/env Rscript
# Usage: Rscript count_edits_per_read.R <bam> <sample_name> <output>
#
# Per-FRAGMENT TRANSCRIPT-LEVEL A>G and C>T substitution counts.
# Forward-stranded paired-end library (R1 sense, R2 antisense).
#
# For each mate, parses MD to count MD-level A>G, C>T, T>C, G>A and the
# reference A/C/T/G in the aligned region. The transcript strand of the
# fragment is derived from the BAM flag:
#   transcript on + strand if  (is_R1 XOR is_reverse_mapped)
# For + strand fragments, MD A>G = transcript A>G and ref A = transcript A.
# For - strand fragments, MD T>C = transcript A>G and ref T = transcript A
# (analogously G>A = transcript C>T, ref G = transcript C). Mates of a pair
# are joined by read name and their per-mate contributions summed.
# R1/R2 overlap is NOT deduplicated — overlap positions are double-counted
# in both numerator and denominator, so per-fragment rates remain unbiased
# but per-fragment counts are inflated relative to true molecule content.
#
# Output: aggregated TSV with columns:
#   sample, region, n_AG, n_CT, n_ref_A, n_ref_C, n_reads
# where the four count columns are TRANSCRIPT-level per FRAGMENT and
# n_reads is the number of fragments.

args        <- commandArgs(trailingOnly = TRUE)
bam         <- args[1]
sample_name <- args[2]
output      <- args[3]

suppressPackageStartupMessages({
  library(Rsamtools)
  library(data.table)
})

#==============================================================================
# PARSE BAM IN CHUNKS (avoids loading all reads into memory at once)
#==============================================================================

CHUNK <- 500000L

flag_filt <- scanBamFlag(
  isUnmappedQuery          = FALSE,
  isSecondaryAlignment     = FALSE,
  isSupplementaryAlignment = FALSE
)
param <- ScanBamParam(
  flag = flag_filt,
  what = c("qname", "flag", "rname", "seq", "cigar"),
  tag  = "MD"
)

#==============================================================================
# COUNT SUBSTITUTIONS PER MATE (MD-level, before strand resolution)
#==============================================================================

# Returns c(n_AG, n_CT, n_TC, n_GA, n_A, n_C, n_T, n_G)
count_one <- function(md_i, seq_i, lsc_i, tsc_i) {

  aln_seq <- substr(seq_i, lsc_i + 1L, nchar(seq_i) - tsc_i)
  n_A_aln <- nchar(gsub("[^A]", "", aln_seq))
  n_C_aln <- nchar(gsub("[^C]", "", aln_seq))
  n_T_aln <- nchar(gsub("[^T]", "", aln_seq))
  n_G_aln <- nchar(gsub("[^G]", "", aln_seq))

  toks <- regmatches(md_i, gregexpr("\\d+|[A-Z]|\\^[A-Z]+", md_i, perl = TRUE))[[1]]
  if (!length(toks)) return(c(0L, 0L, 0L, 0L, n_A_aln, n_C_aln, n_T_aln, n_G_aln))

  is_let <- nchar(toks) == 1L & grepl("^[A-Z]$", toks, perl = TRUE)
  steps  <- ifelse(grepl("^\\d+$", toks, perl = TRUE), as.integer(toks),
                   ifelse(is_let, 1L, 0L))

  let_idx <- which(is_let)
  if (!length(let_idx)) return(c(0L, 0L, 0L, 0L, n_A_aln, n_C_aln, n_T_aln, n_G_aln))

  cum_steps  <- c(0L, cumsum(steps))
  read_pos   <- cum_steps[let_idx] + lsc_i

  ref_bases  <- toks[let_idx]
  read_bases <- substr(rep(seq_i, length(read_pos)), read_pos + 1L, read_pos + 1L)

  n_ag <- sum(ref_bases == "A" & read_bases == "G")
  n_ct <- sum(ref_bases == "C" & read_bases == "T")
  n_tc <- sum(ref_bases == "T" & read_bases == "C")
  n_ga <- sum(ref_bases == "G" & read_bases == "A")

  n_ref_A <- n_A_aln - sum(read_bases == "A") + sum(ref_bases == "A")
  n_ref_C <- n_C_aln - sum(read_bases == "C") + sum(ref_bases == "C")
  n_ref_T <- n_T_aln - sum(read_bases == "T") + sum(ref_bases == "T")
  n_ref_G <- n_G_aln - sum(read_bases == "G") + sum(ref_bases == "G")

  c(n_ag, n_ct, n_tc, n_ga, n_ref_A, n_ref_C, n_ref_T, n_ref_G)
}

process_chunk <- function(raw) {
  qname <- raw$qname
  flag  <- raw$flag
  rname <- as.character(raw$rname)
  seq_v <- as.character(raw$seq)
  cigar <- raw$cigar
  md    <- raw$tag$MD

  keep  <- !is.na(md) & nchar(seq_v) > 0L
  if (!any(keep)) return(NULL)

  qname <- qname[keep]; flag  <- flag[keep]
  rname <- rname[keep]; seq_v <- seq_v[keep]
  cigar <- cigar[keep]; md    <- md[keep]

  lsc <- rep(0L, length(md))
  has_lsc <- grepl("^\\d+S", cigar, perl = TRUE)
  lsc[has_lsc] <- as.integer(sub("^(\\d+)S.*", "\\1", cigar[has_lsc], perl = TRUE))

  tsc <- rep(0L, length(md))
  has_tsc <- grepl("\\d+S$", cigar, perl = TRUE)
  tsc[has_tsc] <- as.integer(sub(".*?(\\d+)S$", "\\1", cigar[has_tsc], perl = TRUE))

  res    <- mapply(count_one, md, seq_v, lsc, tsc, SIMPLIFY = TRUE)
  region <- ifelse(startsWith(rname, "ERCC"), "ercc", "human")

  # Flag bits: 0x40 = R1 (first in pair), 0x10 = reverse-mapped
  is_R1       <- bitwAnd(flag, 64L) != 0L
  is_reverse  <- bitwAnd(flag, 16L) != 0L
  plus_strand <- xor(is_R1, is_reverse)

  n_AG_md <- as.integer(res[1L, ]); n_CT_md <- as.integer(res[2L, ])
  n_TC_md <- as.integer(res[3L, ]); n_GA_md <- as.integer(res[4L, ])
  n_A_md  <- as.integer(res[5L, ]); n_C_md  <- as.integer(res[6L, ])
  n_T_md  <- as.integer(res[7L, ]); n_G_md  <- as.integer(res[8L, ])

  data.table(
    qname   = qname,
    region  = region,
    n_AG    = ifelse(plus_strand, n_AG_md, n_TC_md),
    n_CT    = ifelse(plus_strand, n_CT_md, n_GA_md),
    n_ref_A = ifelse(plus_strand, n_A_md,  n_T_md),
    n_ref_C = ifelse(plus_strand, n_C_md,  n_G_md)
  )
}

bf      <- BamFile(bam, yieldSize = CHUNK)
open(bf)
chunks  <- list()
i       <- 0L
repeat {
  raw <- scanBam(bf, param = param)[[1]]
  if (length(raw$rname) == 0L) break
  i <- i + 1L
  chunks[[i]] <- process_chunk(raw)
}
close(bf)

#==============================================================================
# JOIN MATES BY QNAME AND AGGREGATE
#==============================================================================

dt     <- rbindlist(chunks)
chunks <- NULL; gc()

frag <- dt[, .(region  = region[1L],
               n_AG    = sum(n_AG),
               n_CT    = sum(n_CT),
               n_ref_A = sum(n_ref_A),
               n_ref_C = sum(n_ref_C)),
           by = qname]
dt <- NULL; gc()

result <- frag[, .(n_reads = .N), by = .(region, n_AG, n_CT, n_ref_A, n_ref_C)]
result[, sample := sample_name]
setcolorder(result, c("sample", "region", "n_AG", "n_CT", "n_ref_A", "n_ref_C", "n_reads"))
fwrite(result, output, sep = "\t")
message("Written: ", output)
