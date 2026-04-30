#!/usr/bin/env Rscript
# Usage: Rscript count_edits_per_read.R <bam> <sample_name> <output>
#
# For each mapped primary read: parses MD tag to count A>G and C>T substitutions
# and the total reference A and C bases in the aligned region.
#
# Output: aggregated TSV with columns:
#   sample, region, n_AG, n_CT, n_ref_A, n_ref_C, n_reads
#
# region = "ercc"  for ERCC-prefixed chromosomes (sequencing-error baseline)
#        = "human" for all other chromosomes
#
# n_ref_A / n_ref_C are the denominators (reference A or C bases in the aligned
# region), enabling per-read editing fractions: n_AG / n_ref_A, n_CT / n_ref_C.
#
# Formula: n_ref_A = A_in_aligned_read - A_in_read_at_mismatch_positions
#                   + A_in_MD_mismatch_letters
# Correctly handles: edited A (read=G, MD=A) and G>A error (read=A, MD=G).
#
# Insertions (CIGAR I) are absent from MD and not corrected for;
# they are rare in STAR RNA-seq alignments.

args        <- commandArgs(trailingOnly = TRUE)
bam         <- args[1]
sample_name <- args[2]
output      <- args[3]

suppressPackageStartupMessages({
  library(Rsamtools)
  library(data.table)
})

#==============================================================================
# READ BAM
#==============================================================================

flag_filt <- scanBamFlag(
  isUnmappedQuery          = FALSE,
  isSecondaryAlignment     = FALSE,
  isSupplementaryAlignment = FALSE
)
param <- ScanBamParam(
  flag = flag_filt,
  what = c("rname", "seq", "cigar"),
  tag  = "MD"
)
raw <- scanBam(bam, param = param)[[1]]

rname <- as.character(raw$rname)
seq_v <- as.character(raw$seq)
cigar <- raw$cigar
md    <- raw$tag$MD

keep  <- !is.na(md) & nchar(seq_v) > 0L
rname <- rname[keep]; seq_v <- seq_v[keep]
cigar <- cigar[keep]; md    <- md[keep]

#==============================================================================
# SOFT-CLIP OFFSETS (vectorized)
#==============================================================================

leading_sc <- rep(0L, length(md))
has_lsc    <- grepl("^\\d+S", cigar, perl = TRUE)
leading_sc[has_lsc] <- as.integer(
  sub("^(\\d+)S.*", "\\1", cigar[has_lsc], perl = TRUE)
)

trailing_sc <- rep(0L, length(md))
has_tsc     <- grepl("\\d+S$", cigar, perl = TRUE)
trailing_sc[has_tsc] <- as.integer(
  sub(".*?(\\d+)S$", "\\1", cigar[has_tsc], perl = TRUE)
)

#==============================================================================
# COUNT SUBSTITUTIONS PER READ
#==============================================================================

# MD token types:
#   \d+      -> N matches: advance N aligned positions
#   [A-Z]    -> mismatch: ref base is this letter; advance 1 aligned position
#   \^[A-Z]+ -> deletion from reference: no read positions consumed
#
# Returns c(n_AG, n_CT, n_ref_A, n_ref_C)

count_one <- function(md_i, seq_i, lsc_i, tsc_i) {

  aln_seq <- substr(seq_i, lsc_i + 1L, nchar(seq_i) - tsc_i)
  n_A_aln <- nchar(gsub("[^A]", "", aln_seq))
  n_C_aln <- nchar(gsub("[^C]", "", aln_seq))

  toks <- regmatches(md_i, gregexpr("\\d+|[A-Z]|\\^[A-Z]+", md_i, perl = TRUE))[[1]]
  if (!length(toks)) return(c(0L, 0L, n_A_aln, n_C_aln))

  is_let <- nchar(toks) == 1L & grepl("^[A-Z]$", toks, perl = TRUE)
  steps  <- ifelse(grepl("^\\d+$", toks, perl = TRUE), as.integer(toks),
                   ifelse(is_let, 1L, 0L))

  let_idx <- which(is_let)
  if (!length(let_idx)) return(c(0L, 0L, n_A_aln, n_C_aln))

  cum_steps  <- c(0L, cumsum(steps))
  read_pos   <- cum_steps[let_idx] + lsc_i   # 0-based position in full seq

  ref_bases  <- toks[let_idx]
  read_bases <- substr(rep(seq_i, length(read_pos)), read_pos + 1L, read_pos + 1L)

  n_ag <- sum(ref_bases == "A" & read_bases == "G")
  n_ct <- sum(ref_bases == "C" & read_bases == "T")

  n_ref_A <- n_A_aln - sum(read_bases == "A") + sum(ref_bases == "A")
  n_ref_C <- n_C_aln - sum(read_bases == "C") + sum(ref_bases == "C")

  c(n_ag, n_ct, n_ref_A, n_ref_C)
}

res <- mapply(count_one, md, seq_v, leading_sc, trailing_sc, SIMPLIFY = TRUE)

#==============================================================================
# AGGREGATE AND WRITE
#==============================================================================

region <- ifelse(startsWith(rname, "ERCC"), "ercc", "human")
dt <- data.table(
  sample  = sample_name,
  region  = region,
  n_AG    = as.integer(res[1L, ]),
  n_CT    = as.integer(res[2L, ]),
  n_ref_A = as.integer(res[3L, ]),
  n_ref_C = as.integer(res[4L, ])
)
result <- dt[, .(n_reads = .N), by = .(sample, region, n_AG, n_CT, n_ref_A, n_ref_C)]
fwrite(result, output, sep = "\t")
message("Written: ", output)
