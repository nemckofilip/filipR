#!/usr/bin/env Rscript

# Collapse UMIs from a PRO-seq BAM and count unique molecules per single-
# nucleotide genomic position. Adapted for a configurable UMI length.
#
# Usage:
#   Rscript umi_collapse_proseq.R <bam> <counts.out> <stats.out> <flip.strand> <umi.length>

# ---- Arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Please specify:\n
       [required] 1/ Aligned bam file containing UMIs\n
       [required] 2/ Output count file (.txt)\n
       [required] 3/ Output stats file (.txt)\n
       [required] 4/ Should the strand of the read be flipped? TRUE for PRO-seq\n
       [required] 5/ UMI length (nt)\n")
}

suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(Rsamtools, warn.conflicts = FALSE))
suppressMessages(library(stringdist, warn.conflicts = FALSE))

bam           <- args[1]
counts.output <- args[2]
stats.output  <- args[3]
flip.strand   <- as.logical(args[4])
umi.length    <- as.integer(args[5])

# ---- Import data ----
param <- Rsamtools::ScanBamParam(what = c("qname", "rname", "pos", "strand", "qwidth", "mapq"))
dat <- Rsamtools::scanBam(bam, param = param)[[1]]
dat <- as.data.table(dat)
setnames(dat,
         c("qname", "rname", "pos", "strand", "qwidth", "mapq"),
         c("read", "seqnames", "start", "strand", "width", "mapq"))

# ---- Resize to single nucleotide (3' end) ----
dat[, strand := as.character(strand)]
dat[strand == "-", start := start + width - 1]

# ---- Flip the strand of the reads (PRO-seq: TRUE) ----
if (flip.strand)
  dat[!is.na(strand), strand := ifelse(strand == "+", "-", "+")]

# ---- Save coordinates ----
dat[, coor := paste0(seqnames, ":", start, ":", strand)]

# ---- Extract UMIs ----
dat[!is.na(seqnames), UMI := sub(".*_(.*)$", "\\1", read)]
dat[is.na(seqnames), UMI := strrep("G", umi.length)] # Unaligned reads (kept for total reads)
umi_length <- nchar(dat$UMI)
if (any(umi_length < umi.length)) {
  stop(sprintf("Some UMIs are shorter than %d nt", umi.length))
} else if (any(umi_length > umi.length)) {
  warning(sprintf("Some UMIs were longer than %d nt and will be trimmed.", umi.length))
  dat[umi_length > umi.length, UMI := substr(UMI, 1, umi.length)]
}

# ---- Compute total counts ----
dat <- dat[, .(umi_N = .N), .(coor, UMI)]
dat[, total_counts := sum(umi_N), coor]
setorderv(dat, "umi_N", order = -1)

# ---- Check whether UMIs might be collapsed ----
dat[, collapsed := T, coor]
dat[, idx := .I]
for (i in 1:umi.length) {
  dat[, check := idx[1], .(coor, gsub(paste0("^(.{", i - 1, "})."), "\\1", UMI))]
  potentialDup <- unique(dat[(check < idx), c(check, idx)])
  dat[potentialDup, collapsed := FALSE]
  print(i)
}
dat$idx <- NULL
paste0(sum(dat$collapsed), " / ", nrow(dat), " pre-collapsed")

# ---- UMI collapsing (<= 1 mismatch) ----
while (any(!dat$collapsed)) {
  dat[!(collapsed), c("collapsed", "UMI") := {
    coll <- stringdist(UMI[1],
                       UMI,
                       method = "hamming",
                       nthread = getDTthreads() - 1) <= 1
    UMI[coll] <- UMI[1]
    .(coll, UMI)
  }, coor]
}

# ---- Final collapsing ----
dat <- unique(dat[, .(coor, total_counts, UMI)])
dat <- dat[, .(umi_counts = .N), .(coor, total_counts)]
dat[coor == "NA:NA:NA", umi_counts := NA]

# ---- Save output ----
fwrite(dat, counts.output, sep = "\t", quote = F, na = NA)

# ---- Compute statistics ----
stats <- data.table(total = sum(dat$total_counts),
                    mapped = sum(dat[coor != "NA:NA:NA", total_counts]),
                    umi_counts = sum(dat[coor != "NA:NA:NA", umi_counts]))
fwrite(stats, stats.output, sep = "\t", na = NA)
