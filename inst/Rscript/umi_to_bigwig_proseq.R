#!/usr/bin/env Rscript

# Convert PRO-seq UMI counts into CPM-normalized, per-strand bigwig tracks.
#
# Usage:
#   Rscript umi_to_bigwig_proseq.R <umi.counts> <output.prefix>
# Writes <output.prefix>.ps.bw (plus strand) and <output.prefix>.ns.bw (minus).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Please specify:\n
       [required] 1/ UMI counts file\n
       [required] 2/ Output prefix (.ps.bw; .ns.bw)\n")
}

suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomicRanges, warn.conflicts = FALSE))

countsFile   <- args[1]
outputPrefix <- args[2]

# ---- Import UMI counts ----
counts <- fread(countsFile)
counts <- counts[coor != "NA:NA:NA"]
counts[, c("seqnames", "start", "strand") := tstrsplit(coor, ":", type.convert = T)]

# ---- Expand counts ----
counts <- counts[rep(seq(.N), umi_counts), .(seqnames, start, strand)]
counts[, end := start]

# ---- Resize reads to improve visualization (10 bp window) ----
counts[strand != "-", end := start + 9]
counts[strand == "-", start := end - 9]
setorderv(counts, c("seqnames", "start", "end"))

# ---- Compute coverage and export (CPM-normalized) ----
total_reads <- nrow(counts)
counts[, {
  gr <- GenomicRanges::GRanges(.SD)
  cov <- GenomicRanges::coverage(gr) / total_reads * 1e6
  outputFile <- paste0(outputPrefix, ifelse(strand == "+", ".ps.bw", ".ns.bw"))
  rtracklayer::export(cov, outputFile)
  print(outputFile)
}, strand]
