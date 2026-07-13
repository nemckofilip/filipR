#!/usr/bin/env Rscript

# Count PRO-seq UMI reads over annotated features (strand-specific), with an
# optional blacklist subtraction. Outputs a two-column table (ID, count).
#
# Usage:
#   Rscript count_proseq_reads.R <umi.counts> <annotation.rds> <output.txt> [blacklist.rds]

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3, 4)) {
  stop("Please specify:\n
       [required] 1/ Reference genome UMI count file\n
       [required] 2/ A .rds annotation file\n
       [required] 3/ Output file name\n
       [optional] 4/ A .rds file of blacklisted regions (tRNAs...)\n")
}

suppressMessages(library(data.table))

umi_count  <- args[1]
annotation <- args[2]
outputFile <- args[3]
blacklist  <- if (length(args) == 4) args[4] else NULL

# ---- Import annotation ----
annot <- readRDS(annotation)

# ---- Import UMI counts ----
dat <- fread(umi_count)
dat[, c("seqnames", "start", "strand") := tstrsplit(coor, ":", type.convert = T)]

# ---- Remove reads overlapping blacklisted regions ----
if (!is.null(blacklist)) {
  blacklist <- data.table::copy(readRDS(blacklist))
  ov <- blacklist[dat, .N,
                  on = c("seqnames", "start<=start", "end>=start", "strand"),
                  .EACHI]$N
  stats <- data.table(
    N.removed.reads.clusters = sum(ov > 0, na.rm = T),
    N.removed.umi.counts = sum(dat[ov > 0, umi_counts], na.rm = T),
    perc.removed.umi.counts = sum(dat[ov > 0, umi_counts], na.rm = T) /
      sum(dat$umi_counts, na.rm = T) * 100
  )
  fwrite(stats, gsub(".txt$", "__blacklisted.txt", outputFile), sep = "\t", na = NA)
  dat <- dat[ov == 0]
}

# ---- Compute counts (strand-specific, read position within feature) ----
annot$count <- dat[annot, sum(umi_counts, na.rm = T),
                   on = c("seqnames", "start>=start", "start<=end", "strand"),
                   .EACHI]$V1
fwrite(annot[, .(ID = cluster.id, count)], outputFile, sep = "\t", na = NA)
