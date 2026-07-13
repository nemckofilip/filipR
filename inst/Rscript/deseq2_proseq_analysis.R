#!/usr/bin/env Rscript

# DESeq2 differential expression for PRO-seq count tables. No spike-in: size
# factors come from library size (umi_counts, "libSize") or DESeq2 defaults.
#
# Usage:
#   Rscript deseq2_proseq_analysis.R <counts> <refStats> <names> <conditions> \
#       <controls> <dds.dir> <FC.dir> <PDF.dir> <experiment> <feature> <norm>
# where <counts>, <refStats>, <names>, <conditions>, <controls> are comma-separated.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11) {
  stop("Please specify:\n
       [required]  1/ Comma-separated count tables (ref genome)\n
       [required]  2/ Comma-separated read statistics (ref genome)\n
       [required]  3/ Comma-separated sample names\n
       [required]  4/ Comma-separated conditions\n
       [required]  5/ Comma-separated controls\n
       [required]  6/ dds output folder\n
       [required]  7/ FC tables output folder\n
       [required]  8/ PDF output folder\n
       [required]  9/ Experiment\n
       [required] 10/ Feature\n
       [required] 11/ Normalization: 'default' or 'libSize'\n")
}

suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(DESeq2, warn.conflicts = FALSE))

# ---- Parse variables ----
counts     <- unlist(tstrsplit(args[1], ","))
refStats   <- unlist(tstrsplit(args[2], ","))
names      <- gsub("-", ".", unlist(tstrsplit(args[3], ",")))
conditions <- gsub("-", ".", unlist(tstrsplit(args[4], ",")))
controls   <- gsub("-", ".", unlist(tstrsplit(args[5], ",")))
dds_output_folder <- args[6]
FC_output_folder  <- args[7]
PDF_output_folder <- args[8]
experiment <- args[9]
feature    <- args[10]
norm       <- args[11]

# Names and conditions must not start with a number ----
if (any(grepl("^[0-9]", names)))
  names[grepl("^[0-9]", names)] <- paste0("X", grep("^[0-9]", names, value = T))
if (any(grepl("^[0-9]", conditions)))
  conditions[grepl("^[0-9]", conditions)] <- paste0("X", grep("^[0-9]", conditions, value = T))

# ---- Import count tables ----
dat <- lapply(counts, fread)
names(dat) <- names
dat <- rbindlist(dat, idcol = "condition")
dat[, condition := factor(condition, unique(condition))]
DF <- dcast(dat, ID ~ condition, value.var = "count")
DF <- data.frame(DF[, -1], row.names = DF$ID)

# ---- Remove low count reads ----
DF <- DF[rowSums(DF >= 3) >= 2, ]

# ---- Sample table ----
sampleTable <- data.frame(condition = conditions, row.names = names)

# ---- Reference genome read counts ----
ref <- lapply(refStats, fread)
names(ref) <- names
ref <- rbindlist(ref, idcol = "sample")
ref[, cdition := conditions]

# ---- DESeq2 analysis ----
print(paste("Start", norm, "normalization"))
dds <- DESeqDataSetFromMatrix(countData = DF,
                              colData = sampleTable,
                              design = ~ condition)

# ---- Size factors ----
if (norm == "libSize") {
  sizeFactors(dds) <- ref[, umi_counts / median(umi_counts)]
} else if (norm != "default") {
  stop("normalization should be one of 'default' or 'libSize'")
}

# ---- Compute model and save object ----
dds <- DESeq(dds)
saveRDS(dds, paste0(dds_output_folder, "/", experiment, "_", feature, "_", norm, "_norm_DESeq2.dds"))

# ---- Open pdf to save MA plots ----
outputPdf <- paste0(PDF_output_folder, "/", experiment, "_", feature, "_", norm, "_norm_MAplots.pdf")
pdf(outputPdf, 4, 3)
par(mai = c(.9, 1.5, .9, 1.3), cex.axis = 6 / 12, cex.lab = 7 / 12, las = 1,
    tcl = -0.1, mgp = c(.8, 0.25, 0), font.main = 1, cex.main = 9 / 12)

FC <- list()

# ---- For each control and each condition: compute FC + MA plot ----
for (ctl in unique(controls)) {
  for (cdition in setdiff(unique(conditions), ctl)) {
    res <- results(dds, contrast = c("condition", cdition, ctl))
    res <- as.data.table(as.data.frame(res), keep.rownames = "gene_id")
    res[, diff := fcase(padj < 0.05 & log2FoldChange > log2(1.5), "Up-regulated",
                        padj < 0.05 & log2FoldChange < (-log2(1.5)), "Down-regulated",
                        default = "Unaffected")]
    res[, condition := cdition]
    res[, control := ctl]
    FC[[paste0(cdition, "_vs_", ctl)]] <- data.table::copy(res)

    res[, col := fcase(diff == "Up-regulated", "tomato",
                       diff == "Down-regulated", "cornflowerblue",
                       default = "lightgrey")]
    res <- res[order(col == "lightgrey", decreasing = TRUE)]
    lims <- quantile(res$log2FoldChange, c(0.001, 0.999), na.rm = TRUE)
    res[, pch := ifelse(between(log2FoldChange, lims[1], lims[2]), 16, 17)]
    res[log2FoldChange < lims[1], log2FoldChange := lims[1]]
    res[log2FoldChange > lims[2], log2FoldChange := lims[2]]

    res[, {
      plot(x = log10(baseMean), y = log2FoldChange,
           col = adjustcolor(col, .5), pch = pch,
           ylab = "Fold change (log2)", frame = F, xaxt = "n", cex = .5,
           main = paste(cdition, "vs.", ctl))
      axis(1, padj = -1.45)
      abline(h = 0, lty = 3)
      nUp <- formatC(sum(diff == "Up-regulated"), big.mark = ",")
      nDown <- formatC(sum(diff == "Down-regulated"), big.mark = ",")
      legend(par("usr")[2], par("usr")[4],
             col = adjustcolor(c("tomato", "cornflowerblue"), 0.5),
             legend = c(paste0("Up-regulated (", nUp, ")"),
                        paste0("Down-regulated (", nDown, ")")),
             pch = 16, bty = "n", cex = 6 / 12, border = NA, xpd = NA)
    }]
  }
}
dev.off()

# ---- Save FC file ----
FC <- rbindlist(FC)
setcolorder(FC, c("condition", "control"))
fwrite(FC, paste0(FC_output_folder, "/", experiment, "_", feature, "_", norm, "_norm_DESeq2_FC.txt"),
       sep = "\t", na = NA)
