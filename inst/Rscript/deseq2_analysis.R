#!/usr/bin/env Rscript
# Called internally by fn_deseq2(). Do not run directly.
# Args: <count.files> <sample.names> <conditions> <ctl.conditions>
#       <padj.cutoff> <log2FC.cutoff> <dds.folder> <FC.folder>
#       <MA.folder> <output.prefix> <norm.counts|NULL>

suppressMessages(library(data.table))
suppressMessages(library(DESeq2))

args           <- commandArgs(trailingOnly = TRUE)
counts         <- unlist(tstrsplit(args[1], ","))
sample.names   <- unlist(tstrsplit(args[2], ","))
sample.names   <- gsub("-", ".", sample.names)
conditions     <- unlist(tstrsplit(args[3], ","))
conditions     <- gsub("-", ".", conditions)
controls       <- unique(gsub("-", ".", unlist(tstrsplit(args[4], ","))))
padj.cutoff    <- as.numeric(args[5])
log2FC.cutoff  <- as.numeric(args[6])
dds.folder     <- args[7]
FC.folder      <- args[8]
MA.folder      <- args[9]
output.prefix  <- args[10]
norm.counts    <- if (args[11] == "NULL") NULL else as.numeric(unlist(tstrsplit(args[11], ",")))

if (!any(!conditions %in% controls))
  stop("No non-control conditions found to compare against controls.")

#==============================================================================
# LOAD AND MERGE COUNT FILES
#==============================================================================

dat <- lapply(counts, fread)
names(dat) <- sample.names
dat <- rbindlist(dat, idcol = "sample")
dat[, sample := factor(sample, unique(sample))]
gene.names <- unique(dat, by = "gene_id")[, .(gene_id, gene_name)]
DF <- dcast(dat, gene_id ~ sample, value.var = "count")
DF <- data.frame(DF[, -1], row.names = DF$gene_id)

DF <- DF[rowSums(DF >= 3) >= 2, ]

#==============================================================================
# DESeq2
#==============================================================================

sampleTable <- data.frame(condition = conditions, row.names = make.names(sample.names))
dds <- DESeqDataSetFromMatrix(countData = DF, colData = sampleTable, design = ~ condition)

if (!is.null(norm.counts))
  sizeFactors(dds) <- norm.counts / min(norm.counts)

dds <- DESeq(dds)

dir.create(dds.folder, recursive = TRUE, showWarnings = FALSE)
saveRDS(dds, file.path(dds.folder, paste0(output.prefix, "_DESeq2.dds")))

#==============================================================================
# FC TABLES + MA PLOTS
#==============================================================================

dir.create(FC.folder, recursive = TRUE, showWarnings = FALSE)
dir.create(MA.folder, recursive = TRUE, showWarnings = FALSE)

pdf(file.path(MA.folder, paste0(output.prefix, "_MAplots.pdf")), 4, 3)
par(mai = c(.9, 1.5, .9, 1.3), cex.axis = 6/12, cex.lab = 7/12, las = 1,
    tcl = -0.1, mgp = c(.8, 0.25, 0), font.main = 1, cex.main = 9/12)

FC <- list()

for (ctl in unique(controls)) {
  for (cdition in setdiff(unique(conditions), ctl)) {

    res <- as.data.table(as.data.frame(results(dds, contrast = c("condition", cdition, ctl))),
                         keep.rownames = "gene_id")
    res[, diff := fcase(padj < padj.cutoff & log2FoldChange >   log2FC.cutoff,  "Up-regulated",
                        padj < padj.cutoff & log2FoldChange < (-log2FC.cutoff), "Down-regulated",
                        default = "Unaffected")]
    res[, condition := cdition]
    res[, control   := ctl]
    FC[[paste0(cdition, "_vs_", ctl)]] <- copy(res)

    res[, col := fcase(diff == "Up-regulated",   "tomato",
                       diff == "Down-regulated", "cornflowerblue",
                       default = "lightgrey")]
    res <- res[order(col == "lightgrey", decreasing = TRUE)]

    lims <- quantile(res$log2FoldChange, c(0.001, 0.999), na.rm = TRUE)
    res[, pch := ifelse(between(log2FoldChange, lims[1], lims[2]), 16, 17)]
    res[log2FoldChange < lims[1], log2FoldChange := lims[1]]
    res[log2FoldChange > lims[2], log2FoldChange := lims[2]]

    res[, {
      plot(log10(baseMean), log2FoldChange,
           col  = adjustcolor(col, .5), pch = pch, cex = .5,
           ylab = "Fold change (log2)", frame = FALSE, xaxt = "n",
           main = paste(cdition, "vs.", ctl))
      axis(1, padj = -1.45)
      abline(h = 0, lty = 3)
      nUp   <- formatC(sum(diff == "Up-regulated"),   big.mark = ",")
      nDown <- formatC(sum(diff == "Down-regulated"), big.mark = ",")
      legend(par("usr")[2], par("usr")[4],
             col    = adjustcolor(c("tomato", "cornflowerblue"), 0.5),
             legend = c(paste0("Up-regulated (", nUp, ")"),
                        paste0("Down-regulated (", nDown, ")")),
             pch = 16, bty = "n", cex = 6/12, border = NA, xpd = NA)
    }]
  }
}
dev.off()

FC <- rbindlist(FC)
FC <- merge(FC, gene.names, by = "gene_id", all.x = TRUE, sort = FALSE)
setcolorder(FC, c("condition", "control", "gene_id", "gene_name"))
fwrite(FC, file.path(FC.folder, paste0(output.prefix, "_DESeq2_FC.txt")), sep = "\t", na = "NA")

cat("Done:", output.prefix, "\n")
