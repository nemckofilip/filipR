#!/usr/bin/env Rscript
# Called internally by fn_deseq2(). Do not run directly.
# Args: <count.files> <sample.names> <conditions> <ctl.conditions>
#       <padj.cutoff> <log2FC.cutoff> <dds.folder> <FC.folder>
#       <MA.folder> <output.prefix> <spikein.count.files|NULL>
#       <keep.genes.file|NULL>

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
spikein.files  <- if (args[11] == "NULL") NULL else unlist(tstrsplit(args[11], ","))
keep.genes.file <- if (args[12] == "NULL") NULL else args[12]

if (!any(!conditions %in% controls))
  stop("No non-control conditions found to compare against controls.")

#==============================================================================
# LOAD AND MERGE COUNT FILES
# When spike-in count files are supplied (e.g. per-ERCC counts), their rows are
# appended per sample so DESeq2 can use them as controlGenes for size factors.
#==============================================================================

dat <- lapply(counts, fread)
# Restrict gene rows to a functional keep-list (e.g. drop pseudogenes) before
# testing; spike-in rows are appended afterwards so they are always retained.
if (!is.null(keep.genes.file)) {
  keep <- fread(keep.genes.file, header = FALSE)$V1
  dat  <- lapply(dat, function(d) d[gene_id %in% keep])
}
if (!is.null(spikein.files)) {
  spike     <- lapply(spikein.files, fread)
  spike.ids <- unique(unlist(lapply(spike, function(x) x$gene_id)))
  dat       <- Map(rbind, dat, spike)
} else {
  spike.ids <- character(0)
}
names(dat) <- sample.names
dat <- rbindlist(dat, idcol = "sample")
dat[, sample := factor(sample, unique(sample))]
gene.names <- unique(dat, by = "gene_id")[, .(gene_id, gene_name)]
DF <- dcast(dat, gene_id ~ sample, value.var = "count")
DF <- data.frame(DF[, -1], row.names = DF$gene_id)

DF <- DF[rowSums(DF >= 3) >= 2, ]

# Spike-in rows surviving the detection filter are the control-gene set.
spike.idx <- which(rownames(DF) %in% spike.ids)
if (!is.null(spikein.files) && length(spike.idx) == 0)
  stop("Spike-in count files supplied but no spike-in rows survived filtering.")

#==============================================================================
# DESeq2
# With spike-ins, size factors come from DESeq2's median-of-ratios estimator
# restricted to the spike-in control genes; otherwise the default all-gene
# estimation is used.
#==============================================================================

sampleTable <- data.frame(condition = conditions, row.names = make.names(sample.names))
dds <- DESeqDataSetFromMatrix(countData = DF, colData = sampleTable, design = ~ condition)

if (length(spike.idx) > 0)
  dds <- estimateSizeFactors(dds, controlGenes = spike.idx)

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
    res <- res[!gene_id %in% spike.ids]
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
