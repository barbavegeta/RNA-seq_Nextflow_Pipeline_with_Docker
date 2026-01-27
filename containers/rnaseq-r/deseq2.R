#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(optparse)
  library(ggplot2)
})

opt_list <- list(
  make_option("--counts", type="character"),
  make_option("--design", type="character"),
  make_option("--out_results", type="character"),
  make_option("--out_ma", type="character"),
  make_option("--out_pca", type="character")
)
opt <- parse_args(OptionParser(option_list=opt_list))

stopifnot(file.exists(opt$counts), file.exists(opt$design))

counts <- fread(opt$counts)
gene_id <- counts[[1]]
count_mat <- as.matrix(counts[, -1, with=FALSE])
rownames(count_mat) <- gene_id

design <- fread(opt$design)
design <- as.data.frame(design)
stopifnot(all(design$sample %in% colnames(count_mat)))
rownames(design) <- design$sample
design <- design[colnames(count_mat), , drop=FALSE]

dds <- DESeqDataSetFromMatrix(countData=count_mat, colData=design, design=~ condition)
dds <- DESeq(dds)

res <- results(dds)
res_dt <- as.data.table(res, keep.rownames="gene_id")
fwrite(res_dt, opt$out_results, sep="\t")

png(opt$out_ma, width=900, height=700)
plotMA(res, ylim=c(-5,5))
dev.off()

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

ggsave(opt$out_pca, p, width=7, height=5, dpi=150)
