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
  make_option("--formula", type="character", default="~ condition",
              help="DESeq2 design formula; the tested variable must be 'condition' and come last"),
  make_option("--reference", type="character", default=NULL,
              help="reference level of 'condition' (default: first level alphabetically)"),
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

# every covariate becomes a factor; 'condition' gets an explicit reference level so the
# sign of log2FoldChange is always (other level) vs (reference)
for (v in setdiff(colnames(design), "sample")) design[[v]] <- factor(design[[v]])
if (!is.null(opt$reference)) {
  stopifnot(opt$reference %in% levels(design$condition))
  design$condition <- relevel(design$condition, ref = opt$reference)
}
ref_level <- levels(design$condition)[1]
message("DESeq2 design: ", opt$formula, " | condition reference level: ", ref_level)

dds <- DESeqDataSetFromMatrix(countData=count_mat, colData=design, design=as.formula(opt$formula))
dds <- DESeq(dds)

if (nlevels(design$condition) == 2) {
  test_level <- levels(design$condition)[2]
  res <- results(dds, contrast = c("condition", test_level, ref_level))
  message("Contrast: ", test_level, " vs ", ref_level)
} else {
  res <- results(dds)
}
res_dt <- as.data.table(res, keep.rownames="gene_id")
fwrite(res_dt, opt$out_results, sep="\t")

png(opt$out_ma, width=900, height=700)
plotMA(res, ylim=c(-5,5))
dev.off()

n <- nrow(dds)
if (n < 1000) {
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
} else {
  vsd <- vst(dds, blind = FALSE)
}
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

ggsave(opt$out_pca, p, width=7, height=5, dpi=150)