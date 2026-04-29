library(DESeq2)
library(ggplot2)

cts <- read.delim("/media/pontikos_nas2/HarisQurashi/projects/9_Owen_BulkRNAseq/test/counts.txt", comment.char="#", check.names=FALSE)
counts <- as.matrix(cts[,-(1:6)]); rownames(counts) <- cts$Geneid
colnames(counts) <- sub("\\.markdup\\.sorted\\.bam$", "", basename(colnames(counts)))

meta <- read.delim("/media/pontikos_nas2/HarisQurashi/projects/9_Owen_BulkRNAseq/test/meta.tsv"); rownames(meta) <- meta$sample
setdiff(rownames(meta), colnames(counts)); setdiff(colnames(counts), rownames(meta))

dds <- DESeqDataSetFromMatrix(counts[, rownames(meta), drop=FALSE], meta, ~ condition)

dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized=TRUE)

write.csv(as.data.frame(normalized_counts), "/media/pontikos_nas2/HarisQurashi/projects/9_Owen_BulkRNAseq/test/normalized_counts.csv")

res <- results(dds)
write.csv(as.data.frame(res[order(res$padj),]), "/media/pontikos_nas2/HarisQurashi/projects/9_Owen_BulkRNAseq/test/DE_results.csv")

vsd <- vst(dds, blind=TRUE)
p <- plotPCA(vsd, intgroup="condition") + theme_classic()
ggsave("/media/pontikos_nas2/HarisQurashi/projects/9_Owen_BulkRNAseq/test/PCA.png", p, width=6, height=5, dpi=300)