library(EnhancedVolcano)
library(AnnotationDbi); library(org.Hs.eg.db)

res <- read.csv("/media/pontikos_nas2/HarisQurashi/projects/9_Owen_BulkRNAseq/test/DE_results.csv", row.names=1, check.names=FALSE)
res$padj[is.na(res$padj)] <- 1; res$padj <- pmax(res$padj, .Machine$double.xmin)
res$symbol <- mapIds(org.Hs.eg.db, sub("\\..*$","",rownames(res)), "SYMBOL", "ENSEMBL", multiVals="first")

png("/media/pontikos_nas2/HarisQurashi/projects/9_Owen_BulkRNAseq/test/volcano_enhanced_symbols.png", 3600, 3000, res=300)
print(EnhancedVolcano(res, lab=ifelse(is.na(res$symbol), rownames(res), res$symbol),
                     x="log2FoldChange", y="padj", pCutoff=0.05, FCcutoff=1))
dev.off()