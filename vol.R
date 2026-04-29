library(EnhancedVolcano)

res <- read.csv("/media/pontikos_nas2/HarisQurashi/projects/9_Owen_BulkRNAseq/test/DE_results.csv", row.names=1, check.names=FALSE)
res$padj[is.na(res$padj)] <- 1
res$padj <- pmax(res$padj, .Machine$double.xmin)

p <- EnhancedVolcano(res,
  lab = rownames(res),          
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1
)

png("/media/pontikos_nas2/HarisQurashi/projects/9_Owen_BulkRNAseq/test/volcano_enhanced.png", width=3600, height=3000, res=300)
print(p)
dev.off()