library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(AnnotationDbi)
library(org.Hs.eg.db)

## ---------- Paths ----------
run1_counts <- "/media/pontikos_nas2/HarisQurashi/projects/9_Owen_BulkRNAseq/outputs/1_HAU_controls/star_salmon/salmon.merged.gene_counts.tsv"
run2_counts <- "/media/pontikos_nas2/HarisQurashi/projects/9_Owen_BulkRNAseq/outputs/2_P_T_UT_OF_INTEREST/star_salmon/salmon.merged.gene_counts.tsv"

# directory where you want the 1_..9_.. folders
base_out   <- "/media/pontikos_nas2/HarisQurashi/projects/9_Owen_BulkRNAseq/outputs/2_P_T_UT_OF_INTEREST"
dir.create(base_out, showWarnings = FALSE, recursive = TRUE)

## ---------- Sample -> group mapping (define FIRST) ----------
group_map <- c(
  # controls
  "HAU2001" = "HAU200",
  "HAU2002" = "HAU200",
  "HAU2003" = "HAU200",

  # P1AU210
  "P1AU2101"      = "P1AU210",
  "P1AU2102"      = "P1AU210",
  "P1AU2103_mix"  = "P1AU210",
  "P1AU2104"      = "P1AU210",

  # P2AES210
  "P2AES2101"     = "P2AES210",
  "P2AES2102"     = "P2AES210",
  "P2AES2103"     = "P2AES210",

  # P3AEU210  (P3AEU2102 dropped)
  "P3AEU2101"     = "P3AEU210",
  "P3AEU2103_mix" = "P3AEU210",
  "P3AEU2104_mix" = "P3AEU210",

  # P1AA210
  "P1AA2101_mix"  = "P1AA210",
  "P1AA2102"      = "P1AA210",
  "P1AA2103_mix"  = "P1AA210",
  "P1AA2104_mix"  = "P1AA210",

  # P2AEA210
  "P2AEA2101"     = "P2AEA210",
  "P2AEA2102"     = "P2AEA210",
  "P2AEA2103"     = "P2AEA210",

  # P3AEA210
  "P3AEA2101_mix" = "P3AEA210",
  "P3AEA2102_mix" = "P3AEA210",
  "P3AEA2103_mix" = "P3AEA210",
  "P3AEA2104_mix" = "P3AEA210"
)

## ---------- Load & combine counts (only keep sample columns) ----------
c1 <- read.delim(run1_counts, check.names = FALSE)
c2 <- read.delim(run2_counts, check.names = FALSE)

# same genes in same order
stopifnot(identical(c1[[1]], c2[[1]]))
genes <- c1[[1]]

sample_ids <- names(group_map)

# columns in each file that correspond to our samples of interest
samples_run1 <- intersect(sample_ids, colnames(c1))
samples_run2 <- intersect(sample_ids, colnames(c2))

counts1 <- as.matrix(c1[, samples_run1, drop = FALSE])
counts2 <- as.matrix(c2[, samples_run2, drop = FALSE])

rownames(counts1) <- genes
rownames(counts2) <- genes

counts_all <- cbind(counts1, counts2)
samples <- colnames(counts_all)

## ---------- Metadata ----------
group <- unname(group_map[samples])

meta_all <- data.frame(
  sample   = samples,
  group    = factor(group),
  row.names = samples
)

## ---------- Define the 9 contrasts ----------
comparisons <- list(
  "1_HAU200_vs_P1AU210"    = c("HAU200",  "P1AU210"),
  "2_HAU200_vs_P2AES210"   = c("HAU200",  "P2AES210"),
  "3_HAU200_vs_P3AEU210"   = c("HAU200",  "P3AEU210"),
  "4_P1AU210_vs_P1AA210"   = c("P1AU210", "P1AA210"),
  "5_P2AES210_vs_P2AEA210" = c("P2AES210","P2AEA210"),
  "6_P3AEU210_vs_P3AEA210" = c("P3AEU210","P3AEA210"),
  "7_HAU200_vs_P1AA210"    = c("HAU200",  "P1AA210"),
  "8_HAU200_vs_P2AEA210"   = c("HAU200",  "P2AEA210"),
  "9_HAU200_vs_P3AEA210"   = c("HAU200",  "P3AEA210")
)

## ---------- Helper for volcano (labels already gene symbols) ----------
make_volcano <- function(res, outfile) {
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)  # already gene symbols

  res_df$padj[is.na(res_df$padj)] <- 1
  res_df$padj <- pmax(res_df$padj, .Machine$double.xmin)

  png(outfile, width = 3600, height = 3000, res = 300)
  print(
    EnhancedVolcano(
      res_df,
      lab      = res_df$gene_id,
      x        = "log2FoldChange",
      y        = "padj",
      pCutoff  = 0.05,
      FCcutoff = 1
    )
  )
  dev.off()
}

## ---------- Loop over contrasts ----------
for (comp_name in names(comparisons)) {

  groups_pair <- comparisons[[comp_name]]
  g1 <- groups_pair[1]   # reference
  g2 <- groups_pair[2]   # treatment

  message("Running DESeq2 for: ", comp_name, " (", g2, " vs ", g1, ")")

  # subset samples for these two groups
  keep_samp <- meta_all$group %in% c(g1, g2)
  counts_sub <- counts_all[, keep_samp, drop = FALSE]
  meta_sub   <- droplevels(meta_all[keep_samp, ])

  # drop genes with 0 total counts in this subset
  keep_genes <- rowSums(counts_sub) > 0
  counts_sub <- counts_sub[keep_genes, ]

  # round Salmon estimates to integers for DESeq2
  counts_sub <- round(counts_sub)
  counts_sub[counts_sub < 0] <- 0

  meta_sub$condition <- factor(meta_sub$group, levels = c(g1, g2))

  dds <- DESeqDataSetFromMatrix(countData = counts_sub,
                                colData   = meta_sub,
                                design    = ~ condition)

  dds <- DESeq(dds)

  res <- results(dds, contrast = c("condition", g2, g1))

  out_dir <- file.path(base_out, comp_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  ## ---- Prefix filenames with comparison name ----
  norm_counts <- counts(dds, normalized = TRUE)
  write.csv(
    as.data.frame(norm_counts),
    file = file.path(out_dir, paste0(comp_name, "_normalized_counts.csv"))
  )

  res_ord <- res[order(res$padj), ]
  write.csv(
    as.data.frame(res_ord),
    file = file.path(out_dir, paste0(comp_name, "_DE_results.csv"))
  )

  vsd <- vst(dds, blind = TRUE)
  p_pca <- plotPCA(vsd, intgroup = "condition") + theme_classic()
  ggsave(
    file.path(out_dir, paste0(comp_name, "_PCA.png")),
    p_pca, width = 6, height = 5, dpi = 300
  )

  make_volcano(
    res_ord,
    file.path(out_dir, paste0(comp_name, "_volcano_enhanced_symbols.png"))
  )
}