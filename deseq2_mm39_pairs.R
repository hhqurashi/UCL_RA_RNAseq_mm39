suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(EnhancedVolcano)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

## ---- EDIT THESE TWO ----
STAR_SALMON_DIR <- "/mnt/scratch/hqurashi/10_Sara_RNA_Seq/outputs/star_salmon"
OUT_BASE        <- "/mnt/scratch/hqurashi/10_Sara_RNA_Seq/outputs/DESeq2_pairs_mm39"
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

## ---- pick a counts file (prefer length_scaled) ----
## ---- pick a counts file (FOR DESeq2 prefer raw counts) ----
counts_tsv_counts <- file.path(STAR_SALMON_DIR, "salmon.merged.gene_counts.tsv")
counts_tsv_len    <- file.path(STAR_SALMON_DIR, "salmon.merged.gene_counts_length_scaled.tsv")

counts_tsv <- if (file.exists(counts_tsv_counts)) counts_tsv_counts else counts_tsv_len
if (!file.exists(counts_tsv)) stop("Could not find merged counts TSV in: ", STAR_SALMON_DIR)
message("Using counts file: ", counts_tsv)

## ---- load + keep ONLY expected sample columns ----
counts_raw <- read.delim(counts_tsv, check.names = FALSE)
gene_id <- counts_raw[[1]]
colnames(counts_raw)[1] <- "gene_id"

# keep only your 24 samples (prevents accidental non-sample columns)
samp_pat <- "^(WT|EM2|EM3|EM11)_F_9M_[0-9]+$"
sample_cols <- grep(samp_pat, colnames(counts_raw), value = TRUE)

if (length(sample_cols) == 0) stop("No sample columns matched pattern: ", samp_pat)

counts_mat <- as.matrix(counts_raw[, sample_cols, drop = FALSE])

# coerce safely
counts <- suppressWarnings(apply(counts_mat, 2, as.numeric))
rownames(counts) <- gene_id

# replace NA with 0, enforce non-negative integers
na_n <- sum(is.na(counts))
message("NAs after numeric coercion: ", na_n)
counts[is.na(counts)] <- 0
counts <- round(counts)
counts[counts < 0] <- 0

# sanity: drop any all-zero samples (otherwise size factors become NA)
zero_samp <- colnames(counts)[colSums(counts) == 0]
if (length(zero_samp) > 0) {
  stop("These samples have total count = 0: ", paste(zero_samp, collapse = ", "))
}

samples <- colnames(counts)

## ---- infer group from sample name (strip trailing _<number>) ----
group <- sub("_[0-9]+$", "", samples)

meta <- data.frame(
  sample = samples,
  group  = factor(group),
  row.names = samples,
  stringsAsFactors = FALSE
)

# sanity: expected groups
expected <- c("WT_F_9M", "EM2_F_9M", "EM3_F_9M", "EM11_F_9M")
missing_groups <- setdiff(expected, levels(meta$group))
if (length(missing_groups) > 0) {
  warning("Missing expected groups: ", paste(missing_groups, collapse = ", "),
          "\nFound: ", paste(levels(meta$group), collapse = ", "))
}

# set WT as reference if present
if ("WT_F_9M" %in% levels(meta$group)) {
  meta$group <- relevel(meta$group, "WT_F_9M")
}

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ group)

# drop genes with all zeros
dds <- dds[rowSums(counts(dds)) > 10, ]

# HARD force size factors using poscounts
dds <- estimateSizeFactors(dds, type = "poscounts")

# Now run DESeq without re-estimating size factors
dds <- DESeq(dds)

## ---- overall PCA (all samples) ----
vsd <- vst(dds, blind = FALSE)
p_pca <- plotPCA(vsd, intgroup = "group") + theme_classic()

ggsave(file.path(OUT_BASE, "PCA_all_samples.png"), p_pca, width = 6, height = 5, dpi = 300)
ggsave(file.path(OUT_BASE, "PCA_all_samples.pdf"), p_pca, width = 6, height = 5)

## ---- gene labels: Ensembl -> SYMBOL (mouse) ----
ens <- sub("\\..*$", "", rownames(dds))  # strip version if present
symbol <- mapIds(org.Mm.eg.db,
                 keys = ens,
                 keytype = "ENSEMBL",
                 column = "SYMBOL",
                 multiVals = "first")

## ---- helper: volcano ----
make_volcano <- function(res_df, outfile_png, title) {
  res_df$padj[is.na(res_df$padj)] <- 1
  res_df$padj <- pmax(res_df$padj, .Machine$double.xmin)

  png(outfile_png, width = 3600, height = 3000, res = 300)
  print(
    EnhancedVolcano(
      res_df,
      lab      = res_df$label,
      x        = "log2FoldChange",
      y        = "padj",
      pCutoff  = 0.05,
      FCcutoff = 1,
      title    = title
    )
  )
  dev.off()
}

## ---- define your 6 contrasts (g2 vs g1) ----
comparisons <- list(
  "WT_F_9M_vs_EM2_F_9M"   = c("WT_F_9M",   "EM2_F_9M"),
  "WT_F_9M_vs_EM3_F_9M"   = c("WT_F_9M",   "EM3_F_9M"),
  "WT_F_9M_vs_EM11_F_9M"  = c("WT_F_9M",   "EM11_F_9M"),
  "EM2_F_9M_vs_EM11_F_9M" = c("EM2_F_9M",  "EM11_F_9M"),
  "EM2_F_9M_vs_EM3_F_9M"  = c("EM2_F_9M",  "EM3_F_9M"),
  "EM11_F_9M_vs_EM3_F_9M" = c("EM11_F_9M", "EM3_F_9M")
)

## ---- run contrasts ----
for (nm in names(comparisons)) {
  g1 <- comparisons[[nm]][1]
  g2 <- comparisons[[nm]][2]

  if (!all(c(g1, g2) %in% levels(meta$group))) {
    warning("Skipping ", nm, " because one/both groups not found in data.")
    next
  }

  message("Running: ", nm, " (", g2, " vs ", g1, ")")

  res <- results(dds, contrast = c("group", g2, g1))
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  ens_res <- sub("\\..*$", "", res_df$gene_id)
  res_df$symbol <- unname(symbol[match(ens_res, names(symbol))])
  res_df$label  <- ifelse(is.na(res_df$symbol) | res_df$symbol == "", ens_res, res_df$symbol)

  res_df <- res_df[order(res_df$padj), ]

  out_dir <- file.path(OUT_BASE, nm)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(res_df, file.path(out_dir, paste0(nm, "_DE_results.csv")), row.names = FALSE)

  make_volcano(
    res_df,
    file.path(out_dir, paste0(nm, "_volcano.png")),
    title = nm
  )
}

message("Done. Outputs in: ", OUT_BASE)