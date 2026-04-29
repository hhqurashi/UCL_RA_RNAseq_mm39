#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)          # PCA label style
  library(EnhancedVolcano)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

## -----------------------------
## USER EDITS
## -----------------------------
STAR_SALMON_DIR <- "/mnt/scratch/hqurashi/10_Sara_RNA_Seq/outputs/star_salmon"
OUT_BASE        <- "/mnt/scratch/hqurashi/10_Sara_RNA_Seq/outputs/DESeq2_pairs_mm39_tximport"
GENE_TPM_XLSX <- "/media/pontikos_nas2/HarisQurashi/projects/10_Sara_BulkRNAseq/outputs/star_salmon/gene_tpm_counts.xlsx"
SAMP_PAT <- "^(WT|EM2|EM3|EM11)_F_9M_[0-9]+$"

MIN_TOTAL_COUNT <- 10
ALPHA           <- 0.05
LFC_CUTOFF      <- 1

dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

## -----------------------------
## Locate quant.sf files
## -----------------------------
sample_dirs <- list.dirs(STAR_SALMON_DIR, recursive = FALSE, full.names = TRUE)
if (length(sample_dirs) == 0) stop("No sample subdirectories found in: ", STAR_SALMON_DIR)

pick_quant <- function(d) {
  cands <- c(
    file.path(d, "quant.sf"),
    file.path(d, "salmon", "quant.sf"),
    file.path(d, "logs", "quant.sf")
  )
  cands[file.exists(cands)][1]
}

quant_files <- vapply(sample_dirs, pick_quant, character(1))
sample_names <- basename(sample_dirs)

ok <- !is.na(quant_files) & nzchar(quant_files)
quant_files <- quant_files[ok]
sample_names <- sample_names[ok]

keep <- grepl(SAMP_PAT, sample_names)
quant_files <- quant_files[keep]
sample_names <- sample_names[keep]

if (length(quant_files) == 0) {
  stop("No quant.sf files matched pattern: ", SAMP_PAT,
       "\nCheck STAR_SALMON_DIR and sample directory names.")
}

files <- setNames(quant_files, sample_names)

message("Found ", length(files), " samples with quant.sf.")
print(table(sub("_[0-9]+$", "", names(files))))

## -----------------------------
## Locate tx2gene.tsv
## -----------------------------
tx2gene_path <- file.path(STAR_SALMON_DIR, "tx2gene.tsv")
if (!file.exists(tx2gene_path)) {
  hits <- list.files(STAR_SALMON_DIR, pattern = "^tx2gene\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("Could not find tx2gene.tsv under: ", STAR_SALMON_DIR)
  tx2gene_path <- hits[1]
}
message("Using tx2gene: ", tx2gene_path)

tx2gene <- read.delim(tx2gene_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

if (ncol(tx2gene) < 2) stop("tx2gene.tsv has <2 columns, unexpected.")
tx2gene <- tx2gene[, 1:2]
colnames(tx2gene) <- c("TXNAME", "GENEID")
tx2gene$TXNAME <- sub("\\..*$", "", tx2gene$TXNAME)
tx2gene <- tx2gene[!duplicated(tx2gene$TXNAME), , drop = FALSE]

## -----------------------------
## tximport (gene-level)
## -----------------------------
txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE,
                ignoreAfterBar = TRUE,
                countsFromAbundance = "no")

samples <- names(files)
group <- sub("_[0-9]+$", "", samples)
meta <- data.frame(sample = samples, group = factor(group), row.names = samples)

if ("WT_F_9M" %in% levels(meta$group)) meta$group <- relevel(meta$group, "WT_F_9M")

## -----------------------------
## Helper: gene labels (prefer gene_name from gene_tpm_counts.xlsx)
## -----------------------------

# Build Ensembl->SYMBOL map (optional fallback)
gene_ids <- rownames(txi$counts)
ens_ids  <- sub("\\..*$", "", gene_ids)

sym_map <- mapIds(
  org.Mm.eg.db,
  keys = unique(ens_ids),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

# Read gene_id -> gene_name from your Excel (preferred)
gene_name_map <- NULL
if (file.exists(GENE_TPM_XLSX)) {

  read_gene_map <- function(path) {
    if (requireNamespace("readxl", quietly = TRUE)) {
      df <- readxl::read_xlsx(path, .name_repair = "minimal")
    } else if (requireNamespace("openxlsx", quietly = TRUE)) {
      df <- openxlsx::read.xlsx(path)
    } else {
      stop("Need either 'readxl' or 'openxlsx' installed to read: ", path)
    }

    if (!all(c("gene_id", "gene_name") %in% colnames(df))) {
      stop("gene_tpm_counts file must contain columns: gene_id, gene_name\nFound: ",
           paste(colnames(df), collapse = ", "))
    }

    df <- df[, c("gene_id", "gene_name")]
    df$gene_id <- as.character(df$gene_id)
    df$gene_name <- as.character(df$gene_name)

    # strip Ensembl version so it matches your DESeq2 ids
    df$ens_novers <- sub("\\..*$", "", df$gene_id)

    # keep first mapping if duplicates exist
    df <- df[!duplicated(df$ens_novers) & !is.na(df$ens_novers) & nzchar(df$ens_novers), ]

    setNames(df$gene_name, df$ens_novers)
  }

  gene_name_map <- read_gene_map(GENE_TPM_XLSX)
  message("Loaded gene_name map from: ", GENE_TPM_XLSX,
          " (n=", length(gene_name_map), ")")
} else {
  warning("GENE_TPM_XLSX not found: ", GENE_TPM_XLSX, " — labels will use org.Mm.eg.db fallback only.")
}

make_labels <- function(ids) {
  ens <- sub("\\..*$", "", ids)

  # 1) Prefer gene_name from your Excel
  nm <- if (!is.null(gene_name_map)) unname(gene_name_map[ens]) else NA_character_

  # 2) Fallback to org.Mm.eg.db SYMBOL
  sym <- unname(sym_map[ens])

  # 3) Fallback to Ensembl ID
  out <- ifelse(!is.na(nm) & nzchar(nm), nm,
                ifelse(!is.na(sym) & nzchar(sym), sym, ens))
  out
}

## -----------------------------
## PCA helper
## -----------------------------
`%||%` <- function(a, b) if (!is.null(a) && nzchar(a)) a else b

make_pca_bubbles <- function(vsd_obj, intgroup_col = "group", title = NULL) {
  pca_df <- plotPCA(vsd_obj, intgroup = intgroup_col, returnData = TRUE)
  percentVar <- attr(pca_df, "percentVar")

  if (!("name" %in% colnames(pca_df))) pca_df$name <- rownames(pca_df)
  if (!(intgroup_col %in% colnames(pca_df))) stop("PCA data missing column: ", intgroup_col)

  pca_df[[intgroup_col]] <- factor(pca_df[[intgroup_col]])
  groups <- levels(pca_df[[intgroup_col]])

  hull_list <- list()
  for (g in groups) {
    df_g <- pca_df[pca_df[[intgroup_col]] == g, , drop = FALSE]
    hull_idx <- if (nrow(df_g) <= 2) seq_len(nrow(df_g)) else chull(df_g$PC1, df_g$PC2)
    hull_list[[g]] <- data.frame(PC1=df_g$PC1[hull_idx], PC2=df_g$PC2[hull_idx], group=g)
  }
  hulls <- do.call(rbind, hull_list)
  hulls$group <- factor(hulls$group, levels = levels(pca_df[[intgroup_col]]))

  ggplot() +
    geom_polygon(data=hulls, aes(x=PC1, y=PC2, group=group, fill=group), alpha=0.2, colour=NA) +
    geom_point(data=pca_df, aes(x=PC1, y=PC2, colour=.data[[intgroup_col]]), size=3) +
    geom_text_repel(
      data=pca_df,
      aes(x=PC1, y=PC2, label=name, colour=.data[[intgroup_col]]),
      size=2, box.padding=0.4, point.padding=0.6, max.overlaps=Inf,
      segment.alpha=0.6, show.legend=FALSE
    ) +
    xlab(paste0("PC1 (", round(percentVar[1] * 100, 1), "%)")) +
    ylab(paste0("PC2 (", round(percentVar[2] * 100, 1), "%)")) +
    ggtitle(title %||% "") +
    theme_bw() +
    theme(panel.grid=element_blank(), legend.title=element_blank(),
          plot.title=element_text(hjust=0.5))
}

## -----------------------------
## Overall PCA (all samples)
## -----------------------------
dds_all <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ group)
dds_all <- dds_all[rowSums(counts(dds_all)) > MIN_TOTAL_COUNT, ]
dds_all <- estimateSizeFactors(dds_all, type = "poscounts")
vsd_all <- vst(dds_all, blind = FALSE)

p_all <- make_pca_bubbles(vsd_all, intgroup_col = "group", title = "All samples")
ggsave(file.path(OUT_BASE, "PCA_all_samples.png"), p_all, width = 6, height = 5, dpi = 600)
ggsave(file.path(OUT_BASE, "PCA_all_samples.pdf"), p_all, width = 6, height = 5)

## -----------------------------
## Volcano helper
## -----------------------------
make_volcano <- function(res_df, out_png, title) {
  res_df$padj_plot <- res_df$padj
  res_df$padj_plot[is.na(res_df$padj_plot)] <- 1
  res_df$padj_plot <- pmax(res_df$padj_plot, .Machine$double.xmin)

  png(out_png, width = 3600, height = 3000, res = 300)
  print(EnhancedVolcano(
    res_df, lab=res_df$label, x="log2FoldChange", y="padj_plot",
    pCutoff=ALPHA, FCcutoff=LFC_CUTOFF, title=title
  ))
  dev.off()
}

## -----------------------------
## Define contrasts (g2 vs g1)
## -----------------------------
comparisons <- list(
  "WT_F_9M_vs_EM2_F_9M"   = c("WT_F_9M",   "EM2_F_9M"),
  "WT_F_9M_vs_EM3_F_9M"   = c("WT_F_9M",   "EM3_F_9M"),
  "WT_F_9M_vs_EM11_F_9M"  = c("WT_F_9M",   "EM11_F_9M"),
  "EM2_F_9M_vs_EM11_F_9M" = c("EM2_F_9M",  "EM11_F_9M"),
  "EM2_F_9M_vs_EM3_F_9M"  = c("EM2_F_9M",  "EM3_F_9M"),
  "EM11_F_9M_vs_EM3_F_9M" = c("EM11_F_9M", "EM3_F_9M")
)

summary_rows <- list()

## -----------------------------
## Run DESeq2 per comparison
## -----------------------------
for (nm in names(comparisons)) {

  g1 <- comparisons[[nm]][1]
  g2 <- comparisons[[nm]][2]

  if (!all(c(g1, g2) %in% levels(meta$group))) {
    warning("Skipping ", nm, " because one/both groups not present in meta.")
    next
  }

  message("\n=== ", nm, " (", g2, " vs ", g1, ") ===")

  keep_samp <- meta$group %in% c(g1, g2)
  meta_sub  <- droplevels(meta[keep_samp, , drop = FALSE])

  txi_sub <- txi
  txi_sub$counts    <- txi$counts[, rownames(meta_sub), drop = FALSE]
  txi_sub$abundance <- txi$abundance[, rownames(meta_sub), drop = FALSE]
  txi_sub$length    <- txi$length[, rownames(meta_sub), drop = FALSE]

  dds <- DESeqDataSetFromTximport(txi_sub, colData = meta_sub, design = ~ group)
  dds <- dds[rowSums(counts(dds)) > MIN_TOTAL_COUNT, ]
  dds <- DESeq(dds, sfType = "poscounts")

  # unshrunken stats (padj etc.) in the requested direction
  res <- results(dds, contrast = c("group", g2, g1), alpha = ALPHA)

  # shrink LFC safely even if DESeq2 stored the coef in the reverse direction
  if (requireNamespace("apeglm", quietly = TRUE)) {
    rn <- resultsNames(dds)

    coef_forward <- paste0("group_", g2, "_vs_", g1)
    coef_reverse <- paste0("group_", g1, "_vs_", g2)

    if (coef_forward %in% rn) {
      res_shr <- lfcShrink(dds, coef = coef_forward, type = "apeglm", res = res)
    } else if (coef_reverse %in% rn) {
      # shrink reverse coef, then flip sign to match g2 vs g1
      res_tmp <- lfcShrink(dds, coef = coef_reverse, type = "apeglm")
      res_shr <- res
      res_shr$log2FoldChange <- -res_tmp$log2FoldChange
      if (!is.null(res_tmp$lfcSE)) res_shr$lfcSE <- res_tmp$lfcSE
    } else {
      stop("Could not find coef name for ", g2, " vs ", g1,
           "\nTried: ", coef_forward, ", ", coef_reverse,
           "\nAvailable resultsNames(dds):\n", paste(rn, collapse = "\n"))
    }
  } else {
    warning("apeglm not installed; using unshrunken log2FoldChange for ", nm)
    res_shr <- res
  }

  res_df <- as.data.frame(res_shr)
  res_df$gene_id <- rownames(res_df)
  res_df$label   <- make_labels(res_df$gene_id)
  res_df <- res_df[order(res_df$padj), ]

  n_padj <- sum(res_df$padj < ALPHA, na.rm = TRUE)
  n_both <- sum(res_df$padj < ALPHA & abs(res_df$log2FoldChange) >= LFC_CUTOFF, na.rm = TRUE)

  message("padj < ", ALPHA, " : ", n_padj)
  message("padj < ", ALPHA, " AND |log2FC| >= ", LFC_CUTOFF, " : ", n_both)

  out_dir <- file.path(OUT_BASE, nm)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  norm_counts <- counts(dds, normalized = TRUE)
  write.csv(as.data.frame(norm_counts), file = file.path(out_dir, paste0(nm, "_normalized_counts.csv")))

  write.csv(res_df, file = file.path(out_dir, paste0(nm, "_DE_results.csv")), row.names = FALSE)

  sig_df <- res_df[!is.na(res_df$padj) &
                     (res_df$padj < ALPHA) &
                     (abs(res_df$log2FoldChange) >= LFC_CUTOFF), , drop = FALSE]
  write.csv(sig_df, file = file.path(out_dir, paste0(nm, "_DE_sig_padj_lfc.csv")), row.names = FALSE)

  up_df   <- sig_df[sig_df$log2FoldChange >=  LFC_CUTOFF, , drop = FALSE]
  down_df <- sig_df[sig_df$log2FoldChange <= -LFC_CUTOFF, , drop = FALSE]
  write.csv(up_df,   file = file.path(out_dir, paste0(nm, "_DE_sig_UP_padj_lfc.csv")), row.names = FALSE)
  write.csv(down_df, file = file.path(out_dir, paste0(nm, "_DE_sig_DOWN_padj_lfc.csv")), row.names = FALSE)

  vsd <- vst(dds, blind = TRUE)
  p_pca <- make_pca_bubbles(vsd, intgroup_col = "group", title = nm)
  ggsave(file.path(out_dir, paste0(nm, "_PCA.png")), p_pca, width = 6, height = 5, dpi = 600)
  ggsave(file.path(out_dir, paste0(nm, "_PCA.pdf")), p_pca, width = 6, height = 5)

  make_volcano(res_df, file.path(out_dir, paste0(nm, "_volcano.png")), title = nm)

  res_no_gdpd1 <- res_df[!(res_df$label %in% c("Gdpd1", "GDPD1")) &
                           !(res_df$gene_id %in% c("Gdpd1", "GDPD1")), , drop = FALSE]
  make_volcano(res_no_gdpd1,
               file.path(out_dir, paste0(nm, "_volcano_noGdpd1.png")),
               title = paste0(nm, " (no Gdpd1)"))

  summary_rows[[nm]] <- data.frame(
    comparison = nm, g1 = g1, g2 = g2,
    genes_tested = nrow(res_df),
    n_padj_lt_alpha = n_padj,
    n_padj_and_lfc = n_both,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(OUT_BASE, "summary_counts.csv"), row.names = FALSE)

message("\nDone. Outputs in: ", OUT_BASE)