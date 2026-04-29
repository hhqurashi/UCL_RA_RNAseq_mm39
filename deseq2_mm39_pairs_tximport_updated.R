#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(EnhancedVolcano)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

## -----------------------------
## USER EDITS
## -----------------------------
STAR_SALMON_DIR <- "/mnt/scratch/hqurashi/10_Sara_RNA_Seq/outputs/star_salmon"
OUT_BASE        <- "/mnt/scratch/hqurashi/10_Sara_RNA_Seq/outputs"
GENE_TPM_XLSX   <- "/mnt/scratch/hqurashi/10_Sara_RNA_Seq/outputs/star_salmon/gene_tpm_counts.xlsx"
SAMP_PAT        <- "^(WT|EM2|EM3|EM11)_F_9M_[0-9]+$"

MIN_TOTAL_COUNT <- 10
ALPHA           <- 0.05
LFC_CUTOFF      <- 1

# NEW: genes to force-label + bold in the extra volcano plot
HIGHLIGHT_GENES <- c("GDPD1","Gdpd1","YPEL2","Ypel2")

# (Optional) override the x tick step manually; otherwise auto-picked from data
FORCE_X_TICK_STEP <- NA_real_   # e.g. 2, 5, 10; leave NA to auto

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

quant_files  <- vapply(sample_dirs, pick_quant, character(1))
sample_names <- basename(sample_dirs)

ok <- !is.na(quant_files) & nzchar(quant_files)
quant_files  <- quant_files[ok]
sample_names <- sample_names[ok]

keep <- grepl(SAMP_PAT, sample_names)
quant_files  <- quant_files[keep]
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
  hits <- list.files(STAR_SALMON_DIR, pattern="^tx2gene\\.tsv$", recursive=TRUE, full.names=TRUE)
  if (length(hits) == 0) stop("Could not find tx2gene.tsv under: ", STAR_SALMON_DIR)
  tx2gene_path <- hits[1]
}
message("Using tx2gene: ", tx2gene_path)

tx2gene <- read.delim(tx2gene_path, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
if (ncol(tx2gene) < 2) stop("tx2gene.tsv has <2 columns, unexpected.")
tx2gene <- tx2gene[, 1:2]
colnames(tx2gene) <- c("TXNAME", "GENEID")
tx2gene$TXNAME <- sub("\\..*$", "", tx2gene$TXNAME)
tx2gene <- tx2gene[!duplicated(tx2gene$TXNAME), , drop=FALSE]

## -----------------------------
## tximport (gene-level)
## -----------------------------
txi <- tximport(
  files,
  type="salmon",
  tx2gene=tx2gene,
  ignoreTxVersion=TRUE,
  ignoreAfterBar=TRUE,
  countsFromAbundance="no"
)

samples <- names(files)
group <- sub("_[0-9]+$", "", samples)
meta <- data.frame(sample=samples, group=factor(group), row.names=samples)

if ("WT_F_9M" %in% levels(meta$group)) meta$group <- relevel(meta$group, "WT_F_9M")

## -----------------------------
## Helper: gene labels
## -----------------------------
gene_ids <- rownames(txi$counts)
ens_ids  <- sub("\\..*$", "", gene_ids)

sym_map <- mapIds(
  org.Mm.eg.db,
  keys = unique(ens_ids),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

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

    if (!all(c("gene_id","gene_name") %in% colnames(df))) {
      stop("gene_tpm_counts file must contain columns: gene_id, gene_name\nFound: ",
           paste(colnames(df), collapse=", "))
    }

    df <- df[, c("gene_id","gene_name")]
    df$gene_id   <- as.character(df$gene_id)
    df$gene_name <- as.character(df$gene_name)
    df$ens_novers <- sub("\\..*$", "", df$gene_id)

    df <- df[!duplicated(df$ens_novers) & !is.na(df$ens_novers) & nzchar(df$ens_novers), ]
    setNames(df$gene_name, df$ens_novers)
  }

  gene_name_map <- read_gene_map(GENE_TPM_XLSX)
  message("Loaded gene_name map from: ", GENE_TPM_XLSX, " (n=", length(gene_name_map), ")")
} else {
  warning("GENE_TPM_XLSX not found: ", GENE_TPM_XLSX, " — labels will use org.Mm.eg.db fallback only.")
}

make_labels <- function(ids) {
  ens <- sub("\\..*$", "", ids)
  nm  <- if (!is.null(gene_name_map)) unname(gene_name_map[ens]) else NA_character_
  sym <- unname(sym_map[ens])

  out <- ifelse(!is.na(nm) & nzchar(nm), nm,
                ifelse(!is.na(sym) & nzchar(sym), sym, ens))
  out
}

## -----------------------------
## PCA helper
## -----------------------------
`%||%` <- function(a, b) if (!is.null(a) && nzchar(a)) a else b

make_pca_bubbles <- function(vsd_obj, intgroup_col="group", title=NULL) {
  pca_df <- plotPCA(vsd_obj, intgroup=intgroup_col, returnData=TRUE)
  percentVar <- attr(pca_df, "percentVar")

  if (!("name" %in% colnames(pca_df))) pca_df$name <- rownames(pca_df)
  if (!(intgroup_col %in% colnames(pca_df))) stop("PCA data missing column: ", intgroup_col)

  pca_df[[intgroup_col]] <- factor(pca_df[[intgroup_col]])
  groups <- levels(pca_df[[intgroup_col]])

  hull_list <- list()
  for (g in groups) {
    df_g <- pca_df[pca_df[[intgroup_col]] == g, , drop=FALSE]
    hull_idx <- if (nrow(df_g) <= 2) seq_len(nrow(df_g)) else chull(df_g$PC1, df_g$PC2)
    hull_list[[g]] <- data.frame(PC1=df_g$PC1[hull_idx], PC2=df_g$PC2[hull_idx], group=g)
  }
  hulls <- do.call(rbind, hull_list)
  hulls$group <- factor(hulls$group, levels=levels(pca_df[[intgroup_col]]))

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
dds_all <- DESeqDataSetFromTximport(txi, colData=meta, design=~ group)
dds_all <- dds_all[rowSums(counts(dds_all)) > MIN_TOTAL_COUNT, ]
dds_all <- estimateSizeFactors(dds_all, type="poscounts")
vsd_all <- vst(dds_all, blind=FALSE)

p_all <- make_pca_bubbles(vsd_all, intgroup_col="group", title="All samples")
ggsave(file.path(OUT_BASE, "PCA_all_samples.png"), p_all, width=6, height=5, dpi=600)
ggsave(file.path(OUT_BASE, "PCA_all_samples.pdf"), p_all, width=6, height=5)

## -----------------------------
## NEW: PCA driver genes (PC1/PC2 loadings) + PCA excluding selected samples
## -----------------------------
EXCLUDE_SAMPLES <- c("EM3_F_9M_5", "EM3_F_9M_1")  # run both "leave-one-out" PCAs
PCA_NTOP_VAR   <- 500    # match plotPCA default behaviour
PCA_TOP_GENES  <- 50     # how many driver genes to report per PC

# row variance helper (use matrixStats if available)
rowVars_fast <- function(m) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::rowVars(m, na.rm = TRUE))
  } else {
    return(apply(m, 1, var, na.rm = TRUE))
  }
}

# returns a single table with top drivers for PC1 + PC2
get_pc_driver_table <- function(vsd_obj, ntop = 500, top_n = 50) {
  mat <- assay(vsd_obj)

  rv <- rowVars_fast(mat)
  ntop <- min(ntop, length(rv))
  sel <- order(rv, decreasing = TRUE)[seq_len(ntop)]

  # PCA on samples x genes, so rotation is genes x PCs
  pca <- prcomp(t(mat[sel, , drop = FALSE]), center = TRUE, scale. = FALSE)
  rot <- pca$rotation

  mk <- function(pc_idx, pc_name) {
    v <- rot[, pc_idx]
    df <- data.frame(
      gene_id = rownames(rot),
      loading = as.numeric(v),
      abs_loading = abs(as.numeric(v)),
      stringsAsFactors = FALSE
    )
    df <- df[order(-df$abs_loading), , drop = FALSE]
    df <- head(df, top_n)
    df$label <- make_labels(df$gene_id)
    df$PC <- pc_name
    df
  }

  rbind(mk(1, "PC1"), mk(2, "PC2"))
}

write_pc_drivers <- function(vsd_obj, tag) {
  drivers <- get_pc_driver_table(vsd_obj, ntop = PCA_NTOP_VAR, top_n = PCA_TOP_GENES)
  out_csv <- file.path(
    OUT_BASE,
    paste0("PCA_driver_genes_", tag, "_top", PCA_TOP_GENES, "_ntop", PCA_NTOP_VAR, ".csv")
  )
  write.csv(drivers, out_csv, row.names = FALSE)
  message("Wrote PCA driver genes: ", out_csv)
}

run_pca_excluding <- function(dds_obj, drop_samples) {
  drop_samples <- unique(drop_samples)

  missing <- setdiff(drop_samples, colnames(dds_obj))
  if (length(missing) > 0) {
    warning("Some samples not found for exclusion: ", paste(missing, collapse = ", "))
  }

  keep_cols <- setdiff(colnames(dds_obj), drop_samples)
  if (length(keep_cols) < 3) {
    warning("Too few samples left after exclusion; skipping PCA.")
    return(invisible(NULL))
  }

  dds_no <- dds_obj[, keep_cols]
  dds_no$group <- droplevels(dds_no$group)

  dds_no <- estimateSizeFactors(dds_no, type = "poscounts")
  vsd_no <- vst(dds_no, blind = FALSE)

  tag <- paste0("no_", paste(drop_samples, collapse = "_"))

  p_no <- make_pca_bubbles(
    vsd_no, intgroup_col = "group",
    title = paste0("All samples (excluding ", paste(drop_samples, collapse = " + "), ")")
  )

  ggsave(file.path(OUT_BASE, paste0("PCA_all_samples_", tag, ".png")),
         p_no, width = 6, height = 5, dpi = 600)
  ggsave(file.path(OUT_BASE, paste0("PCA_all_samples_", tag, ".pdf")),
         p_no, width = 6, height = 5)

  write_pc_drivers(vsd_no, tag = tag)
  invisible(NULL)
}

# Driver genes for all-samples PCA
write_pc_drivers(vsd_all, tag = "all_samples")

# Leave-one-out PCAs (your existing behaviour)
for (s in EXCLUDE_SAMPLES) run_pca_excluding(dds_all, drop_samples = s)

# NEW: exclude BOTH together
run_pca_excluding(dds_all, drop_samples = EXCLUDE_SAMPLES)

## -----------------------------
## Volcano helpers (NEW: consistent x-axis + bold forced labels)
## -----------------------------
pick_tick_step <- function(max_abs) {
  if (is.na(max_abs) || max_abs <= 0) return(1)
  if (!is.na(FORCE_X_TICK_STEP)) return(FORCE_X_TICK_STEP)
  if (max_abs <= 4)  return(1)
  if (max_abs <= 10) return(2)
  if (max_abs <= 25) return(5)
  return(10)
}

make_xscale <- function(max_abs) {
  step <- pick_tick_step(max_abs)
  lim  <- ceiling(max_abs / step) * step
  lim  <- max(lim, step)
  breaks <- seq(-lim, lim, by=step)
  list(xlim=c(-lim, lim), breaks=breaks, step=step, lim=lim)
}

prep_padj_plot <- function(df) {
  df$padj_plot <- df$padj
  df$padj_plot[is.na(df$padj_plot)] <- 1
  df$padj_plot <- pmax(df$padj_plot, .Machine$double.xmin)
  df
}

save_volcano <- function(p, out_png) {
  png(out_png, width=3600, height=3000, res=300)
  print(p)
  dev.off()
}

make_volcano_standard <- function(res_df, title, xlim_common, xbreaks_common, lab_vec=NULL) {
  df <- prep_padj_plot(res_df)
  if (is.null(lab_vec)) lab_vec <- df$label

  p <- EnhancedVolcano(
    df,
    lab = lab_vec,
    x = "log2FoldChange",
    y = "padj_plot",
    pCutoff = ALPHA,
    FCcutoff = LFC_CUTOFF,
    title = title
  ) +
    scale_x_continuous(limits=xlim_common, breaks=xbreaks_common) +
    theme(plot.title = element_text(hjust = 0.5))

  p
}

make_volcano_highlight <- function(res_df, title, xlim_common, xbreaks_common, highlight_genes) {
  df <- prep_padj_plot(res_df)

  # Base plot labels are normal BUT we hide the highlight genes so we can re-add them in bold once.
  base_lab <- df$label
  base_lab[tolower(base_lab) %in% tolower(highlight_genes)] <- ""

  p <- make_volcano_standard(df, title, xlim_common, xbreaks_common, lab_vec=base_lab)

  # Add bold labels for GDPD1 / YPEL2 regardless of significance
  hl <- df[tolower(df$label) %in% tolower(highlight_genes), , drop=FALSE]
  if (nrow(hl) > 0) {
    hl$yval <- -log10(hl$padj_plot)

    p <- p +
      ggrepel::geom_text_repel(
        data = hl,
        aes(x = log2FoldChange, y = yval, label = label),
        inherit.aes = FALSE,
        fontface = "bold",
        size = 4,
        box.padding = 0.6,
        point.padding = 0.6,
        max.overlaps = Inf,
        segment.alpha = 0.6
      )
  }

  p
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

# NEW: keep results in memory so we can compute a global x-axis, then plot consistently
res_list <- list()

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
  meta_sub  <- droplevels(meta[keep_samp, , drop=FALSE])

  txi_sub <- txi
  txi_sub$counts    <- txi$counts[, rownames(meta_sub), drop=FALSE]
  txi_sub$abundance <- txi$abundance[, rownames(meta_sub), drop=FALSE]
  txi_sub$length    <- txi$length[, rownames(meta_sub), drop=FALSE]

  dds <- DESeqDataSetFromTximport(txi_sub, colData=meta_sub, design=~ group)
  dds <- dds[rowSums(counts(dds)) > MIN_TOTAL_COUNT, ]
  dds <- DESeq(dds, sfType="poscounts")

  res <- results(dds, contrast=c("group", g2, g1), alpha=ALPHA)

  if (requireNamespace("apeglm", quietly=TRUE)) {
    rn <- resultsNames(dds)
    coef_forward <- paste0("group_", g2, "_vs_", g1)
    coef_reverse <- paste0("group_", g1, "_vs_", g2)

    if (coef_forward %in% rn) {
      res_shr <- lfcShrink(dds, coef=coef_forward, type="apeglm", res=res)
    } else if (coef_reverse %in% rn) {
      res_tmp <- lfcShrink(dds, coef=coef_reverse, type="apeglm")
      res_shr <- res
      res_shr$log2FoldChange <- -res_tmp$log2FoldChange
      if (!is.null(res_tmp$lfcSE)) res_shr$lfcSE <- res_tmp$lfcSE
    } else {
      stop("Could not find coef name for ", g2, " vs ", g1,
           "\nTried: ", coef_forward, ", ", coef_reverse,
           "\nAvailable resultsNames(dds):\n", paste(rn, collapse="\n"))
    }
  } else {
    warning("apeglm not installed; using unshrunken log2FoldChange for ", nm)
    res_shr <- res
  }

  res_df <- as.data.frame(res_shr)
  res_df$gene_id <- rownames(res_df)
  res_df$label   <- make_labels(res_df$gene_id)
  res_df <- res_df[order(res_df$padj), ]

  n_padj <- sum(res_df$padj < ALPHA, na.rm=TRUE)
  n_both <- sum(res_df$padj < ALPHA & abs(res_df$log2FoldChange) >= LFC_CUTOFF, na.rm=TRUE)

  message("padj < ", ALPHA, " : ", n_padj)
  message("padj < ", ALPHA, " AND |log2FC| >= ", LFC_CUTOFF, " : ", n_both)

  out_dir <- file.path(OUT_BASE, nm)
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

  norm_counts <- counts(dds, normalized=TRUE)
  write.csv(as.data.frame(norm_counts), file=file.path(out_dir, paste0(nm, "_normalized_counts.csv")))

  write.csv(res_df, file=file.path(out_dir, paste0(nm, "_DE_results.csv")), row.names=FALSE)

  sig_df <- res_df[!is.na(res_df$padj) &
                     (res_df$padj < ALPHA) &
                     (abs(res_df$log2FoldChange) >= LFC_CUTOFF), , drop=FALSE]
  write.csv(sig_df, file=file.path(out_dir, paste0(nm, "_DE_sig_padj_lfc.csv")), row.names=FALSE)

  up_df   <- sig_df[sig_df$log2FoldChange >=  LFC_CUTOFF, , drop=FALSE]
  down_df <- sig_df[sig_df$log2FoldChange <= -LFC_CUTOFF, , drop=FALSE]
  write.csv(up_df,   file=file.path(out_dir, paste0(nm, "_DE_sig_UP_padj_lfc.csv")), row.names=FALSE)
  write.csv(down_df, file=file.path(out_dir, paste0(nm, "_DE_sig_DOWN_padj_lfc.csv")), row.names=FALSE)

  vsd <- vst(dds, blind=TRUE)
  p_pca <- make_pca_bubbles(vsd, intgroup_col="group", title=nm)
  ggsave(file.path(out_dir, paste0(nm, "_PCA.png")), p_pca, width=6, height=5, dpi=600)
  ggsave(file.path(out_dir, paste0(nm, "_PCA.pdf")), p_pca, width=6, height=5)

  # NEW: store for consistent volcano plotting after we know global x-axis
  res_list[[nm]] <- res_df

  summary_rows[[nm]] <- data.frame(
    comparison = nm, g1 = g1, g2 = g2,
    genes_tested = nrow(res_df),
    n_padj_lt_alpha = n_padj,
    n_padj_and_lfc = n_both,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(OUT_BASE, "summary_counts.csv"), row.names=FALSE)

## -----------------------------
## NEW: Make ALL volcano plots with common x-axis limits + tick breaks
## -----------------------------
if (length(res_list) > 0) {
  max_abs_lfc <- max(vapply(res_list, function(df) max(abs(df$log2FoldChange), na.rm=TRUE), numeric(1)), na.rm=TRUE)
  xs <- make_xscale(max_abs_lfc)

  message("\nUsing common volcano x-axis: [", xs$xlim[1], ", ", xs$xlim[2], "] with tick step ", xs$step)

  for (nm in names(res_list)) {
    res_df <- res_list[[nm]]
    out_dir <- file.path(OUT_BASE, nm)

    # 1) Standard volcano (consistent x-axis)
    p1 <- make_volcano_standard(res_df, title=nm, xlim_common=xs$xlim, xbreaks_common=xs$breaks)
    save_volcano(p1, file.path(out_dir, paste0(nm, "_volcano.png")))

    # 2) Volcano excluding GDPD1 (consistent x-axis)
    res_no_gdpd1 <- res_df[!(tolower(res_df$label) %in% c("gdpd1")) &
                             !(tolower(res_df$gene_id) %in% c("gdpd1")), , drop=FALSE]
    p2 <- make_volcano_standard(res_no_gdpd1, title=paste0(nm, " (no Gdpd1)"),
                                xlim_common=xs$xlim, xbreaks_common=xs$breaks)
    save_volcano(p2, file.path(out_dir, paste0(nm, "_volcano_noGdpd1.png")))

    # 3) NEW: Standard volcano + forced bold labels for GDPD1 + YPEL2 (consistent x-axis)
    p3 <- make_volcano_highlight(res_df,
                                 title=paste0(nm, " (highlight GDPD1 + YPEL2)"),
                                 xlim_common=xs$xlim,
                                 xbreaks_common=xs$breaks,
                                 highlight_genes=HIGHLIGHT_GENES)
    save_volcano(p3, file.path(out_dir, paste0(nm, "_volcano_highlight_GDPD1_YPEL2.png")))
  }
}

message("\nDone. Outputs in: ", OUT_BASE)