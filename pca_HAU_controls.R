# PCA for HAU controls from nf-core/rnaseq star_salmon output
# with semi-transparent bubbles and non-overlapping labels

library(ggplot2)
library(ggrepel)

## -------- Paths --------
base_dir   <- "/media/pontikos_nas2/HarisQurashi/projects/9_Owen_BulkRNAseq/outputs/1_HAU_controls"
counts_tsv <- file.path(base_dir, "star_salmon", "salmon.merged.gene_counts_length_scaled.tsv")
meta_csv   <- "/media/pontikos_nas2/HarisQurashi/projects/9_Owen_BulkRNAseq/data/samplesheet.csv"
out_dir    <- file.path(base_dir, "z1_PCA")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## -------- Load data --------
counts_raw <- read.delim(counts_tsv, check.names = FALSE)

gene_ids <- counts_raw[[1]]
colnames(counts_raw)[1] <- "gene_id"

meta <- read.csv(meta_csv, stringsAsFactors = FALSE)

sample_cols <- intersect(colnames(counts_raw), meta$sample)
if (length(sample_cols) == 0) stop("No overlapping sample names between counts file and samplesheet.")

counts <- counts_raw[, sample_cols, drop = FALSE]
rownames(counts) <- gene_ids

meta <- meta[match(sample_cols, meta$sample), , drop = FALSE]
rownames(meta) <- meta$sample

## -------- Filter + transform --------
keep <- rowSums(counts) > 0
counts <- counts[keep, ]

log_counts <- log2(counts + 1)

## -------- PCA --------
p <- prcomp(t(log_counts), center = TRUE, scale. = TRUE)
var_expl <- p$sdev^2 / sum(p$sdev^2)

pca_df <- data.frame(
  sample = rownames(p$x),
  PC1    = p$x[, 1],
  PC2    = p$x[, 2],
  stringsAsFactors = FALSE
)

## -------- Top loadings for PC1 / PC2 --------
loadings <- p$rotation   # rows = genes, cols = PCs

get_top_loadings <- function(pc, n = 30) {
  vals <- loadings[, pc]
  ord  <- order(abs(vals), decreasing = TRUE)
  data.frame(
    gene        = rownames(loadings)[ord][1:n],
    loading     = vals[ord][1:n],
    abs_loading = abs(vals[ord][1:n]),
    stringsAsFactors = FALSE
  )
}

top_pc1 <- get_top_loadings("PC1", n = 30)
top_pc2 <- get_top_loadings("PC2", n = 30)

write.table(top_pc1,
            file = file.path(out_dir, "PC1_top30_loadings.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(top_pc2,
            file = file.path(out_dir, "PC2_top30_loadings.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
            
cat("Top PC1 loadings:\n")
print(head(top_pc1, 10))
cat("\nTop PC2 loadings:\n")
print(head(top_pc2, 10))

pca_df <- merge(pca_df, meta, by.x = "sample", by.y = "sample", all.x = TRUE)

## -------- Groups --------
group_map <- c(
  "HAU2001"      = "HAU200x",
  "HAU2002"      = "HAU200x",
  "HAU2003"      = "HAU200x",
  "HAU2901"      = "HAU290x",
  "HAU2902"      = "HAU290x",
  "HAU2903"      = "HAU290x",
  "HAU701_mix"   = "HAU70x",
  "HAU702_mix"   = "HAU70x",
  "HAU703_mix"   = "HAU70x",
  "HAU704"       = "HAU70x",
  "HAU1001_mix"  = "HAU100x",
  "HAU1002_mix"  = "HAU100x",
  "HAU1003_mix"  = "HAU100x",
  "HAU1004_mix"  = "HAU100x",
  "HAU1251"      = "HAU125x",
  "HAU1252_mix"  = "HAU125x",
  "HAU1253"      = "HAU125x",
  "HAU1254_mix"  = "HAU125x",
  "HAU1601_mix"  = "HAU160x",
  "HAU1602"      = "HAU160x",
  "HAU1603_mix"  = "HAU160x",
  "HAU1604_mix"  = "HAU160x"
)

pca_df$group <- group_map[pca_df$sample]
pca_df$group[is.na(pca_df$group)] <- "Other"

desired_levels <- c("HAU70x", "HAU100x", "HAU125x", "HAU160x", "HAU200x", "HAU290x", "Other")
pca_df$group <- factor(pca_df$group,
                       levels = desired_levels[desired_levels %in% unique(pca_df$group)])

## -------- Convex hulls --------
groups <- sort(unique(pca_df$group))
hull_list <- list()

for (g in groups) {
  df_g <- pca_df[pca_df$group == g, , drop = FALSE]
  if (nrow(df_g) == 1) {
    hull_idx <- 1
  } else if (nrow(df_g) == 2) {
    hull_idx <- 1:2
  } else {
    hull_idx <- chull(df_g$PC1, df_g$PC2)
  }
  hull_list[[g]] <- data.frame(
    PC1   = df_g$PC1[hull_idx],
    PC2   = df_g$PC2[hull_idx],
    group = g
  )
}

hulls <- do.call(rbind, hull_list)
hulls$group <- factor(hulls$group, levels = levels(pca_df$group))

## -------- Plot --------
p_pca <- ggplot() +
  # bubbles
  geom_polygon(
    data = hulls,
    aes(x = PC1, y = PC2, group = group, fill = group),
    alpha = 0.2,
    colour = NA
  ) +
  # points
  geom_point(
    data = pca_df,
    aes(x = PC1, y = PC2, colour = group),
    size = 3
  ) +
  # non-overlapping labels, above points
  geom_text_repel(
    data = pca_df,
    aes(x = PC1, y = PC2, label = sample, colour = group),
    size = 2,
    box.padding   = 0.4,
    point.padding = 0.6,
    max.overlaps  = Inf,
    segment.alpha = 0.6,
    show.legend   = FALSE
  ) +
  xlab(paste0("PC1 (", round(var_expl[1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(var_expl[2] * 100, 1), "%)")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank()
  )

ggsave(file.path(out_dir, "pca_groups_bubbles.pdf"), plot = p_pca, width = 6, height = 5)
ggsave(file.path(out_dir, "pca_groups_bubbles.png"), plot = p_pca, width = 6, height = 5, dpi = 600)