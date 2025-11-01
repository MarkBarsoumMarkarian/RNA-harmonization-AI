# ==============================================================================
# Simplified EDA - TCGA PAAD (Troubleshooting Version)
# ==============================================================================

library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

cat("Starting EDA...\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
cat("Loading data...\n")
tcga_data <- readRDS("C:/Users/Bmark/Desktop/rna-harmonization-ai/data/raw/tcga_paad_rnaseq.rds")

# Extract components
counts <- assay(tcga_data, "unstranded")
clinical <- as.data.frame(colData(tcga_data))

cat("✓ Loaded:", nrow(counts), "genes x", ncol(counts), "samples\n\n")

# ==============================================================================
# 2. FILTER TO PRIMARY TUMORS
# ==============================================================================
cat("Filtering to primary tumors...\n")
tumor_idx <- which(clinical$sample_type == "Primary Tumor")
counts_tumor <- counts[, tumor_idx]
clinical_tumor <- clinical[tumor_idx, ]

cat("✓ Filtered:", ncol(counts_tumor), "tumor samples\n\n")

# ==============================================================================
# 3. FILTER LOW-COUNT GENES
# ==============================================================================
cat("Filtering genes...\n")
keep <- rowSums(counts_tumor >= 10) >= 10
counts_filtered <- counts_tumor[keep, ]

cat("✓ Kept", nrow(counts_filtered), "genes\n\n")

# ==============================================================================
# 4. VST NORMALIZATION
# ==============================================================================
cat("Normalizing with VST...\n")
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = clinical_tumor,
  design = ~ 1
)

vst_data <- vst(dds, blind = TRUE)
vst_matrix <- assay(vst_data)

cat("✓ VST complete\n\n")

# ==============================================================================
# 5. PCA
# ==============================================================================
cat("Running PCA...\n")
pca_result <- prcomp(t(vst_matrix), scale. = FALSE)
var_exp <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

cat("✓ PC1:", round(var_exp[1], 2), "%\n")
cat("✓ PC2:", round(var_exp[2], 2), "%\n\n")

# ==============================================================================
# 6. PLOT 1: PCA
# ==============================================================================
cat("Creating PCA plot...\n")
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)

p1 <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.7, color = "steelblue") +
  labs(
    title = "PCA - TCGA PAAD (BEFORE Batch Correction)",
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)")
  ) +
  theme_minimal(base_size = 12)

ggsave("figures/01_pca_before_correction.png", p1, 
       width = 8, height = 6, dpi = 300)
cat("✓ Saved: figures/01_pca_before_correction.png\n\n")

# ==============================================================================
# 7. PLOT 2: CORRELATION HEATMAP (SIMPLIFIED)
# ==============================================================================
cat("Creating correlation heatmap...\n")

# Use top 500 most variable genes - explicitly use stats::var
gene_vars <- apply(vst_matrix, 1, stats::var)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:500])
sample_cor <- stats::cor(vst_matrix[top_genes, ])

# Clean up row/column names to avoid any naming conflicts
rownames(sample_cor) <- paste0("S", 1:nrow(sample_cor))
colnames(sample_cor) <- paste0("S", 1:ncol(sample_cor))

# Simple heatmap without clustering annotations
png("figures/02_sample_correlation_heatmap.png", 
    width = 10, height = 10, units = "in", res = 300)

pheatmap(
  mat = sample_cor,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Sample Correlation Heatmap",
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 10,
  clustering_method = "complete"
)

dev.off()
cat("✓ Saved: figures/02_sample_correlation_heatmap.png\n\n")

# ==============================================================================
# 8. PLOT 3: LIBRARY SIZES
# ==============================================================================
cat("Creating library size plot...\n")

lib_sizes <- colSums(counts_tumor)

p3 <- ggplot(data.frame(size = lib_sizes), aes(x = size)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  labs(
    title = "Library Size Distribution",
    x = "Library Size (total counts)",
    y = "Number of Samples"
  ) +
  theme_minimal()

ggsave("figures/03_library_size_distribution.png", p3, 
       width = 8, height = 5, dpi = 300)
cat("✓ Saved: figures/03_library_size_distribution.png\n\n")

# ==============================================================================
# 9. SAVE PROCESSED DATA
# ==============================================================================
cat("Saving processed data...\n")

saveRDS(counts_filtered, "data/processed/tcga_paad_counts_filtered.rds")
saveRDS(vst_matrix, "data/processed/tcga_paad_vst.rds")
saveRDS(clinical_tumor, "data/processed/tcga_paad_clinical_filtered.rds")

cat("✓ All data saved\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("==============================================================================\n")
cat("EDA COMPLETE!\n")
cat("==============================================================================\n")
cat("Processed:\n")
cat("  -", nrow(counts_filtered), "genes\n")
cat("  -", ncol(counts_filtered), "samples\n")
cat("\nFigures created:\n")
cat("  1. PCA plot\n")
cat("  2. Sample correlation heatmap\n")
cat("  3. Library size distribution\n")
cat("\nNext: Download GEO validation data (script 02)\n")
cat("==============================================================================\n")