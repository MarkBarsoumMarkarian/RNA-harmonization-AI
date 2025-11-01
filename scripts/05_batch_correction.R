# ==============================================================================
# Batch Correction Comparison - ComBat, Harmony, and MNN
# ==============================================================================
# This script compares three batch correction methods and evaluates their
# effectiveness at removing technical variation while preserving biology
# ==============================================================================

library(tidyverse)
library(sva)          # ComBat
library(harmony)      # Harmony
library(batchelor)    # MNN
library(ggplot2)
library(patchwork)

cat("==============================================================================\n")
cat("BATCH CORRECTION COMPARISON\n")
cat("==============================================================================\n\n")

# ==============================================================================
# 1. LOAD PREPROCESSED DATA
# ==============================================================================
cat("Step 1: Loading preprocessed data...\n")

combined_expr <- readRDS("data/processed/combined_expression_raw.rds")
sample_metadata <- readRDS("data/processed/sample_metadata.rds")

cat("✓ Data loaded:\n")
cat("  Genes:", nrow(combined_expr), "\n")
cat("  Samples:", ncol(combined_expr), "\n")
cat("  Datasets:", length(unique(sample_metadata$dataset)), "\n\n")

# ==============================================================================
# 2. PCA ON UNCORRECTED DATA (BASELINE)
# ==============================================================================
cat("Step 2: Running PCA on uncorrected data...\n")

pca_uncorrected <- prcomp(t(combined_expr), scale. = FALSE)
var_exp_uncorrected <- (pca_uncorrected$sdev^2 / sum(pca_uncorrected$sdev^2)) * 100

pca_df_uncorrected <- data.frame(
  PC1 = pca_uncorrected$x[, 1],
  PC2 = pca_uncorrected$x[, 2],
  dataset = sample_metadata$dataset,
  platform = sample_metadata$platform
)

cat("✓ PCA complete\n")
cat("  PC1:", round(var_exp_uncorrected[1], 2), "%\n")
cat("  PC2:", round(var_exp_uncorrected[2], 2), "%\n\n")

# ==============================================================================
# 3. METHOD 1: COMBAT BATCH CORRECTION
# ==============================================================================
cat("Step 3: Applying ComBat batch correction...\n")

# ComBat requires batch as a factor
batch <- factor(sample_metadata$dataset)

# Run ComBat
expr_combat <- ComBat(
  dat = combined_expr,
  batch = batch,
  mod = NULL,  # No model matrix (unsupervised)
  par.prior = TRUE,
  prior.plots = FALSE
)

cat("✓ ComBat complete\n\n")

# PCA on ComBat-corrected data
pca_combat <- prcomp(t(expr_combat), scale. = FALSE)
var_exp_combat <- (pca_combat$sdev^2 / sum(pca_combat$sdev^2)) * 100

pca_df_combat <- data.frame(
  PC1 = pca_combat$x[, 1],
  PC2 = pca_combat$x[, 2],
  dataset = sample_metadata$dataset,
  platform = sample_metadata$platform
)

# ==============================================================================
# 4. METHOD 2: HARMONY BATCH CORRECTION
# ==============================================================================
cat("Step 4: Applying Harmony batch correction...\n")

# Harmony works on PCA space
pca_for_harmony <- prcomp(t(combined_expr), scale. = FALSE)

# Run Harmony
harmony_obj <- HarmonyMatrix(
  data_mat = pca_for_harmony$x,
  meta_data = sample_metadata,
  vars_use = "dataset",
  do_pca = FALSE,
  verbose = FALSE
)

# Harmony returns corrected PCA coordinates
pca_df_harmony <- data.frame(
  PC1 = harmony_obj[, 1],
  PC2 = harmony_obj[, 2],
  dataset = sample_metadata$dataset,
  platform = sample_metadata$platform
)

var_exp_harmony <- (apply(harmony_obj, 2, stats::var) / sum(apply(harmony_obj, 2, stats::var))) * 100

cat("✓ Harmony complete\n\n")

# ==============================================================================
# 5. METHOD 3: MNN (MUTUAL NEAREST NEIGHBORS)
# ==============================================================================
cat("Step 5: Applying MNN batch correction...\n")

# Split data by batch
batch_list <- lapply(unique(sample_metadata$dataset), function(b) {
  combined_expr[, sample_metadata$dataset == b]
})
names(batch_list) <- unique(sample_metadata$dataset)

# Run MNN
mnn_result <- do.call(fastMNN, c(batch_list, list(auto.merge = TRUE)))

# Extract corrected expression
expr_mnn <- assay(mnn_result, "reconstructed")
colnames(expr_mnn) <- colnames(combined_expr)

cat("✓ MNN complete\n\n")

# PCA on MNN-corrected data
pca_mnn <- prcomp(t(expr_mnn), scale. = FALSE)
var_exp_mnn <- (pca_mnn$sdev^2 / sum(pca_mnn$sdev^2)) * 100

pca_df_mnn <- data.frame(
  PC1 = pca_mnn$x[, 1],
  PC2 = pca_mnn$x[, 2],
  dataset = sample_metadata$dataset,
  platform = sample_metadata$platform
)

# ==============================================================================
# 6. VISUALIZE COMPARISONS
# ==============================================================================
cat("Step 6: Creating comparison visualizations...\n")

# Function to create PCA plot
plot_pca <- function(pca_df, var_exp, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = dataset)) +
    geom_point(size = 2, alpha = 0.7) +
    labs(
      title = title,
      x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2], 1), "%)")
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
}

# Create all plots
p1 <- plot_pca(pca_df_uncorrected, var_exp_uncorrected, "A. Uncorrected")
p2 <- plot_pca(pca_df_combat, var_exp_combat, "B. ComBat")
p3 <- plot_pca(pca_df_harmony, var_exp_harmony, "C. Harmony")
p4 <- plot_pca(pca_df_mnn, var_exp_mnn, "D. MNN")

# Combine plots
combined_plot <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Batch Correction Comparison",
    subtitle = "Effect of different methods on dataset integration"
  )

ggsave("figures/05_batch_correction_comparison.png", combined_plot,
       width = 14, height = 12, dpi = 300)

cat("✓ Saved: figures/05_batch_correction_comparison.png\n\n")

# ==============================================================================
# 7. QUANTITATIVE EVALUATION
# ==============================================================================
cat("Step 7: Quantitative evaluation of batch mixing...\n\n")

# Function to calculate silhouette score (batch vs. other)
calc_batch_silhouette <- function(pca_coords, batch_labels) {
  dist_matrix <- dist(pca_coords[, 1:10])  # Use first 10 PCs
  sil <- cluster::silhouette(as.integer(factor(batch_labels)), dist_matrix)
  mean(sil[, 3])
}

# Calculate metrics
metrics <- data.frame(
  Method = c("Uncorrected", "ComBat", "Harmony", "MNN"),
  Batch_Silhouette = c(
    calc_batch_silhouette(pca_uncorrected$x, sample_metadata$dataset),
    calc_batch_silhouette(pca_combat$x, sample_metadata$dataset),
    calc_batch_silhouette(harmony_obj, sample_metadata$dataset),
    calc_batch_silhouette(pca_mnn$x, sample_metadata$dataset)
  )
)

cat("Batch Mixing Metrics (lower is better):\n")
print(metrics)
cat("\nInterpretation: Lower silhouette scores indicate better batch mixing\n\n")

# ==============================================================================
# 8. SAVE CORRECTED DATA
# ==============================================================================
cat("Step 8: Saving batch-corrected data...\n")

# Save all corrected versions
saveRDS(expr_combat, "data/processed/combined_expression_combat.rds")
saveRDS(expr_mnn, "data/processed/combined_expression_mnn.rds")
saveRDS(harmony_obj, "data/processed/combined_pca_harmony.rds")

# Save evaluation metrics
write.csv(metrics, "results/batch_correction_metrics.csv", row.names = FALSE)

# Save PCA results for all methods
pca_results <- list(
  uncorrected = pca_df_uncorrected,
  combat = pca_df_combat,
  harmony = pca_df_harmony,
  mnn = pca_df_mnn
)
saveRDS(pca_results, "results/pca_all_methods.rds")

cat("✓ All corrected data saved\n\n")

# ==============================================================================
# 9. SELECT BEST METHOD FOR DOWNSTREAM ANALYSIS
# ==============================================================================
cat("Step 9: Selecting best method...\n")

# Based on lowest batch silhouette
best_method_idx <- which.min(metrics$Batch_Silhouette[-1]) + 1
best_method <- metrics$Method[best_method_idx]

cat("Recommended method:", best_method, "\n")
cat("(Lowest batch silhouette score)\n\n")

# Save best method for next steps
if (best_method == "ComBat") {
  final_expr <- expr_combat
} else if (best_method == "MNN") {
  final_expr <- expr_mnn
} else {
  # For Harmony, we need to use ComBat as fallback since Harmony returns PCA coords
  cat("Note: Harmony returns PCA coordinates. Using ComBat for ML pipeline.\n")
  final_expr <- expr_combat
}

saveRDS(final_expr, "data/processed/combined_expression_corrected.rds")
cat("✓ Saved final corrected expression matrix\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("==============================================================================\n")
cat("BATCH CORRECTION COMPLETE\n")
cat("==============================================================================\n\n")

cat("METHODS COMPARED:\n")
cat("  1. ComBat (empirical Bayes)\n")
cat("  2. Harmony (PCA-based)\n")
cat("  3. MNN (mutual nearest neighbors)\n\n")

cat("RESULTS:\n")
print(metrics)
cat("\nBEST METHOD:", best_method, "\n\n")

cat("FILES CREATED:\n")
cat("  - data/processed/combined_expression_combat.rds\n")
cat("  - data/processed/combined_expression_mnn.rds\n")
cat("  - data/processed/combined_pca_harmony.rds\n")
cat("  - data/processed/combined_expression_corrected.rds (best method)\n")
cat("  - figures/05_batch_correction_comparison.png\n")
cat("  - results/batch_correction_metrics.csv\n\n")

cat("NEXT STEP:\n")
cat("  Run scripts/06_ml_models.R to train classifiers\n")
cat("==============================================================================\n")