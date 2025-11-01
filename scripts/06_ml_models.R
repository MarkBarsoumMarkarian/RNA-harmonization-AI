# ==============================================================================
# Machine Learning Models for Cross-Cohort Classification
# ==============================================================================
# This script:
#   1. Trains RandomForest and XGBoost classifiers on TCGA data
#   2. Tests on GEO validation cohorts (true external validation)
#   3. Evaluates with and without batch correction
#   4. Extracts feature importance
# ==============================================================================

library(tidyverse)
library(randomForest)
library(xgboost)
library(caret)
library(pROC)

cat("==============================================================================\n")
cat("MACHINE LEARNING PIPELINE\n")
cat("==============================================================================\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
cat("Step 1: Loading data...\n")

# Load both uncorrected and corrected data
expr_raw <- readRDS("data/processed/combined_expression_raw.rds")
expr_corrected <- readRDS("data/processed/combined_expression_corrected.rds")
sample_metadata <- readRDS("data/processed/sample_metadata.rds")

# Load TCGA clinical data for labels
tcga_clinical <- readRDS("data/processed/tcga_paad_clinical_filtered.rds")

cat("✓ Data loaded\n")
cat("  Genes:", nrow(expr_raw), "\n")
cat("  Samples:", ncol(expr_raw), "\n\n")

# ==============================================================================
# 2. PREPARE LABELS FOR TCGA DATA  
# ==============================================================================
cat("Step 2: Preparing classification labels...\n")

# Get TCGA samples from expression data
tcga_samples <- colnames(expr_raw)[sample_metadata$dataset == "TCGA"]

cat("TCGA samples found:", length(tcga_samples), "\n")

# Extract vital status (most reliable label)
if ("vital_status" %in% colnames(tcga_clinical)) {
labels <- tcga_clinical[tcga_samples, "vital_status"]
names(labels) <- tcga_samples

# IMPORTANT: as.character() removes names, so save them first
label_names <- names(labels)
labels <- as.character(labels)
names(labels) <- label_names
  
  # Remove NAs
  labels <- labels[!is.na(labels)]
  
  cat("✓ Using vital status classification\n")
  cat("  Total samples with labels:", length(labels), "\n")
  cat("  Classes:", paste(unique(labels), collapse = " vs "), "\n")
  cat("  Distribution:\n")
  print(table(labels))
  cat("\n")
  
  label_type <- "vital_status"
  
  # Convert to factor AFTER getting labeled_samples
  # Keep the names for now!
  
} else {
  stop("vital_status not found in clinical data")
}

# Verify we have enough samples
if (length(labels) < 20) {
  stop("Insufficient samples with labels: ", length(labels))
}

# Verify we have 2 classes
if (length(unique(labels)) < 2) {
  stop("Need at least 2 classes for classification")
}

cat("Final label summary:\n")
cat("  Samples:", length(labels), "\n")
cat("  Classes:", paste(unique(labels), collapse = " vs "), "\n")
cat("  Class sizes:", paste(table(labels), collapse = " / "), "\n\n")

# ==============================================================================
# 3. FEATURE SELECTION (TOP VARIABLE GENES)
# ==============================================================================
cat("Step 3: Selecting top variable features...\n")

# Get samples with labels
labeled_samples <- names(labels)

cat("  Labeled samples:", length(labeled_samples), "\n")

# Subset expression data to labeled samples only
cat("  Subsetting expression matrices...\n")
expr_raw_labeled <- expr_raw[, labeled_samples, drop = FALSE]
expr_corrected_labeled <- expr_corrected[, labeled_samples, drop = FALSE]

cat("  Raw expression:", nrow(expr_raw_labeled), "x", ncol(expr_raw_labeled), "\n")
cat("  Corrected expression:", nrow(expr_corrected_labeled), "x", ncol(expr_corrected_labeled), "\n")

# Calculate variance for each gene
cat("  Calculating gene variance...\n")
gene_vars <- apply(expr_corrected_labeled, 1, stats::var)

# Select top variable genes (or all if fewer than 500)
n_features <- min(500, length(gene_vars))
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:n_features])

cat("✓ Selected", n_features, "most variable genes\n")
cat("  Top genes:", paste(head(top_genes, 3), collapse = ", "), "\n\n")

# ==============================================================================
# 4. PREPARE TRAINING DATA
# ==============================================================================
cat("Step 4: Preparing training data...\n")

# Subset to top genes first
cat("  Subsetting to top genes...\n")
expr_raw_top <- expr_raw_labeled[top_genes, , drop = FALSE]
expr_corrected_top <- expr_corrected_labeled[top_genes, , drop = FALSE]

cat("  Creating training matrices...\n")
# Training: Use only the labeled TCGA samples and selected genes
X_train_raw <- t(expr_raw_top)
X_train_corrected <- t(expr_corrected_top)

# NOW convert labels to factor (after we've used the names)
y_train <- factor(labels[rownames(X_train_raw)])

cat("✓ Training data prepared:\n")
cat("  Samples:", nrow(X_train_raw), "\n")
cat("  Features:", ncol(X_train_raw), "\n")
cat("  Classes:", paste(levels(y_train), "=", table(y_train), collapse = " / "), "\n\n")

# ==============================================================================
# 5. PREPARE TEST DATA (GEO SAMPLES)
# ==============================================================================
cat("Step 5: Preparing test data (GEO validation)...\n")

# Test: All non-TCGA samples
geo_mask <- !colnames(expr_raw) %in% labeled_samples
geo_samples <- colnames(expr_raw)[geo_mask]

if (length(geo_samples) > 0) {
  X_test_raw <- t(expr_raw[top_genes, geo_samples])
  X_test_corrected <- t(expr_corrected[top_genes, geo_samples])
  test_datasets <- sample_metadata[geo_samples, "dataset"]
  
  cat("✓ Test data prepared:\n")
  cat("  Samples:", nrow(X_test_raw), "\n")
  cat("  Datasets:", paste(unique(test_datasets), collapse = ", "), "\n\n")
} else {
  cat("  No GEO samples available for testing\n\n")
  X_test_raw <- NULL
  X_test_corrected <- NULL
  test_datasets <- NULL
}

# ==============================================================================
# 6. MODEL 1: RANDOM FOREST
# ==============================================================================
cat("Step 6: Training Random Forest models...\n\n")

set.seed(42)

# RF on uncorrected data
cat("  Training RF on uncorrected data...\n")
rf_raw <- randomForest(
  x = X_train_raw,
  y = y_train,
  ntree = 500,
  importance = TRUE
)

# RF on corrected data
cat("  Training RF on corrected data...\n")
rf_corrected <- randomForest(
  x = X_train_corrected,
  y = y_train,
  ntree = 500,
  importance = TRUE
)

cat("✓ Random Forest models trained\n\n")

# ==============================================================================
# 7. MODEL 2: XGBOOST
# ==============================================================================
cat("Step 7: Training XGBoost models...\n\n")

# Convert labels to numeric (0/1)
y_train_numeric <- as.numeric(y_train) - 1

# XGBoost parameters
xgb_params <- list(
  objective = "binary:logistic",
  eval_metric = "logloss",
  max_depth = 6,
  eta = 0.3,
  nrounds = 100
)

# XGBoost on uncorrected data
cat("  Training XGBoost on uncorrected data...\n")
dtrain_raw <- xgb.DMatrix(data = X_train_raw, label = y_train_numeric)
xgb_raw <- xgb.train(
  params = xgb_params,
  data = dtrain_raw,
  nrounds = 100,
  verbose = 0
)

# XGBoost on corrected data
cat("  Training XGBoost on corrected data...\n")
dtrain_corrected <- xgb.DMatrix(data = X_train_corrected, label = y_train_numeric)
xgb_corrected <- xgb.train(
  params = xgb_params,
  data = dtrain_corrected,
  nrounds = 100,
  verbose = 0
)

cat("✓ XGBoost models trained\n\n")

# ==============================================================================
# 7. EVALUATION ON TEST SET (GEO COHORTS)
# ==============================================================================
cat("Step 7: Evaluating on GEO validation cohorts...\n\n")

# Note: We can't evaluate without test labels
# For now, we'll demonstrate the prediction workflow
# In a real scenario, you'd need to map GEO samples to outcomes

cat("NOTE: External validation requires labels for GEO samples.\n")
cat("For demonstration, showing prediction probabilities:\n\n")

# Predictions on test set (uncorrected)
pred_rf_raw <- predict(rf_raw, X_test_raw, type = "prob")
pred_xgb_raw <- predict(xgb_raw, xgb.DMatrix(X_test_raw))

# Predictions on test set (corrected)
pred_rf_corrected <- predict(rf_corrected, X_test_corrected, type = "prob")
pred_xgb_corrected <- predict(xgb_corrected, xgb.DMatrix(X_test_corrected))

# Show distribution of predictions by dataset
results_summary <- data.frame(
  dataset = test_datasets,
  rf_raw_prob = pred_rf_raw[, 2],
  rf_corrected_prob = pred_rf_corrected[, 2],
  xgb_raw_prob = pred_xgb_raw,
  xgb_corrected_prob = pred_xgb_corrected
)

cat("Prediction summary (mean probability per dataset):\n")
print(results_summary %>%
  group_by(dataset) %>%
  summarise(
    RF_raw = mean(rf_raw_prob),
    RF_corrected = mean(rf_corrected_prob),
    XGB_raw = mean(xgb_raw_prob),
    XGB_corrected = mean(xgb_corrected_prob)
  ))
cat("\n")

# ==============================================================================
# 9. FEATURE IMPORTANCE
# ==============================================================================
cat("Step 9: Extracting feature importance...\n")

# RF feature importance
rf_importance <- importance(rf_corrected)
rf_top_features <- head(rf_importance[order(-rf_importance[, "MeanDecreaseGini"]), ], 20)

# XGBoost feature importance
xgb_importance <- xgb.importance(model = xgb_corrected)
xgb_top_features <- head(xgb_importance, 20)

# Plot feature importance
png("figures/06_feature_importance_rf.png", width = 10, height = 8, units = "in", res = 300)
varImpPlot(rf_corrected, n.var = 20, main = "Random Forest - Top 20 Features")
dev.off()

png("figures/06_feature_importance_xgb.png", width = 10, height = 8, units = "in", res = 300)
xgb.plot.importance(xgb_top_features, main = "XGBoost - Top 20 Features")
dev.off()

cat("✓ Feature importance plots saved\n\n")

# ==============================================================================
# 10. SAVE MODELS AND RESULTS
# ==============================================================================
cat("Step 10: Saving models and results...\n")

# Save models
saveRDS(rf_raw, "results/model_rf_raw.rds")
saveRDS(rf_corrected, "results/model_rf_corrected.rds")
saveRDS(xgb_raw, "results/model_xgb_raw.rds")
saveRDS(xgb_corrected, "results/model_xgb_corrected.rds")

# Save feature importance
write.csv(rf_top_features, "results/rf_feature_importance.csv")
write.csv(xgb_top_features, "results/xgb_feature_importance.csv")

# Save predictions
if (!is.null(results_summary)) {
  saveRDS(results_summary, "results/test_predictions.rds")
}

cat("✓ All models and results saved\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("==============================================================================\n")
cat("MACHINE LEARNING PIPELINE COMPLETE\n")
cat("==============================================================================\n\n")

cat("MODELS TRAINED:\n")
cat("  1. Random Forest (uncorrected)\n")
cat("  2. Random Forest (batch-corrected)\n")
cat("  3. XGBoost (uncorrected)\n")
cat("  4. XGBoost (batch-corrected)\n\n")

cat("CLASSIFICATION TASK:\n")
cat("  Type:", label_type, "\n")
cat("  Classes:", paste(levels(labels), collapse = " vs "), "\n")
cat("  Training samples:", length(y_train), "\n")
cat("  Test samples:", nrow(X_test_raw), "\n\n")

cat("TOP 5 IMPORTANT GENES (Random Forest):\n")
print(head(rownames(rf_top_features), 5))
cat("\n")

cat("FILES CREATED:\n")
cat("  - results/model_rf_raw.rds\n")
cat("  - results/model_rf_corrected.rds\n")
cat("  - results/model_xgb_raw.rds\n")
cat("  - results/model_xgb_corrected.rds\n")
cat("  - results/rf_feature_importance.csv\n")
cat("  - results/xgb_feature_importance.csv\n")
cat("  - figures/06_feature_importance_rf.png\n")
cat("  - figures/06_feature_importance_xgb.png\n\n")

cat("NEXT STEP:\n")
cat("  Review feature importance and prepare manuscript figures\n")
cat("  Consider: scripts/07_shiny_app.R for interactive tool\n")
cat("==============================================================================\n")