# ==============================================================================
# RNA Harmonization AI Project - Environment Setup
# ==============================================================================
# This script installs all required R packages for the project
# Run this ONCE at the beginning of the project
# ==============================================================================

cat("Starting package installation...\n")
cat("This may take 15-30 minutes depending on your system.\n\n")

# ==============================================================================
# 1. CORE DATA MANIPULATION PACKAGES
# ==============================================================================
cat("Installing core packages...\n")
core_packages <- c(
  "tidyverse",      # Data manipulation and visualization
  "data.table",     # Fast data handling
  "readr",          # Reading data files
  "here"            # Project-relative paths
)

install.packages(core_packages)

# ==============================================================================
# 2. BIOCONDUCTOR SETUP
# ==============================================================================
cat("\nSetting up Bioconductor...\n")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# ==============================================================================
# 3. BIOCONDUCTOR PACKAGES - DATA ACCESS
# ==============================================================================
cat("\nInstalling Bioconductor data access packages...\n")
BiocManager::install(c(
  "TCGAbiolinks",          # TCGA data download and processing
  "GEOquery",              # GEO data access
  "SummarizedExperiment"   # Data structure for genomic data
), update = FALSE, ask = FALSE)

# ==============================================================================
# 4. BIOCONDUCTOR PACKAGES - NORMALIZATION & PREPROCESSING
# ==============================================================================
cat("\nInstalling normalization packages...\n")
BiocManager::install(c(
  "DESeq2",       # RNA-seq normalization (VST, rlog)
  "edgeR",        # Alternative normalization (TMM, CPM)
  "limma"         # Additional preprocessing tools
), update = FALSE, ask = FALSE)

# ==============================================================================
# 5. BATCH CORRECTION METHODS
# ==============================================================================
cat("\nInstalling batch correction packages...\n")
BiocManager::install(c(
  "sva",          # ComBat and ComBat-seq
  "harmony",      # Harmony batch correction
  "batchelor"     # MNN and other methods
), update = FALSE, ask = FALSE)

# ==============================================================================
# 6. MACHINE LEARNING PACKAGES
# ==============================================================================
cat("\nInstalling machine learning packages...\n")
ml_packages <- c(
  "randomForest",   # Random Forest classifier
  "xgboost",        # XGBoost classifier
  "caret",          # ML workflow and cross-validation
  "pROC",           # ROC curves and AUC
  "glmnet"          # Elastic net (optional)
)

install.packages(ml_packages)

# ==============================================================================
# 7. VISUALIZATION PACKAGES
# ==============================================================================
cat("\nInstalling visualization packages...\n")
viz_packages <- c(
  "ggplot2",        # Core plotting (included in tidyverse)
  "pheatmap",       # Heatmaps
  "ggrepel",        # Better label placement
  "RColorBrewer",   # Color palettes
  "viridis",        # Perceptually uniform colors
  "patchwork"       # Combine plots
)

install.packages(viz_packages)

BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE)

# ==============================================================================
# 8. DIMENSIONALITY REDUCTION
# ==============================================================================
cat("\nInstalling dimensionality reduction packages...\n")
install.packages(c("umap", "Rtsne"))

# ==============================================================================
# 9. SHINY APP PACKAGES
# ==============================================================================
cat("\nInstalling Shiny packages...\n")
shiny_packages <- c(
  "shiny",
  "shinydashboard",
  "shinyWidgets",
  "plotly",         # Interactive plots
  "DT"              # Interactive tables
)

install.packages(shiny_packages)

# ==============================================================================
# 10. UTILITY PACKAGES
# ==============================================================================
cat("\nInstalling utility packages...\n")
utility_packages <- c(
  "parallel",       # Parallel processing
  "doParallel",     # Parallel backend
  "foreach",        # Parallel loops
  "tictoc"          # Timing code execution
)

install.packages(utility_packages)

# ==============================================================================
# VERIFICATION
# ==============================================================================
cat("\n\n==============================================================================\n")
cat("INSTALLATION COMPLETE - VERIFYING PACKAGES\n")
cat("==============================================================================\n\n")

# List of critical packages to verify
critical_packages <- c(
  "tidyverse", "TCGAbiolinks", "GEOquery", "DESeq2", 
  "sva", "randomForest", "xgboost", "caret", "shiny"
)

all_loaded <- TRUE
for (pkg in critical_packages) {
  if (require(pkg, quietly = TRUE, character.only = TRUE)) {
    cat("✓", pkg, "successfully installed\n")
  } else {
    cat("✗", pkg, "FAILED to install\n")
    all_loaded <- FALSE
  }
}

cat("\n")
if (all_loaded) {
  cat("SUCCESS! All critical packages are ready.\n")
  cat("You can proceed to data download.\n")
} else {
  cat("WARNING: Some packages failed to install.\n")
  cat("Please review error messages above and reinstall failed packages.\n")
}

cat("\n==============================================================================\n")
cat("Next steps:\n")
cat("1. Run scripts/01_download_tcga_data.R to download TCGA data\n")
cat("2. Run scripts/02_download_geo_data.R to download GEO data\n")
cat("==============================================================================\n")