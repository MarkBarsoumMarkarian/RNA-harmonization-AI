# ==============================================================================
# Download GEO Pancreatic Cancer Validation Datasets
# ==============================================================================
# This script downloads independent pancreatic cancer datasets from GEO
# These will be used for external validation of the harmonization pipeline
# ==============================================================================

library(GEOquery)
library(tidyverse)

cat("==============================================================================\n")
cat("DOWNLOADING GEO PANCREATIC CANCER VALIDATION COHORTS\n")
cat("==============================================================================\n\n")

# ==============================================================================
# PANCREATIC CANCER DATASETS FROM GEO
# ==============================================================================
# Selected datasets with sufficient sample sizes and clinical annotations
# These are real pancreatic cancer cohorts from different centers

pancreatic_datasets <- list(
  # Dataset 1: GSE62452 - PDAC samples (Affymetrix)
  GSE62452 = list(
    description = "Pancreatic ductal adenocarcinoma (PDAC) cohort",
    platform = "Affymetrix",
    samples = "~130 samples"
  ),
  
  # Dataset 2: GSE28735 - Pancreatic cancer (Affymetrix)
  GSE28735 = list(
    description = "Pancreatic adenocarcinoma vs normal tissue",
    platform = "Affymetrix",
    samples = "~90 samples"
  ),
  
  # Dataset 3: GSE71729 - PDAC RNA-seq
  GSE71729 = list(
    description = "RNA-seq of pancreatic cancer (PDX models)",
    platform = "RNA-seq",
    samples = "~50 samples"
  )
)

cat("Target GEO Datasets:\n")
for (gse_id in names(pancreatic_datasets)) {
  info <- pancreatic_datasets[[gse_id]]
  cat(sprintf("- %s: %s (%s, %s)\n", 
              gse_id, info$description, info$platform, info$samples))
}
cat("\n")

# ==============================================================================
# FUNCTION: DOWNLOAD AND PROCESS GEO DATASET
# ==============================================================================
download_geo_dataset <- function(gse_id) {
  
  cat("==============================================================================\n")
  cat("Processing:", gse_id, "\n")
  cat("==============================================================================\n")
  
  tryCatch({
    # Download from GEO
    cat("Downloading", gse_id, "from GEO...\n")
    gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)
    
    # GEO can return multiple platform objects - take the first one
    if (length(gse) > 1) {
      cat("Multiple platforms detected. Using the first platform.\n")
    }
    gse <- gse[[1]]
    
    # Extract expression data
    expr_data <- exprs(gse)
    cat("✓ Expression matrix:", nrow(expr_data), "features x", 
        ncol(expr_data), "samples\n")
    
    # Extract phenotype/clinical data
    pheno_data <- pData(gse)
    cat("✓ Phenotype data:", nrow(pheno_data), "samples x", 
        ncol(pheno_data), "variables\n")
    
    # Extract feature data (gene annotations)
    feature_data <- fData(gse)
    cat("✓ Feature annotations:", nrow(feature_data), "features\n")
    
    # Create organized output
    output <- list(
      gse_id = gse_id,
      expression = expr_data,
      phenotype = pheno_data,
      features = feature_data,
      platform = annotation(gse)
    )
    
    # Save individual dataset
    output_file <- file.path("data", "raw", paste0("geo_", tolower(gse_id), ".rds"))
    saveRDS(output, file = output_file)
    cat("✓ Saved to:", output_file, "\n")
    
    # Save phenotype data as CSV for easy viewing
    csv_file <- file.path("data", "raw", paste0("geo_", tolower(gse_id), "_pheno.csv"))
    write.csv(pheno_data, file = csv_file, row.names = TRUE)
    cat("✓ Saved phenotype CSV to:", csv_file, "\n\n")
    
    return(list(success = TRUE, data = output))
    
  }, error = function(e) {
    cat("✗ Error downloading", gse_id, ":\n")
    cat(as.character(e), "\n\n")
    return(list(success = FALSE, error = e))
  })
}

# ==============================================================================
# DOWNLOAD ALL DATASETS
# ==============================================================================
cat("\nStarting downloads...\n")
cat("Note: Each dataset takes 2-5 minutes to download and process.\n\n")

results <- list()
for (gse_id in names(pancreatic_datasets)) {
  result <- download_geo_dataset(gse_id)
  results[[gse_id]] <- result
  Sys.sleep(2)  # Brief pause between downloads to be nice to GEO servers
}

# ==============================================================================
# SUMMARY REPORT
# ==============================================================================
cat("==============================================================================\n")
cat("DOWNLOAD SUMMARY\n")
cat("==============================================================================\n\n")

successful <- sum(sapply(results, function(x) x$success))
failed <- length(results) - successful

cat("Results:\n")
cat("- Successfully downloaded:", successful, "datasets\n")
cat("- Failed:", failed, "datasets\n\n")

if (successful > 0) {
  cat("Successfully downloaded datasets:\n")
  for (gse_id in names(results)) {
    if (results[[gse_id]]$success) {
      cat("✓", gse_id, "\n")
    }
  }
  cat("\n")
}

if (failed > 0) {
  cat("Failed datasets:\n")
  for (gse_id in names(results)) {
    if (!results[[gse_id]]$success) {
      cat("✗", gse_id, "\n")
    }
  }
  cat("\n")
}

cat("FILES CREATED:\n")
cat("Each dataset has:\n")
cat("- .rds file: Full data object (expression + phenotype + features)\n")
cat("- _pheno.csv file: Phenotype/clinical data in table format\n\n")

cat("NEXT STEPS:\n")
cat("1. Review phenotype CSV files to understand available clinical variables\n")
cat("2. Run 03_eda.R to explore data and identify batch effects\n")
cat("3. Identify common genes across TCGA and GEO datasets\n")
cat("==============================================================================\n")

cat("\n✓ GEO download script completed!\n")

# ==============================================================================
# BONUS: QUICK DATA INSPECTION
# ==============================================================================
cat("\n==============================================================================\n")
cat("QUICK DATA INSPECTION\n")
cat("==============================================================================\n\n")

for (gse_id in names(results)) {
  if (results[[gse_id]]$success) {
    data <- results[[gse_id]]$data
    cat(gse_id, ":\n")
    cat("  Platform:", data$platform, "\n")
    cat("  Samples:", ncol(data$expression), "\n")
    cat("  Features:", nrow(data$expression), "\n")
    cat("  Clinical variables:", ncol(data$phenotype), "\n\n")
  }
}