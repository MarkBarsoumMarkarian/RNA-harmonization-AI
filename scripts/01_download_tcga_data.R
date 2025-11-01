# ==============================================================================
# Download TCGA Pancreatic Cancer (PAAD) Data
# ==============================================================================
# This script downloads RNA-seq data from TCGA-PAAD project
# Data includes: gene expression counts, clinical data, and sample metadata
# ==============================================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

cat("==============================================================================\n")
cat("DOWNLOADING TCGA-PAAD (PANCREATIC ADENOCARCINOMA) DATA\n")
cat("==============================================================================\n\n")

# ==============================================================================
# 1. DEFINE QUERY PARAMETERS
# ==============================================================================
cat("Step 1: Setting up query parameters...\n")

project <- "TCGA-PAAD"
data_category <- "Transcriptome Profiling"
data_type <- "Gene Expression Quantification"
workflow_type <- "STAR - Counts"  # STAR aligner, raw counts

# ==============================================================================
# 2. QUERY TCGA DATABASE
# ==============================================================================
cat("Step 2: Querying TCGA database...\n")
cat("This may take a few minutes...\n\n")

query_rnaseq <- GDCquery(
  project = project,
  data.category = data_category,
  data.type = data_type,
  workflow.type = workflow_type,
  experimental.strategy = "RNA-Seq"
)

cat("Query Summary:\n")
cat("- Project:", project, "\n")
cat("- Samples found:", length(getResults(query_rnaseq)$cases), "\n")
cat("- Data type:", data_type, "\n\n")

# ==============================================================================
# 3. DOWNLOAD DATA
# ==============================================================================
cat("Step 3: Downloading data from GDC...\n")
cat("WARNING: This will download ~2-3 GB of data.\n")
cat("Download time: 10-30 minutes depending on internet speed.\n\n")

tryCatch({
  GDCdownload(
    query = query_rnaseq,
    method = "api",
    files.per.chunk = 10
  )
  cat("✓ Download completed successfully!\n\n")
}, error = function(e) {
  cat("✗ Download failed with error:\n")
  cat(as.character(e), "\n")
  cat("Try running this script again or check your internet connection.\n")
  stop(e)
})

# ==============================================================================
# 4. PREPARE DATA (PROCESS INTO R OBJECT)
# ==============================================================================
cat("Step 4: Processing downloaded data into SummarizedExperiment object...\n")
cat("This may take 5-10 minutes...\n\n")

tcga_paad_data <- GDCprepare(
  query = query_rnaseq,
  save = FALSE,
  summarizedExperiment = TRUE
)

cat("✓ Data processing complete!\n\n")

# ==============================================================================
# 5. EXPLORE BASIC DATA STRUCTURE
# ==============================================================================
cat("Step 5: Exploring data structure...\n\n")

cat("Data Dimensions:\n")
cat("- Genes:", nrow(tcga_paad_data), "\n")
cat("- Samples:", ncol(tcga_paad_data), "\n\n")

cat("Sample Types:\n")
sample_types <- table(tcga_paad_data$sample_type)
print(sample_types)
cat("\n")

cat("Available Clinical Variables:\n")
clinical_vars <- colnames(colData(tcga_paad_data))
cat("Total variables:", length(clinical_vars), "\n")
cat("Key variables:\n")
key_vars <- c("sample_type", "gender", "age_at_diagnosis", "vital_status", 
              "tumor_stage", "paper_Subtype", "ajcc_pathologic_stage")
available_key_vars <- key_vars[key_vars %in% clinical_vars]
for (var in available_key_vars) {
  cat("  -", var, "\n")
}
cat("\n")

# ==============================================================================
# 6. SAVE PROCESSED DATA
# ==============================================================================
cat("Step 6: Saving processed data...\n")

# Save full object
output_file <- file.path("data", "raw", "tcga_paad_rnaseq.rds")
saveRDS(tcga_paad_data, file = output_file)
cat("✓ Saved full data to:", output_file, "\n")

# Extract and save clinical data separately for easy access
clinical_data <- as.data.frame(colData(tcga_paad_data))
clinical_file <- file.path("data", "raw", "tcga_paad_clinical.csv")
write.csv(clinical_data, file = clinical_file, row.names = TRUE)
cat("✓ Saved clinical data to:", clinical_file, "\n\n")

# ==============================================================================
# 7. GENERATE SUMMARY REPORT
# ==============================================================================
cat("==============================================================================\n")
cat("DOWNLOAD SUMMARY\n")
cat("==============================================================================\n")
cat("Downloaded and processed TCGA-PAAD RNA-seq data\n\n")

cat("DATASETS CREATED:\n")
cat("1. tcga_paad_rnaseq.rds - Full SummarizedExperiment object\n")
cat("2. tcga_paad_clinical.csv - Clinical data in table format\n\n")

cat("DATA STRUCTURE:\n")
cat("- Genes (rows):", nrow(tcga_paad_data), "\n")
cat("- Samples (columns):", ncol(tcga_paad_data), "\n")
cat("- Primary Solid Tumors:", sum(tcga_paad_data$sample_type == "Primary Solid Tumor"), "\n")
cat("- Solid Tissue Normal:", sum(tcga_paad_data$sample_type == "Solid Tissue Normal"), "\n\n")

cat("NEXT STEPS:\n")
cat("1. Run 02_download_geo_data.R to get validation cohorts\n")
cat("2. Run 03_eda.R to explore the data\n")
cat("==============================================================================\n")

# Clean up
cat("\nCleaning up temporary files...\n")
# GDC creates temp folders - optional cleanup
# unlink("GDCdata", recursive = TRUE)

cat("\n✓ Script completed successfully!\n")