# ==============================================================================
# Simplified Preprocessing - RNA-seq Only (TCGA + GSE71729)
# ==============================================================================
# This version uses only RNA-seq datasets for better gene overlap
# ==============================================================================

library(tidyverse)
library(biomaRt)

cat("==============================================================================\n")
cat("PREPROCESSING - RNA-SEQ DATASETS ONLY\n")
cat("==============================================================================\n\n")

# ==============================================================================
# 1. LOAD TCGA DATA
# ==============================================================================
cat("Step 1: Loading TCGA data...\n")
tcga_vst <- readRDS("data/processed/tcga_paad_vst.rds")
tcga_clinical <- readRDS("data/processed/tcga_paad_clinical_filtered.rds")

cat("✓ TCGA:", nrow(tcga_vst), "genes x", ncol(tcga_vst), "samples\n\n")

# ==============================================================================
# 2. LOAD GSE71729 (RNA-SEQ)
# ==============================================================================
cat("Step 2: Loading GSE71729 (RNA-seq)...\n")
geo_71729 <- readRDS("data/raw/geo_gse71729.rds")

cat("✓ GSE71729:", nrow(geo_71729$expression), "features x", 
    ncol(geo_71729$expression), "samples\n\n")

# ==============================================================================
# 3. CONVERT TCGA ENSEMBL TO GENE SYMBOLS
# ==============================================================================
cat("Step 3: Converting TCGA Ensembl IDs to gene symbols...\n")

tcga_genes <- rownames(tcga_vst)
tcga_genes_clean <- sub("\\..*", "", tcga_genes)

tryCatch({
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  gene_mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = tcga_genes_clean,
    mart = ensembl
  )
  
  gene_mapping <- gene_mapping[gene_mapping$hgnc_symbol != "", ]
  
  tcga_idx <- match(gene_mapping$ensembl_gene_id, tcga_genes_clean)
  tcga_mapped <- tcga_vst[tcga_idx[!is.na(tcga_idx)], ]
  tcga_symbols <- gene_mapping$hgnc_symbol[!is.na(tcga_idx)]
  
  # Collapse duplicates
  tcga_collapsed <- do.call(rbind, lapply(unique(tcga_symbols), function(gene) {
    idx <- tcga_symbols == gene
    if (sum(idx) == 1) {
      tcga_mapped[idx, ]
    } else {
      colMeans(tcga_mapped[idx, , drop = FALSE])
    }
  }))
  rownames(tcga_collapsed) <- unique(tcga_symbols)
  
  cat("✓ Mapped", nrow(tcga_collapsed), "TCGA genes\n\n")
  
}, error = function(e) {
  cat("✗ BioMart failed. Using Ensembl IDs.\n\n")
  tcga_collapsed <- tcga_vst
  rownames(tcga_collapsed) <- tcga_genes_clean
})

# ==============================================================================
# 4. GET GSE71729 GENE SYMBOLS
# ==============================================================================
cat("Step 4: Extracting GSE71729 gene symbols...\n")

geo_genes <- rownames(geo_71729$expression)

cat("✓ GSE71729:", length(geo_genes), "genes\n")
cat("  Example genes:", paste(head(geo_genes, 5), collapse = ", "), "\n\n")

# ==============================================================================
# 5. FIND COMMON GENES
# ==============================================================================
cat("Step 5: Finding common genes...\n")

common_genes <- intersect(rownames(tcga_collapsed), geo_genes)

cat("✓ Common genes:", length(common_genes), "\n\n")

if (length(common_genes) < 100) {
  stop("Too few common genes. Check gene ID formats.")
}

# ==============================================================================
# 6. SUBSET TO COMMON GENES
# ==============================================================================
cat("Step 6: Subsetting to common genes...\n")

tcga_common <- tcga_collapsed[common_genes, ]
geo_common <- geo_71729$expression[common_genes, ]

# Check if GEO data needs log transformation
if (max(geo_common, na.rm = TRUE) > 100) {
  cat("  Log-transforming GSE71729...\n")
  geo_common <- log2(geo_common + 1)
}

cat("✓ Both datasets prepared\n\n")

# ==============================================================================
# 7. CREATE COMBINED MATRIX
# ==============================================================================
cat("Step 7: Creating combined expression matrix...\n")

combined_expr <- cbind(tcga_common, geo_common)

sample_metadata <- data.frame(
  sample_id = colnames(combined_expr),
  dataset = c(rep("TCGA", ncol(tcga_common)), 
              rep("GSE71729", ncol(geo_common))),
  platform = c(rep("RNA-seq_Illumina", ncol(tcga_common)),
               rep("RNA-seq_Illumina", ncol(geo_common))),
  row.names = colnames(combined_expr)
)

cat("✓ Combined matrix:\n")
cat("  Genes:", nrow(combined_expr), "\n")
cat("  Samples:", ncol(combined_expr), "\n")
cat("  TCGA:", ncol(tcga_common), "samples\n")
cat("  GSE71729:", ncol(geo_common), "samples\n\n")

# ==============================================================================
# 8. SAVE DATA
# ==============================================================================
cat("Step 8: Saving preprocessed data...\n")

saveRDS(combined_expr, "data/processed/combined_expression_raw.rds")
saveRDS(sample_metadata, "data/processed/sample_metadata.rds")
saveRDS(tcga_common, "data/processed/tcga_common_genes.rds")
saveRDS(list(GSE71729 = geo_common), "data/processed/geo_common_genes.rds")

cat("✓ All data saved\n\n")

# ==============================================================================
# 9. VISUALIZATION
# ==============================================================================
cat("Step 9: Creating visualization...\n")

set.seed(42)
sample_genes <- sample(rownames(combined_expr), min(1000, nrow(combined_expr)))

plot_data <- data.frame(
  value = as.vector(combined_expr[sample_genes, ]),
  dataset = rep(sample_metadata$dataset, each = length(sample_genes))
)

p <- ggplot(plot_data, aes(x = value, fill = dataset)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Expression Distribution - RNA-seq Datasets",
    subtitle = "Before Batch Correction",
    x = "Expression Value (log2)",
    y = "Density"
  ) +
  theme_minimal()

ggsave("figures/04_expression_distribution_rnaseq.png", p, 
       width = 10, height = 6, dpi = 300)

cat("✓ Saved visualization\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("==============================================================================\n")
cat("PREPROCESSING COMPLETE (RNA-SEQ ONLY)\n")
cat("==============================================================================\n\n")

cat("STRATEGY: Using only RNA-seq datasets for maximum gene overlap\n\n")

cat("DATASETS:\n")
cat("  TCGA-PAAD:", ncol(tcga_common), "samples\n")
cat("  GSE71729:", ncol(geo_common), "samples\n\n")

cat("FEATURES:\n")
cat("  Common genes:", nrow(combined_expr), "\n")
cat("  (", round(nrow(combined_expr)/nrow(tcga_vst)*100, 1), 
    "% of original TCGA genes)\n\n")

cat("FILES CREATED:\n")
cat("  - data/processed/combined_expression_raw.rds\n")
cat("  - data/processed/sample_metadata.rds\n")
cat("  - figures/04_expression_distribution_rnaseq.png\n\n")

cat("NEXT STEPS:\n")
cat("  1. Run scripts/05_batch_correction.R\n")
cat("  2. Run scripts/06_ml_models.R\n")
cat("==============================================================================\n")