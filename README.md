# Batch-Harmonized AI for Multi-Center Pancreatic Cancer RNA Data

**A reproducible pipeline for cross-cohort RNA biomarker discovery using batch correction and machine learning**

---

## Project Overview

This project develops an AI-powered framework to harmonize RNA expression data across multiple cancer datasets (TCGA + GEO) and build robust machine learning classifiers for pancreatic cancer subtyping. The goal is to create a reproducible, clinically-translatable pipeline that addresses batch effects, a major barrier in multi-center biomarker studies.

### Why Pancreatic Cancer?
- **High unmet need**: 5-year survival rate <10%
- **Understudied**: Less AI/ML work compared to breast/lung cancer
- **Clear subtypes**: Classical, Basal-like, with distinct prognosis
- **Novel contribution**: Limited prior work on cross-cohort harmonization

---

## Research Objectives

1. **Harmonize RNA data** across TCGA and GEO datasets using multiple batch correction methods
2. **Build ML classifiers** (RandomForest, XGBoost) for cancer subtype prediction
3. **Validate cross-cohort performance** with true external validation
4. **Deploy interactive tool** via Shiny app for reproducible analysis

---

## Project Structure

```
rna-harmonization-ai/
├── data/
│   ├── raw/              # Original TCGA & GEO downloads
│   └── processed/        # Cleaned, normalized, harmonized data
├── scripts/
│   ├── 00_setup_environment.R      # Install all packages
│   ├── 01_download_tcga_data.R     # Download TCGA-PAAD
│   ├── 02_download_geo_data.R      # Download GEO validation sets
│   ├── 03_eda_initial.R            # Exploratory data analysis
│   ├── 04_preprocessing.R          # Normalization & feature alignment
│   ├── 05_batch_correction.R       # ComBat, Harmony, MNN comparison
│   ├── 06_ml_models.R              # Train ML classifiers
│   └── 07_evaluation.R             # Cross-cohort validation
├── figures/              # All generated plots
├── results/              # Model outputs, performance metrics
├── shiny_app/            # Interactive web application
└── docs/                 # Manuscript drafts, supplementary materials
```

---

## Quick Start

### 1. Setup Environment (First Time Only)
```r
# In RStudio, run:
source("scripts/00_setup_environment.R")
```
This installs all required packages (~20-30 minutes).

### 2. Download Data
```r
# Download TCGA pancreatic cancer data (~2-3 GB, 15-30 min)
source("scripts/01_download_tcga_data.R")

# Download GEO validation datasets (~5-10 min)
source("scripts/02_download_geo_data.R")
```

### 3. Run Initial EDA
```r
# Explore data and identify batch effects
source("scripts/03_eda_initial.R")

# Check generated figures in figures/ folder
```

---

## Datasets

### Primary Dataset: TCGA-PAAD
- **Project**: The Cancer Genome Atlas - Pancreatic Adenocarcinoma
- **Samples**: 178 primary tumors
- **Data type**: RNA-seq (STAR aligned, raw counts)
- **Platform**: Illumina HiSeq
- **Clinical data**: Vital status, stage, demographics

### Validation Dataset: GEO
- **GSE71729**: RNA-seq PDAC cohort (357 samples, Illumina)

**Note**: We focused on RNA-seq datasets only to maximize gene overlap (14,137 common genes vs. <100 with microarray inclusion). This provides better feature coverage for robust machine learning models.

---

## Methods Pipeline

### 1. Preprocessing
- Normalize RNA counts (DESeq2 VST or log2 CPM)
- Filter low-count genes
- Align gene features across platforms
- Handle missing values

### 2. Batch Correction (Comparison Study)
- **ComBat** (empirical Bayes)
- **Harmony** (PCA-based)
- **MNN** (mutual nearest neighbors)
- Evaluate: PCA mixing, silhouette scores, kBET

### 3. Machine Learning
- **RandomForest**: Interpretable, robust to outliers
- **XGBoost**: High performance gradient boosting
- **Neural Network** (optional): Deep learning approach

### 4. Evaluation Strategy
- **Training**: TCGA-PAAD (in-distribution)
- **Testing**: GEO datasets (out-of-distribution)
- **Metrics**: Accuracy, F1-score, AUROC
- **Interpretability**: Feature importance, SHAP values

---

## Expected Outcomes

### Key Results Demonstrated
1. **Gene coverage**: 14,137 common genes between TCGA and GEO RNA-seq
2. **Batch correction**: ComBat successfully harmonizes datasets (PCA visualization)
3. **Cross-cohort classification**: Models trained on TCGA (178 samples) predict on GEO (357 samples)
4. **Stable biomarkers**: 5 top genes (LAMC2, DKK1, ITGB6, GPRC5A, MAL2) with biological relevance
5. **Survival prediction**: Alive vs Dead classification (85 vs 93 samples, balanced classes)

### Figures for Manuscript
- PCA before/after batch correction
- Cross-cohort classification performance
- Feature importance and biological interpretation
- Shiny app screenshots

---

## Technology Stack

- **Language**: R 4.0+
- **Data Access**: TCGAbiolinks, GEOquery
- **Preprocessing**: DESeq2, edgeR, limma
- **Batch Correction**: sva (ComBat), harmony, batchelor (MNN)
- **Machine Learning**: randomForest, xgboost, caret
- **Visualization**: ggplot2, pheatmap, ComplexHeatmap
- **Deployment**: Shiny, shinydashboard

---

## Manuscript Roadmap

### Target Journal
- **bioRxiv** (preprint) → **Bioinformatics** or **Briefings in Bioinformatics**

### Estimated Timeline
- **Weeks 1-2**: Data download & EDA
- **Weeks 3-4**: Batch correction comparison
- **Weeks 5-6**: ML model development
- **Weeks 7-8**: Cross-cohort validation
- **Weeks 9-10**: Shiny app development
- **Weeks 11-12**: Manuscript writing

### Novelty Statement
*"We present the first open-source R pipeline that integrates batch correction (ComBat) with machine learning for robust, cross-cohort RNA biomarker discovery in pancreatic cancer. By focusing on RNA-seq platforms, we achieved high gene coverage (14,137 genes) enabling identification of 5 novel prognostic biomarkers. Our Shiny application enables reproducible analysis and accelerates clinical translation."*

---

## Contributing

This is a research project. For questions or collaborations, please open an issue or contact the author.

---

## Key References

1. Johnson et al. (2007) - ComBat batch correction
2. Korsunsky et al. (2019) - Harmony algorithm
3. Haghverdi et al. (2018) - MNN method
4. TCGA Research Network - Pancreatic cancer characterization

---

## License

MIT License - Open source for academic and commercial use.

---

**Project Status**: Active Development  
**Last Updated**: October 2025  
**Author**: Mark Barsoum Markarian
