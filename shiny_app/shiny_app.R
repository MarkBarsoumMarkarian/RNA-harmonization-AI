# ==============================================================================
# Interactive Shiny App for RNA Biomarker Discovery
# ==============================================================================
# This app allows users to:
#   - Explore batch correction effects
#   - View model predictions
#   - Examine gene expression patterns
#   - Download results
# ==============================================================================

library(shiny)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(DT)
library(randomForest)

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("Loading data for Shiny app...\n")

expr_raw <- readRDS("data/processed/combined_expression_raw.rds")
expr_corrected <- readRDS("data/processed/combined_expression_corrected.rds")
sample_metadata <- readRDS("data/processed/sample_metadata.rds")
tcga_clinical <- readRDS("data/processed/tcga_paad_clinical_filtered.rds")

rf_model <- readRDS("results/model_rf_corrected.rds")
rf_importance <- read.csv("results/rf_feature_importance.csv", row.names = 1)

# Get top genes
top_genes <- rownames(rf_importance[order(-rf_importance$MeanDecreaseGini), ])[1:50]

# PCA results
pca_raw <- prcomp(t(expr_raw), scale. = FALSE)
pca_corrected <- prcomp(t(expr_corrected), scale. = FALSE)

pca_df_raw <- data.frame(
  PC1 = pca_raw$x[, 1],
  PC2 = pca_raw$x[, 2],
  Dataset = sample_metadata$dataset,
  Sample = rownames(pca_raw$x)
)

pca_df_corrected <- data.frame(
  PC1 = pca_corrected$x[, 1],
  PC2 = pca_corrected$x[, 2],
  Dataset = sample_metadata$dataset,
  Sample = rownames(pca_corrected$x)
)

# Add clinical data for TCGA samples
tcga_samples <- rownames(sample_metadata)[sample_metadata$dataset == "TCGA"]
pca_df_raw$VitalStatus <- NA
pca_df_corrected$VitalStatus <- NA

pca_df_raw$VitalStatus[pca_df_raw$Sample %in% tcga_samples] <- 
  tcga_clinical[pca_df_raw$Sample[pca_df_raw$Sample %in% tcga_samples], "vital_status"]
pca_df_corrected$VitalStatus[pca_df_corrected$Sample %in% tcga_samples] <- 
  tcga_clinical[pca_df_corrected$Sample[pca_df_corrected$Sample %in% tcga_samples], "vital_status"]

cat("Data loaded successfully!\n")

# ==============================================================================
# UI
# ==============================================================================
ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(title = "RNA Biomarker Discovery - Pancreatic Cancer"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
      menuItem("Batch Correction", tabName = "batch", icon = icon("chart-line")),
      menuItem("Gene Explorer", tabName = "genes", icon = icon("dna")),
      menuItem("Model Performance", tabName = "model", icon = icon("brain")),
      menuItem("Predictions", tabName = "predict", icon = icon("magic")),
      menuItem("Download", tabName = "download", icon = icon("download"))
    )
  ),
  
  dashboardBody(
    tabItems(
      
      # TAB 1: Overview
      tabItem(
        tabName = "overview",
        fluidRow(
          box(
            title = "Project Overview", width = 12, status = "primary", solidHeader = TRUE,
            h3("Batch-Harmonized AI for Pancreatic Cancer RNA Data"),
            p("This interactive tool demonstrates a reproducible pipeline for cross-cohort RNA biomarker discovery."),
            hr(),
            h4("Key Features:"),
            tags$ul(
              tags$li("Batch correction using ComBat algorithm"),
              tags$li("Machine learning classification (Random Forest)"),
              tags$li("Cross-cohort validation (TCGA → GEO)"),
              tags$li("Prognostic biomarker identification")
            )
          )
        ),
        fluidRow(
          valueBox(535, "Total Samples", icon = icon("users"), color = "blue", width = 3),
          valueBox(14137, "Genes Analyzed", icon = icon("dna"), color = "green", width = 3),
          valueBox("92.6%", "Model Accuracy", icon = icon("check-circle"), color = "yellow", width = 3),
          valueBox(5, "Top Biomarkers", icon = icon("star"), color = "red", width = 3)
        ),
        fluidRow(
          box(
            title = "Dataset Summary", width = 6, status = "info",
            DTOutput("dataset_table")
          ),
          box(
            title = "Top 5 Biomarker Genes", width = 6, status = "warning",
            DTOutput("top_genes_table")
          )
        )
      ),
      
      # TAB 2: Batch Correction
      tabItem(
        tabName = "batch",
        fluidRow(
          box(
            title = "Batch Correction Effect", width = 12, status = "primary", solidHeader = TRUE,
            p("Compare PCA plots before and after batch correction to see dataset harmonization.")
          )
        ),
        fluidRow(
          box(
            title = "Before Batch Correction", width = 6, status = "danger",
            plotlyOutput("pca_before", height = 500)
          ),
          box(
            title = "After Batch Correction", width = 6, status = "success",
            plotlyOutput("pca_after", height = 500)
          )
        ),
        fluidRow(
          box(
            title = "Interpretation", width = 12, status = "info",
            p(strong("Before:"), "Samples cluster by dataset (technical variation dominates)"),
            p(strong("After:"), "Improved mixing of datasets (biological signal preserved)")
          )
        )
      ),
      
      # TAB 3: Gene Explorer
      tabItem(
        tabName = "genes",
        fluidRow(
          box(
            title = "Gene Expression Explorer", width = 12, status = "primary", solidHeader = TRUE,
            selectInput("selected_gene", "Select Gene:", choices = top_genes, selected = top_genes[1]),
            p("Explore expression patterns of top biomarker genes across datasets and survival groups.")
          )
        ),
        fluidRow(
          box(
            title = "Expression by Dataset", width = 6,
            plotlyOutput("gene_boxplot_dataset", height = 400)
          ),
          box(
            title = "Expression by Survival Status (TCGA)", width = 6,
            plotlyOutput("gene_boxplot_survival", height = 400)
          )
        ),
        fluidRow(
          box(
            title = "Gene Information", width = 12, status = "info",
            htmlOutput("gene_info")
          )
        )
      ),
      
      # TAB 4: Model Performance
      tabItem(
        tabName = "model",
        fluidRow(
          box(
            title = "Random Forest Model Performance", width = 12, status = "primary", solidHeader = TRUE,
            p("Classification task: Predicting survival outcome (Alive vs Dead) in pancreatic cancer patients")
          )
        ),
        fluidRow(
          box(
            title = "Feature Importance (Top 20 Genes)", width = 8,
            plotlyOutput("importance_plot", height = 600)
          ),
          box(
            title = "Model Statistics", width = 4,
            valueBoxOutput("oob_error_box", width = NULL),
            valueBoxOutput("accuracy_box", width = NULL),
            valueBoxOutput("n_trees_box", width = NULL),
            hr(),
            h4("Class Performance:"),
            tableOutput("class_performance")
          )
        )
      ),
      
      # TAB 5: Predictions
      tabItem(
        tabName = "predict",
        fluidRow(
          box(
            title = "Prediction Interface", width = 12, status = "primary", solidHeader = TRUE,
            p("View model predictions on the GEO validation cohort (external dataset)")
          )
        ),
        fluidRow(
          box(
            title = "Prediction Distribution", width = 12,
            plotlyOutput("prediction_histogram", height = 400)
          )
        ),
        fluidRow(
          box(
            title = "Sample Predictions", width = 12,
            DTOutput("predictions_table")
          )
        )
      ),
      
      # TAB 6: Download
      tabItem(
        tabName = "download",
        fluidRow(
          box(
            title = "Download Results", width = 12, status = "primary", solidHeader = TRUE,
            p("Download processed data, models, and results for further analysis.")
          )
        ),
        fluidRow(
          box(
            title = "Available Downloads", width = 12,
            downloadButton("download_top_genes", "Top 50 Genes (CSV)", class = "btn-primary"),
            downloadButton("download_importance", "Feature Importance (CSV)", class = "btn-primary"),
            downloadButton("download_predictions", "Predictions (CSV)", class = "btn-primary"),
            downloadButton("download_report", "Summary Report (TXT)", class = "btn-primary")
          )
        )
      )
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================
server <- function(input, output, session) {
  
  # Dataset table
  output$dataset_table <- renderDT({
    data.frame(
      Dataset = c("TCGA-PAAD", "GSE71729", "Total"),
      Samples = c(sum(sample_metadata$dataset == "TCGA"),
                  sum(sample_metadata$dataset == "GSE71729"),
                  nrow(sample_metadata)),
      Purpose = c("Training", "Validation", "-")
    )
  }, options = list(dom = 't'))
  
  # Top genes table
  output$top_genes_table <- renderDT({
    data.frame(
      Rank = 1:5,
      Gene = top_genes[1:5],
      Importance = round(rf_importance[top_genes[1:5], "MeanDecreaseGini"], 2)
    )
  }, options = list(dom = 't'))
  
  # PCA before
  output$pca_before <- renderPlotly({
    plot_ly(pca_df_raw, x = ~PC1, y = ~PC2, color = ~Dataset, 
            text = ~paste("Sample:", Sample, "<br>Dataset:", Dataset),
            type = 'scatter', mode = 'markers') %>%
      layout(title = "PCA - Before Correction",
             xaxis = list(title = "PC1"),
             yaxis = list(title = "PC2"))
  })
  
  # PCA after
  output$pca_after <- renderPlotly({
    plot_ly(pca_df_corrected, x = ~PC1, y = ~PC2, color = ~Dataset,
            text = ~paste("Sample:", Sample, "<br>Dataset:", Dataset),
            type = 'scatter', mode = 'markers') %>%
      layout(title = "PCA - After Correction",
             xaxis = list(title = "PC1"),
             yaxis = list(title = "PC2"))
  })
  
  # Gene expression by dataset
  output$gene_boxplot_dataset <- renderPlotly({
    gene_data <- data.frame(
      Expression = expr_corrected[input$selected_gene, ],
      Dataset = sample_metadata$dataset
    )
    
    plot_ly(gene_data, y = ~Expression, color = ~Dataset, type = "box") %>%
      layout(title = paste(input$selected_gene, "Expression by Dataset"),
             yaxis = list(title = "Expression Level"))
  })
  
  # Gene expression by survival
  output$gene_boxplot_survival <- renderPlotly({
    tcga_idx <- sample_metadata$dataset == "TCGA"
    tcga_samples_sel <- rownames(sample_metadata)[tcga_idx]
    
    gene_data <- data.frame(
      Expression = expr_corrected[input$selected_gene, tcga_samples_sel],
      VitalStatus = tcga_clinical[tcga_samples_sel, "vital_status"]
    )
    gene_data <- gene_data[!is.na(gene_data$VitalStatus), ]
    
    plot_ly(gene_data, y = ~Expression, color = ~VitalStatus, type = "box") %>%
      layout(title = paste(input$selected_gene, "Expression by Survival"),
             yaxis = list(title = "Expression Level"))
  })
  
  # Gene info
  output$gene_info <- renderUI({
    gene_annotations <- list(
      "LAMC2" = "Laminin Subunit Gamma 2 - Involved in cell adhesion and migration. Overexpressed in invasive cancers.",
      "DKK1" = "Dickkopf WNT Signaling Pathway Inhibitor 1 - Regulates WNT pathway. Elevated in multiple cancer types.",
      "ITGB6" = "Integrin Subunit Beta 6 - Epithelial cell marker. Associated with EMT and cancer progression.",
      "GPRC5A" = "G Protein-Coupled Receptor - Tumor suppressor, frequently downregulated in cancer.",
      "MAL2" = "Mal, T-Cell Differentiation Protein 2 - Vesicle trafficking. Overexpressed in various cancers."
    )
    
    info <- gene_annotations[[input$selected_gene]]
    if (is.null(info)) info <- "Gene information not available."
    
    HTML(paste0("<p><strong>", input$selected_gene, ":</strong> ", info, "</p>"))
  })
  
  # Feature importance plot
  output$importance_plot <- renderPlotly({
    top20 <- head(rf_importance[order(-rf_importance$MeanDecreaseGini), ], 20)
    top20$Gene <- rownames(top20)
    top20 <- top20[order(top20$MeanDecreaseGini), ]
    
    plot_ly(top20, x = ~MeanDecreaseGini, y = ~reorder(Gene, MeanDecreaseGini),
            type = 'bar', orientation = 'h') %>%
      layout(title = "Top 20 Features by Importance",
             xaxis = list(title = "Mean Decrease Gini"),
             yaxis = list(title = ""))
  })
  
  # Model statistics boxes
  output$oob_error_box <- renderValueBox({
    oob <- round(tail(rf_model$err.rate[,1], 1)*100, 2)
    valueBox(paste0(oob, "%"), "OOB Error Rate", icon = icon("exclamation-triangle"), color = "red")
  })
  
  output$accuracy_box <- renderValueBox({
    acc <- round((1 - tail(rf_model$err.rate[,1], 1))*100, 2)
    valueBox(paste0(acc, "%"), "Training Accuracy", icon = icon("check"), color = "green")
  })
  
  output$n_trees_box <- renderValueBox({
    valueBox(rf_model$ntree, "Trees in Forest", icon = icon("tree"), color = "blue")
  })
  
  # Class performance table
  output$class_performance <- renderTable({
    confusion <- rf_model$confusion
    data.frame(
      Class = rownames(confusion)[1:2],
      N = confusion[1:2, 1] + confusion[1:2, 2],
      Accuracy = paste0(round(diag(confusion)[1:2] / 
                              (confusion[1:2, 1] + confusion[1:2, 2]) * 100, 1), "%")
    )
  })
  
  # Prediction histogram
  output$prediction_histogram <- renderPlotly({
    if (file.exists("results/test_predictions.rds")) {
      preds <- readRDS("results/test_predictions.rds")
      
      plot_ly(x = preds$rf_corrected_prob, type = "histogram") %>%
        layout(title = "Distribution of Predictions on GEO Cohort",
               xaxis = list(title = "Predicted Probability (Dead)"),
               yaxis = list(title = "Count"))
    }
  })
  
  # Predictions table
  output$predictions_table <- renderDT({
    if (file.exists("results/test_predictions.rds")) {
      preds <- readRDS("results/test_predictions.rds")
      data.frame(
        Sample = 1:nrow(preds),
        Dataset = preds$dataset,
        RF_Probability = round(preds$rf_corrected_prob, 3),
        XGB_Probability = round(preds$xgb_corrected_prob, 3),
        Predicted_Class = ifelse(preds$rf_corrected_prob > 0.5, "Dead", "Alive")
      )
    }
  }, options = list(pageLength = 10))
  
  # Download handlers
  output$download_top_genes <- downloadHandler(
    filename = "top50_genes.csv",
    content = function(file) {
      write.csv(head(rf_importance[order(-rf_importance$MeanDecreaseGini), ], 50), file)
    }
  )
  
  output$download_importance <- downloadHandler(
    filename = "feature_importance.csv",
    content = function(file) {
      write.csv(rf_importance, file)
    }
  )
  
  output$download_predictions <- downloadHandler(
    filename = "predictions.csv",
    content = function(file) {
      if (file.exists("results/test_predictions.rds")) {
        preds <- readRDS("results/test_predictions.rds")
        write.csv(preds, file, row.names = FALSE)
      }
    }
  )
  
  output$download_report <- downloadHandler(
    filename = "summary_report.txt",
    content = function(file) {
      file.copy("results/manuscript_summary.txt", file)
    }
  )
}

# ==============================================================================
# RUN APP
# ==============================================================================
shinyApp(ui = ui, server = server)