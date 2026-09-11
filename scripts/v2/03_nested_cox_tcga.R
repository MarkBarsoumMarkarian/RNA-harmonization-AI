options(stringsAsFactors = FALSE)
source("scripts/v2/lib_survival.R")

dir.create("results/v2", recursive = TRUE, showWarnings = FALSE)
data <- readRDS("data/v2/tcga_cbioportal.rds")
expr <- data$expression
clinical <- data$clinical

if (!identical(colnames(expr), clinical$sample_id)) stop("Expression/clinical mismatch")
if (any(clinical$specimen_class != "primary_tumor")) stop("Non-primary sample detected")
if (anyDuplicated(clinical$patient_id)) stop("Duplicated patients detected")
if (any(!is.finite(clinical$time) | clinical$time <= 0)) stop("Invalid follow-up time")
if (any(!clinical$event %in% c(0L, 1L))) stop("Invalid event encoding")

outer_folds <- 5L
repeats <- as.integer(Sys.getenv("V2_REPEATS", "5"))
top_mad <- as.integer(Sys.getenv("V2_TOP_MAD", "1000"))
feature_counts <- c(5L, 10L, 20L)
theta_values <- c(0.1, 1, 10)

predictions <- list()
selections <- list()
counter <- 1L

for (repeat_id in seq_len(repeats)) {
  message("Nested CV repeat ", repeat_id, "/", repeats)
  fold <- stratified_folds(clinical$event, outer_folds, 4100L + repeat_id)
  for (fold_id in seq_len(outer_folds)) {
    test <- fold == fold_id
    train <- !test
    parameters <- tune_ridge_cox(
      expr[, train, drop = FALSE],
      clinical[train, , drop = FALSE],
      feature_counts = feature_counts,
      theta_values = theta_values,
      top_mad = top_mad,
      seed = 51000L + 100L * repeat_id + fold_id
    )
    fitted <- fit_one_split(
      expr[, train, drop = FALSE],
      expr[, test, drop = FALSE],
      clinical[train, , drop = FALSE],
      parameters,
      top_mad = top_mad
    )
    predictions[[counter]] <- data.frame(
      sample_id = clinical$sample_id[test],
      patient_id = clinical$patient_id[test],
      time = clinical$time[test],
      event = clinical$event[test],
      repeat_id = repeat_id,
      fold_id = fold_id,
      risk = fitted$risk,
      risk_z = fitted$risk_z,
      n_features = parameters$n_features,
      theta = parameters$theta,
      inner_cindex = parameters$mean_inner_cindex,
      stringsAsFactors = FALSE
    )
    selections[[counter]] <- data.frame(
      repeat_id = repeat_id,
      fold_id = fold_id,
      gene = fitted$genes,
      coefficient = as.numeric(fitted$coefficient),
      stringsAsFactors = FALSE
    )
    counter <- counter + 1L
  }
}

predictions <- do.call(rbind, predictions)
selections <- do.call(rbind, selections)

# Repeated outer-CV predictions are averaged per patient after putting each
# split's linear predictor on its training-risk scale.
oof <- aggregate(risk_z ~ sample_id + patient_id + time + event, predictions, mean)
oof <- oof[match(clinical$sample_id, oof$sample_id), ]

c_index <- cindex_value(oof$time, oof$event, oof$risk_z)
ci <- bootstrap_cindex(oof$time, oof$event, oof$risk_z, repetitions = 1000L, seed = 8675309L)

performance <- data.frame(
  analysis = "Repeated nested 5-fold CV; expression-only ridge Cox",
  validation = "Internal TCGA only",
  n_patients = nrow(oof),
  n_events = sum(oof$event),
  c_index = c_index,
  ci_95_low = ci[1L],
  ci_95_high = ci[2L],
  repeats = repeats,
  outer_folds = outer_folds,
  top_mad_inside_training = top_mad,
  stringsAsFactors = FALSE
)

selection_frequency <- aggregate(
  list(selected_outer_fits = selections$gene),
  list(gene = selections$gene),
  length
)
selection_frequency$selection_frequency <-
  selection_frequency$selected_outer_fits / (repeats * outer_folds)
selection_frequency <- selection_frequency[
  order(-selection_frequency$selection_frequency, selection_frequency$gene),
]

original_five <- intersect(c("LAMC2", "DKK1", "ITGB6", "GPRC5A", "MAL2"), rownames(expr))
candidate_results <- lapply(original_five, function(gene) {
  value <- as.numeric(scale(expr[gene, ]))
  fit <- survival::coxph(survival::Surv(clinical$time, clinical$event) ~ value)
  coefficient <- stats::coef(fit)[1L]
  se <- sqrt(stats::vcov(fit)[1L, 1L])
  data.frame(
    gene = gene,
    hazard_ratio_per_sd = exp(coefficient),
    ci_95_low = exp(coefficient - 1.96 * se),
    ci_95_high = exp(coefficient + 1.96 * se),
    p_value = summary(fit)$coefficients[1L, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
})
candidate_results <- do.call(rbind, candidate_results)
candidate_results$fdr_bh_within_five <- stats::p.adjust(candidate_results$p_value, "BH")

write.csv(performance, "results/v2/tcga_nested_cv_performance.csv", row.names = FALSE)
write.csv(predictions, "results/v2/tcga_all_outer_predictions.csv", row.names = FALSE)
write.csv(oof, "results/v2/tcga_oof_predictions.csv", row.names = FALSE)
write.csv(selection_frequency, "results/v2/tcga_selected_gene_stability.csv", row.names = FALSE)
write.csv(candidate_results, "results/v2/original_five_gene_cox.csv", row.names = FALSE)

manifest <- c(
  paste("Run time:", format(Sys.time(), tz = "UTC"), "UTC"),
  paste("R:", R.version.string),
  paste("survival:", as.character(utils::packageVersion("survival"))),
  paste("Input source:", data$provenance$source),
  paste("Input purpose:", data$provenance$purpose),
  paste("Patients:", nrow(clinical)),
  paste("Events:", sum(clinical$event)),
  paste("Genes:", nrow(expr)),
  paste("Outer folds:", outer_folds),
  paste("Repeats:", repeats),
  paste("Training-only MAD pool:", top_mad),
  paste("Candidate feature counts:", paste(feature_counts, collapse = ",")),
  paste("Candidate ridge theta:", paste(theta_values, collapse = ",")),
  "Status: provisional internal validation; not external validation"
)
writeLines(manifest, "results/v2/run_manifest.txt")

message("Nested TCGA analysis complete")
print(performance, row.names = FALSE)
message("Original five-gene univariate Cox results (exploratory, not external validation):")
print(candidate_results, row.names = FALSE)
