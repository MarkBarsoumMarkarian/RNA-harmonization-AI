options(stringsAsFactors = FALSE)
source("scripts/v2/lib_survival.R")

input_file <- Sys.getenv("V2_COHORTS_FILE", "data/v2/cohorts.rds")
result_dir <- Sys.getenv("V2_RESULTS_DIR", "results/v2")
cohorts <- readRDS(input_file)
if (!all(c("TCGA_PAAD", "GSE85916") %in% names(cohorts))) {
  stop("This fixed-panel test requires TCGA_PAAD and GSE85916")
}

train <- cohorts$TCGA_PAAD
test <- cohorts$GSE85916
candidates <- c("LAMC2", "DKK1", "ITGB6", "GPRC5A", "MAL2")
common_genes <- intersect(rownames(train$expression), rownames(test$expression))
if (!all(candidates %in% common_genes)) stop("One or more panel genes are not shared")

# The transformation is sample-wise and fixed before outcome modelling. Test
# samples do not contribute to training centering, scaling, coefficients, or
# any parameter choice.
train_rank <- sample_percentile_rank(train$expression[common_genes, , drop = FALSE])
test_rank <- sample_percentile_rank(test$expression[common_genes, , drop = FALSE])
rownames(train_rank) <- common_genes
rownames(test_rank) <- common_genes

x_train <- t(train_rank[candidates, , drop = FALSE])
x_test <- t(test_rank[candidates, , drop = FALSE])
center <- colMeans(x_train)
spread <- apply(x_train, 2L, stats::sd)
spread[!is.finite(spread) | spread == 0] <- 1
x_train <- sweep(sweep(x_train, 2L, center, "-"), 2L, spread, "/")
x_test <- sweep(sweep(x_test, 2L, center, "-"), 2L, spread, "/")

fit <- survival::coxph(
  survival::Surv(train$clinical$time, train$clinical$event) ~ x_train,
  ties = "efron",
  x = TRUE
)
coefficient <- stats::coef(fit)
train_risk <- drop(x_train %*% coefficient)
test_risk <- drop(x_test %*% coefficient)

test_cindex <- cindex_value(test$clinical$time, test$clinical$event, test_risk)
test_ci <- bootstrap_cindex(
  test$clinical$time, test$clinical$event, test_risk,
  repetitions = 1000L, seed = 71207L
)
train_cindex <- cindex_value(train$clinical$time, train$clinical$event, train_risk)

performance <- data.frame(
  model = "Original five genes; multivariable Cox; TCGA-fitted coefficients",
  training_cohort = "TCGA_PAAD",
  held_out_cohort = "GSE85916",
  held_out_patients = nrow(test$clinical),
  held_out_events = sum(test$clinical$event),
  apparent_training_c_index = train_cindex,
  held_out_c_index = test_cindex,
  held_out_ci_95_low = test_ci[1L],
  held_out_ci_95_high = test_ci[2L],
  stringsAsFactors = FALSE
)
coefficients <- data.frame(
  gene = sub("^x_train", "", names(coefficient)),
  coefficient = as.numeric(coefficient),
  hazard_ratio_per_tcga_sd = exp(as.numeric(coefficient)),
  stringsAsFactors = FALSE
)
predictions <- data.frame(
  sample_id = test$clinical$sample_id,
  patient_id = test$clinical$patient_id,
  time = test$clinical$time,
  event = test$clinical$event,
  risk = test_risk,
  stringsAsFactors = FALSE
)

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(performance, file.path(result_dir, "fixed_five_external_performance.csv"), row.names = FALSE)
write.csv(coefficients, file.path(result_dir, "fixed_five_tcga_coefficients.csv"), row.names = FALSE)
write.csv(predictions, file.path(result_dir, "fixed_five_external_predictions.csv"), row.names = FALSE)
message("Fixed five-gene external test complete")
print(performance, row.names = FALSE)
