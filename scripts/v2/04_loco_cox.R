options(stringsAsFactors = FALSE)
source("scripts/v2/lib_survival.R")

input_file <- Sys.getenv("V2_COHORTS_FILE", "data/v2/cohorts.rds")
result_dir <- Sys.getenv("V2_RESULTS_DIR", "results/v2")
if (!file.exists(input_file)) {
  stop(
    "Missing ", input_file, ". See scripts/v2/README.md for the required schema. ",
    "GSE71729 must not be inserted until auditable patient-level survival time ",
    "and event status have been joined to primary tumours."
  )
}

cohorts <- readRDS(input_file)
if (!is.list(cohorts) || length(cohorts) < 2L || is.null(names(cohorts))) {
  stop("LOCO requires a named list containing at least two cohorts")
}

validate_cohort <- function(x, name) {
  if (!is.matrix(x$expression)) stop(name, ": expression must be a matrix")
  required <- c("sample_id", "patient_id", "time", "event", "specimen_class")
  if (!all(required %in% names(x$clinical))) stop(name, ": missing clinical fields")
  if (!identical(colnames(x$expression), x$clinical$sample_id)) {
    stop(name, ": expression/clinical order mismatch")
  }
  if (anyDuplicated(rownames(x$expression))) stop(name, ": duplicated genes")
  if (anyDuplicated(x$clinical$patient_id)) stop(name, ": duplicated patients")
  if (any(x$clinical$specimen_class != "primary_tumor")) {
    stop(name, ": non-primary specimens are forbidden in prognostic validation")
  }
  if (any(!is.finite(x$clinical$time) | x$clinical$time <= 0)) {
    stop(name, ": survival times must be positive and finite")
  }
  if (any(!x$clinical$event %in% c(0L, 1L))) stop(name, ": invalid event encoding")
  if (sum(x$clinical$event) < 10L) stop(name, ": fewer than ten events")
  invisible(TRUE)
}

for (name in names(cohorts)) validate_cohort(cohorts[[name]], name)

# This is a platform feature dictionary, not a fitted transformation. No test
# values or outcomes are used. For a locked final model, replace it with a
# predeclared assay panel before opening the external cohort.
common_genes <- Reduce(intersect, lapply(cohorts, function(x) rownames(x$expression)))
if (length(common_genes) < 100L) stop("Fewer than 100 shared, uniquely mapped genes")

# Each sample is transformed independently. Unlike joint ComBat or cohort-wise
# z-scoring, this does not borrow the held-out cohort distribution.
for (name in names(cohorts)) {
  cohorts[[name]]$expression <- sample_percentile_rank(
    cohorts[[name]]$expression[common_genes, , drop = FALSE]
  )
  rownames(cohorts[[name]]$expression) <- common_genes
  cohorts[[name]]$clinical$cohort <- name
}

predictions <- list()
selections <- list()
for (held_out in names(cohorts)) {
  message("Holding out cohort: ", held_out)
  train_names <- setdiff(names(cohorts), held_out)
  train_expr <- do.call(cbind, lapply(cohorts[train_names], `[[`, "expression"))
  train_clinical <- do.call(rbind, lapply(cohorts[train_names], `[[`, "clinical"))
  test_expr <- cohorts[[held_out]]$expression
  test_clinical <- cohorts[[held_out]]$clinical

  parameters <- tune_ridge_cox(
    train_expr,
    train_clinical,
    feature_counts = c(5L, 10L, 20L),
    theta_values = c(0.1, 1, 10),
    top_mad = min(1000L, length(common_genes)),
    seed = 9000L + match(held_out, names(cohorts))
  )
  fitted <- fit_one_split(
    train_expr,
    test_expr,
    train_clinical,
    parameters,
    top_mad = min(1000L, length(common_genes))
  )
  predictions[[held_out]] <- data.frame(
    held_out_cohort = held_out,
    sample_id = test_clinical$sample_id,
    patient_id = test_clinical$patient_id,
    time = test_clinical$time,
    event = test_clinical$event,
    risk = fitted$risk,
    risk_z = fitted$risk_z,
    n_features = parameters$n_features,
    theta = parameters$theta,
    stringsAsFactors = FALSE
  )
  selections[[held_out]] <- data.frame(
    held_out_cohort = held_out,
    gene = fitted$genes,
    coefficient = as.numeric(fitted$coefficient),
    stringsAsFactors = FALSE
  )
}

predictions <- do.call(rbind, predictions)
selections <- do.call(rbind, selections)
performance <- do.call(rbind, lapply(split(predictions, predictions$held_out_cohort), function(x) {
  estimate <- cindex_value(x$time, x$event, x$risk_z)
  ci <- bootstrap_cindex(x$time, x$event, x$risk_z, repetitions = 1000L, seed = 10101L)
  data.frame(
    held_out_cohort = x$held_out_cohort[1L],
    n_patients = nrow(x),
    n_events = sum(x$event),
    c_index = estimate,
    ci_95_low = ci[1L],
    ci_95_high = ci[2L],
    stringsAsFactors = FALSE
  )
}))

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(predictions, file.path(result_dir, "loco_predictions.csv"), row.names = FALSE)
write.csv(performance, file.path(result_dir, "loco_performance.csv"), row.names = FALSE)
write.csv(selections, file.path(result_dir, "loco_selected_genes.csv"), row.names = FALSE)
message("LOCO validation complete")
print(performance, row.names = FALSE)
