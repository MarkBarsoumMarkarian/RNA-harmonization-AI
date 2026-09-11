options(stringsAsFactors = FALSE)

input_file <- Sys.getenv("V2_COHORTS_FILE", "data/v2/cohorts.rds")
result_dir <- Sys.getenv("V2_RESULTS_DIR", "results/v2")
if (!file.exists(input_file)) stop("Missing ", input_file)
cohorts <- readRDS(input_file)
candidates <- c("LAMC2", "DKK1", "ITGB6", "GPRC5A", "MAL2")

cohort_estimates <- list()
counter <- 1L
for (cohort_name in names(cohorts)) {
  item <- cohorts[[cohort_name]]
  clinical <- item$clinical
  if (any(clinical$specimen_class != "primary_tumor")) {
    stop(cohort_name, ": non-primary specimen detected")
  }
  if (anyDuplicated(clinical$patient_id)) stop(cohort_name, ": duplicated patient")
  if (!identical(colnames(item$expression), clinical$sample_id)) {
    stop(cohort_name, ": expression/clinical mismatch")
  }
  available <- intersect(candidates, rownames(item$expression))
  for (gene in available) {
    value <- as.numeric(scale(item$expression[gene, ]))
    fit <- survival::coxph(survival::Surv(clinical$time, clinical$event) ~ value)
    cohort_estimates[[counter]] <- data.frame(
      cohort = cohort_name,
      gene = gene,
      n = nrow(clinical),
      events = sum(clinical$event),
      log_hr = stats::coef(fit)[1L],
      se = sqrt(stats::vcov(fit)[1L, 1L]),
      stringsAsFactors = FALSE
    )
    counter <- counter + 1L
  }
}
cohort_estimates <- do.call(rbind, cohort_estimates)

random_effects <- function(beta, se) {
  weight_fixed <- 1 / se^2
  fixed <- sum(weight_fixed * beta) / sum(weight_fixed)
  q <- sum(weight_fixed * (beta - fixed)^2)
  df <- length(beta) - 1L
  c_value <- sum(weight_fixed) - sum(weight_fixed^2) / sum(weight_fixed)
  tau2 <- if (df > 0L && c_value > 0) max(0, (q - df) / c_value) else 0
  weight_random <- 1 / (se^2 + tau2)
  pooled <- sum(weight_random * beta) / sum(weight_random)
  pooled_se <- sqrt(1 / sum(weight_random))
  i2 <- if (q > 0 && df > 0L) max(0, (q - df) / q) * 100 else 0
  c(beta = pooled, se = pooled_se, tau2 = tau2, i2 = i2, q = q, df = df)
}

meta <- do.call(rbind, lapply(split(cohort_estimates, cohort_estimates$gene), function(x) {
  estimate <- random_effects(x$log_hr, x$se)
  data.frame(
    gene = x$gene[1L],
    cohorts = nrow(x),
    patients = sum(x$n),
    events = sum(x$events),
    hazard_ratio_per_cohort_sd = exp(estimate["beta"]),
    ci_95_low = exp(estimate["beta"] - 1.96 * estimate["se"]),
    ci_95_high = exp(estimate["beta"] + 1.96 * estimate["se"]),
    p_value = 2 * stats::pnorm(-abs(estimate["beta"] / estimate["se"])),
    i2_percent = estimate["i2"],
    tau2 = estimate["tau2"],
    stringsAsFactors = FALSE
  )
}))
meta$fdr_bh_within_five <- stats::p.adjust(meta$p_value, "BH")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(cohort_estimates, file.path(result_dir, "candidate_cohort_cox.csv"), row.names = FALSE)
write.csv(meta, file.path(result_dir, "candidate_random_effects_meta.csv"), row.names = FALSE)
message("Candidate meta-analysis complete")
print(meta, row.names = FALSE)
