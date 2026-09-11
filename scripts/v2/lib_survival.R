if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The recommended R package 'survival' is required")
}

cindex_value <- function(time, event, risk) {
  keep <- is.finite(time) & !is.na(event) & is.finite(risk)
  if (sum(keep) < 3L || sum(event[keep]) < 2L) return(NA_real_)
  fit <- survival::concordance(
    survival::Surv(time[keep], event[keep]) ~ risk[keep],
    reverse = TRUE
  )
  as.numeric(fit$concordance)
}

stratified_folds <- function(event, k = 5L, seed = 1L) {
  set.seed(seed)
  fold <- integer(length(event))
  for (status in sort(unique(event))) {
    idx <- sample(which(event == status))
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold
}

sample_percentile_rank <- function(expr) {
  ranked <- apply(expr, 2L, rank, ties.method = "average", na.last = "keep")
  if (is.null(dim(ranked))) ranked <- matrix(ranked, ncol = 1L)
  ranked / (nrow(expr) + 1)
}

median_impute <- function(train_x, test_x) {
  med <- apply(train_x, 2L, stats::median, na.rm = TRUE)
  med[!is.finite(med)] <- 0
  for (j in seq_len(ncol(train_x))) {
    train_x[!is.finite(train_x[, j]), j] <- med[j]
    test_x[!is.finite(test_x[, j]), j] <- med[j]
  }
  list(train = train_x, test = test_x)
}

screen_training_features <- function(
  train_expr,
  test_expr,
  time,
  event,
  cohort = NULL,
  top_mad = 1000L
) {
  if (!identical(rownames(train_expr), rownames(test_expr))) {
    stop("Training/test genes are not aligned")
  }
  x_train <- t(train_expr)
  x_test <- t(test_expr)

  mad_values <- apply(x_train, 2L, stats::mad, na.rm = TRUE)
  mad_values[!is.finite(mad_values)] <- -Inf
  keep_n <- min(as.integer(top_mad), sum(mad_values > 0))
  if (keep_n < 2L) stop("Too few non-constant training genes")
  keep <- names(sort(mad_values, decreasing = TRUE))[seq_len(keep_n)]
  x_train <- x_train[, keep, drop = FALSE]
  x_test <- x_test[, keep, drop = FALSE]

  imputed <- median_impute(x_train, x_test)
  x_train <- imputed$train
  x_test <- imputed$test

  center <- colMeans(x_train)
  scale <- apply(x_train, 2L, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  x_train <- sweep(sweep(x_train, 2L, center, "-"), 2L, scale, "/")
  x_test <- sweep(sweep(x_test, 2L, center, "-"), 2L, scale, "/")

  if (is.null(cohort) || length(unique(cohort)) == 1L) {
    null_fit <- survival::coxph(survival::Surv(time, event) ~ 1)
  } else {
    cohort_factor <- factor(cohort)
    null_fit <- survival::coxph(
      survival::Surv(time, event) ~ survival::strata(cohort_factor)
    )
  }
  martingale <- residuals(null_fit, type = "martingale")
  score <- abs(drop(crossprod(x_train, martingale)))
  ranking <- names(sort(score, decreasing = TRUE))

  list(train = x_train, test = x_test, ranking = ranking)
}

fit_ridge_cox <- function(x, time, event, theta, cohort = NULL) {
  if (is.null(cohort) || length(unique(cohort)) == 1L) {
    fit <- survival::coxph(
      survival::Surv(time, event) ~ survival::ridge(x, theta = theta),
      ties = "efron",
      singular.ok = TRUE
    )
  } else {
    cohort_factor <- factor(cohort)
    fit <- survival::coxph(
      survival::Surv(time, event) ~
        survival::ridge(x, theta = theta) + survival::strata(cohort_factor),
      ties = "efron",
      singular.ok = TRUE
    )
  }
  fit
}

tune_ridge_cox <- function(
  expr,
  clinical,
  feature_counts = c(5L, 10L, 20L),
  theta_values = c(0.1, 1, 10),
  top_mad = 1000L,
  seed = 1L
) {
  k_inner <- min(3L, sum(clinical$event == 1L), sum(clinical$event == 0L))
  if (k_inner < 2L) stop("Insufficient events/censoring for inner validation")
  folds <- stratified_folds(clinical$event, k_inner, seed)
  grid <- expand.grid(
    n_features = as.integer(feature_counts),
    theta = as.numeric(theta_values),
    stringsAsFactors = FALSE
  )
  scores <- matrix(NA_real_, nrow(grid), k_inner)

  for (fold_id in seq_len(k_inner)) {
    validation <- folds == fold_id
    training <- !validation
    prepared <- screen_training_features(
      expr[, training, drop = FALSE],
      expr[, validation, drop = FALSE],
      clinical$time[training],
      clinical$event[training],
      clinical$cohort[training],
      top_mad = top_mad
    )
    for (g in seq_len(nrow(grid))) {
      genes <- head(prepared$ranking, grid$n_features[g])
      fit <- try(
        fit_ridge_cox(
          prepared$train[, genes, drop = FALSE],
          clinical$time[training],
          clinical$event[training],
          grid$theta[g],
          clinical$cohort[training]
        ),
        silent = TRUE
      )
      if (inherits(fit, "try-error")) next
      risk <- drop(prepared$test[, genes, drop = FALSE] %*% stats::coef(fit))
      scores[g, fold_id] <- cindex_value(
        clinical$time[validation], clinical$event[validation], risk
      )
    }
  }
  grid$mean_inner_cindex <- rowMeans(scores, na.rm = TRUE)
  grid$mean_inner_cindex[!is.finite(grid$mean_inner_cindex)] <- -Inf
  # Deterministic parsimony tie-break: best C-index, fewer genes, stronger penalty.
  ordering <- order(-grid$mean_inner_cindex, grid$n_features, -grid$theta)
  grid[ordering[1L], , drop = FALSE]
}

fit_one_split <- function(
  train_expr,
  test_expr,
  train_clinical,
  parameters,
  top_mad = 1000L
) {
  prepared <- screen_training_features(
    train_expr,
    test_expr,
    train_clinical$time,
    train_clinical$event,
    train_clinical$cohort,
    top_mad = top_mad
  )
  genes <- head(prepared$ranking, parameters$n_features)
  fit <- fit_ridge_cox(
    prepared$train[, genes, drop = FALSE],
    train_clinical$time,
    train_clinical$event,
    parameters$theta,
    train_clinical$cohort
  )
  coefficient <- stats::coef(fit)
  train_risk <- drop(prepared$train[, genes, drop = FALSE] %*% coefficient)
  test_risk <- drop(prepared$test[, genes, drop = FALSE] %*% coefficient)
  train_sd <- stats::sd(train_risk)
  if (!is.finite(train_sd) || train_sd == 0) train_sd <- 1
  risk_z <- (test_risk - mean(train_risk)) / train_sd
  list(risk = test_risk, risk_z = risk_z, genes = genes, coefficient = coefficient)
}

bootstrap_cindex <- function(time, event, risk, repetitions = 1000L, seed = 1L) {
  set.seed(seed)
  estimates <- replicate(repetitions, {
    index <- sample.int(length(time), replace = TRUE)
    cindex_value(time[index], event[index], risk[index])
  })
  stats::quantile(estimates[is.finite(estimates)], c(0.025, 0.975), names = FALSE)
}
