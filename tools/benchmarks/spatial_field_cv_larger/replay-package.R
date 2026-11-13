pkgload::load_all(compile = FALSE, quiet = TRUE)
source("tools/benchmarks/spatial_field_cv_larger/run.R")
path <- "tools/benchmarks/spatial_field_cv_larger/results"
cases <- read.csv(file.path(path, "design.csv"))
scenarios <- read.csv(file.path(path, "scenarios.csv"))
predictions <- read.csv(file.path(path, "predictions.csv"))
result <- lapply(seq_len(nrow(cases)), function(i) {
  case <- cases[i, , drop = FALSE]
  if (case$cohort != "larger") return(NULL)
  d <- simulate_larger(case)
  x <- multiScaleR:::.msr_field_laplace_block(d$train_y, d$test_y, d$train_eta, d$test_eta,
    d$train_xy, d$test_xy,
    case$covariance, case$field_sd, d$range, case$family, d$dispersion,
    tilted_log_moments = TRUE)
  reference <- predictions[predictions$id == case$id, ]
  data.frame(id = case$id, score_difference = abs(x$joint$neg_log_score -
    scenarios$fitted_loss[match(case$id, scenarios$id)]),
    mean_scaled_difference = max(abs(x$prediction - reference$fitted_mean) /
      pmax(1, abs(reference$fitted_mean))),
    variance_scaled_difference = max(abs(x$observation_predictive_SE^2 - reference$fitted_variance) /
      pmax(1, abs(reference$fitted_variance))))
})
table <- do.call(rbind, result)
stopifnot(nrow(table) == 28L, all(is.finite(as.matrix(table[-1L]))),
  max(as.matrix(table[-1L])) < 1e-10)
print(vapply(table[-1L], max, numeric(1)))
