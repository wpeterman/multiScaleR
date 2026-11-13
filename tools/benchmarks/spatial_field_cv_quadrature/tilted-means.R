# Follow-up comparison on the original simulation design and references.
# This study evaluates the first moment only and changes no package behavior.
source("tools/benchmarks/spatial_field_cv_quadrature/run.R")
args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1L] else
  "tools/benchmarks/spatial_field_cv_quadrature/results"
pkgload::load_all(compile = FALSE, quiet = TRUE)
cases <- read.csv(file.path(path, "scenarios.csv"), stringsAsFactors = FALSE)
predictions <- read.csv(file.path(path, "predictions.csv"), stringsAsFactors = FALSE)
cases <- subset(cases, cohort == "Laplace" & family %in% c("poisson", "nbinom2", "Gamma"))
predictions <- split(predictions[predictions$id %in% cases$id, ], predictions$id[predictions$id %in% cases$id])
results <- vector("list", nrow(cases))
for (i in seq_len(nrow(cases))) {
  case <- cases[i, , drop = FALSE]
  d <- simulate_case(case)
  ti <- seq_along(d$train_y)
  ai <- length(ti) + seq_along(d$test_y)
  ti_nodes <- seq_len(d$n_train)
  test_node <- d$groups[ai]
  train_lower <- d$lower[ti_nodes, ti_nodes, drop = FALSE]
  run <- capture(function() {
    if (!is.finite(case$evaluated_analysis)) stop("The original analysis likelihood is unavailable.")
    # For E[exp(t'Z) | y], completing the normal square moves the
    # marginal-likelihood integral to the shifted prior N(t, I).
    # Reexpress it as N(0, I) with an offset to the training predictors.
    nodes <- unique(test_node)
    shifted <- vapply(nodes, function(node) {
      offset <- drop(train_lower %*% d$lower[node, ti_nodes])
      fit <- multiScaleR:::.msr_field_laplace(d$train_y,
        d$train_eta + offset[d$groups[ti]], d$groups[ti], train_lower,
        d$family, d$dispersion)
      fit$log_density - case$evaluated_analysis
    }, 0)
    mean <- exp(d$test_eta + case$field_sd^2/2 + shifted[match(test_node, nodes)])
    if (any(!is.finite(mean))) stop("Tilted response mean is not finite.")
    mean
  })
  ref <- predictions[[case$id]]
  stopifnot(identical(as.numeric(d$test_y), as.numeric(ref$observed)) ||
    isTRUE(all.equal(as.numeric(d$test_y), as.numeric(ref$observed), tolerance = 1e-12)))
  compare <- case$reference_status == "converged" && !is.null(run$value)
  results[[i]] <- cbind(case[rep(1L, length(ai)), c("id", "family", "field_sd", "rho",
    "visits", "layout", "covariance", "reference_status")],
    data.frame(assessment_row = ref$assessment_row, reference_mean = ref$reference_mean,
      gaussian_mean = ref$evaluated_mean, gaussian_relative_error = ref$scaled_abs_mean_error,
      tilted_mean = if (is.null(run$value)) NA_real_ else run$value,
      tilted_relative_error = if (compare) abs(run$value/ref$reference_mean - 1) else NA_real_,
      tilted_status = if (is.null(run$value)) "failed" else "returned", tilted_message = run$message,
      tilted_warnings = run$warnings))
}
result <- do.call(rbind, results)
rownames(result) <- NULL
write.csv(result, file.path(path, "tilted-means.csv"), row.names = FALSE)
summary <- lapply(split(result, list(result$family, result$field_sd), drop = TRUE), function(x) {
  paired <- x[is.finite(x$gaussian_relative_error) & is.finite(x$tilted_relative_error), ]
  q <- function(error) unname(quantile(error, c(.5, .95, 1)))
  g <- q(paired$gaussian_relative_error); t <- q(paired$tilted_relative_error)
  data.frame(family = x$family[1L], field_sd = x$field_sd[1L], rows = nrow(x),
    paired = nrow(paired), failed = length(unique(x$id[x$tilted_status == "failed"])),
    gaussian_p50 = g[1L], gaussian_p95 = g[2L], gaussian_max = g[3L],
    tilted_p50 = t[1L], tilted_p95 = t[2L], tilted_max = t[3L])
})
summary <- do.call(rbind, summary)
rownames(summary) <- NULL
write.csv(summary, file.path(path, "tilted-summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
