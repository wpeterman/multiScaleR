# Post-study investigation of the largest Gaussian annual joint-score difference.
# Invoke with the run.R qualification OUTPUT assess arguments.
source("tools/benchmarks/bayesian_scales/qualification_20261112/run.R")
stopifnot(mode == "qualification", selector == "assess")
eligible <- summary[summary$family == "gaussian" & summary$design == "annual" &
  summary$status == "eligible", ]
refs <- eligible[eligible$method == "nuts", ]
alt <- eligible[eligible$method != "nuts", ]
distance <- abs(alt$joint_log_score - refs$joint_log_score[match(alt$case, refs$case)])
case <- plan[match(alt$case[which.max(distance)], plan$case), ]
input <- simulate_case(case)
training <- input$data[input$analysis, , drop = FALSE]

normal_quadrature <- function(n) {
  jacobi <- matrix(0, n, n)
  pairs <- cbind(seq_len(n - 1L), seq.int(2L, n))
  jacobi[pairs] <- jacobi[pairs[, 2:1]] <- sqrt(seq_len(n - 1L))
  eig <- eigen(jacobi, symmetric = TRUE)
  list(z = eig$values, w = eig$vectors[1L, ]^2)
}
logsum <- function(x) {
  maximum <- max(x)
  maximum + log(sum(exp(x - maximum)))
}
joint_loss <- function(fit, nodes) {
  if (is.null(fit$importance)) {
    draw <- fit$draws
    w <- rep(1 / nrow(draw$beta), nrow(draw$beta))
  } else {
    w <- fit$importance$weights
    draw <- ns$.msr_bayes_parameter_blocks(fit$importance$candidates, fit$parameters)
  }
  use <- which(w > 0)
  draw <- lapply(draw, function(x) x[use, , drop = FALSE])
  w <- w[use] / sum(w[use])
  ids <- input$assessment
  spec <- fit$specification
  sm <- spec$scales$habitat
  term <- sm$random[[1L]]
  stopifnot(length(sm$beta_index) == 1L, length(spec$response$random) == 0L,
    input$inputs$kernel == "gaussian", length(ids) == 3L)
  base <- draw$scale_beta[, sm$beta_index]
  sd <- exp(draw$log_scale_sd[, term$sd_index])
  year <- match(as.character(input$data$year[ids[1L]]), term$levels)
  stopifnot(!is.na(year), length(unique(input$data$year[ids[-1L]])) == 1L,
    !as.character(input$data$year[ids[2L]]) %in% term$levels)
  existing_scale <- base + draw$scale_z[, term$z_index[year]] * sd
  d <- input$inputs$binned$bin_rep
  d2 <- d^2 - min(d)^2
  V <- input$inputs$binned$vsum$habitat[ids, , drop = FALSE]
  C <- input$inputs$binned$csum$habitat[ids, , drop = FALSE]
  one <- function(log_scale, j) {
    weights <- exp(-outer(exp(-2 * log_scale), d2) / 2)
    raw <- drop(weights %*% V[j, ]) / drop(weights %*% C[j, ])
    x <- (raw - input$reference[["mean"]]) / input$reference[["sd"]]
    mu <- draw$beta[, 1L] + draw$beta[, 2L] * x + input$data$offset[ids[j]]
    dnorm(input$data$y[ids[j]], mu, exp(draw$log_dispersion[, 1L]), log = TRUE)
  }
  known <- one(existing_scale, 1L) + log(w)
  integrand <- vapply(nodes$z, function(z) {
    shared_scale <- base + z * sd
    logsum(known + one(shared_scale, 2L) + one(shared_scale, 3L))
  }, numeric(1))
  -logsum(integrand + log(nodes$w))
}
nodes64 <- normal_quadrature(64L)
nodes128 <- normal_quadrature(128L)
nodes256 <- normal_quadrature(256L)
res <- lapply(c(methods, "nuts_independent"), function(method) {
  is_nuts <- method %in% c("nuts", "nuts_independent")
  ind <- method == "nuts_independent"
  fit <- if (is_nuts) multiScale_bayes(input$formula, input$inputs, training,
    scales = input$scales, family = input$family, row_id = "key", priors = input$priors,
    method = "sample", chains = 4L, iter = if (ind) 8000L else 3000L,
    warmup = if (ind) 4000L else 1500L, seed = case$seed + if (ind) 11L else 1L,
    control = list(adapt_delta = .99, max_treedepth = 12L))
  else multiScale_bayes(input$formula, input$inputs, training,
    scales = input$scales, family = input$family, row_id = "key", priors = input$priors,
    method = "marginal", seed = case$seed + 2L, laplace_draws = 256L,
    control = list(correction = "importance", conditional_beta = method == "conditional_beta",
      importance_min = 32768L, importance_max = 65536L, importance_df = 5,
      importance_refinements = 1L))
  passed <- if (is_nuts) diagnose_bayes_scales(fit)$passed else fit$diagnostics$passed
  primary <- eligible$joint_log_score[eligible$case == case$case & eligible$method == method]
  data.frame(case = case$case, method = method, passed,
    primary_mc_score = if (length(primary)) primary else NA_real_,
    quadrature_64 = joint_loss(fit, nodes64),
    quadrature_128 = joint_loss(fit, nodes128), quadrature_256 = joint_loss(fit, nodes256))
})
res <- do.call(rbind, res)
utils::write.csv(res, file.path(output, "tail-diagnostic.csv"), row.names = FALSE)
print(res, row.names = FALSE)
