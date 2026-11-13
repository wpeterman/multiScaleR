# Resolve the largest repeated-site Poisson predictive-score discrepancy.
# Source the original frozen runner and retain the primary comparisons.
source("tools/benchmarks/bayesian_scales/qualification_20261112/run.R")
stopifnot(mode == "qualification", selector == "assess")
case <- plan[match("poisson-repeated-20", plan$case), ]
input <- simulate_case(case)
train <- input$data[input$analysis, , drop = FALSE]
future <- input$data[input$assessment, , drop = FALSE]
stopifnot(length(unique(future$year[-1L])) == 1L,
  length(unique(future$site[-1L])) == 1L, nrow(future) == 3L)

normal_rule <- function(n) {
  J <- matrix(0, n, n)
  i <- cbind(seq_len(n - 1L), seq.int(2L, n))
  J[i] <- J[i[, 2:1]] <- sqrt(seq_len(n - 1L))
  E <- eigen(J, symmetric = TRUE)
  list(x = E$values, w = E$vectors[1L, ]^2)
}
log_sum <- function(x) {
  peak <- max(x)
  peak + log(sum(exp(x - peak)))
}

score <- function(fit, order) {
  sm <- fit$specification$scales$habitat
  st <- sm$random[[1L]]
  rt <- fit$specification$response$random[[1L]]
  iy <- match(as.character(future$year[1L]), st$levels)
  is <- match(as.character(future$site[1L]), rt$levels)
  stopifnot(!anyNA(c(iy, is)), length(sm$beta_index) == 2L)
  if (is.null(fit$importance)) {
    draw <- fit$draws
    weight <- rep(1 / nrow(draw$beta), nrow(draw$beta))
  } else {
    draw <- ns$.msr_bayes_parameter_blocks(fit$importance$candidates, fit$parameters)
    weight <- fit$importance$weights
  }
  idx <- weight > 0
  draw <- lapply(draw, function(x) x[idx, , drop = FALSE])
  weight <- weight[idx] / sum(weight[idx])
  intercept <- draw$scale_beta[, sm$beta_index[1L]]
  trend <- draw$scale_beta[, sm$beta_index[2L]]
  scale_sd <- exp(draw$log_scale_sd[, st$sd_index])
  response_sd <- exp(draw$log_response_sd[, rt$sd_index])
  a <- draw$beta[, 1L]
  b <- draw$beta[, 2L]
  d <- input$inputs$binned$bin_rep
  d <- d - min(d)
  V <- input$inputs$binned$vsum$habitat[input$assessment, , drop = FALSE]
  C <- input$inputs$binned$csum$habitat[input$assessment, , drop = FALSE]
  fixed_eta <- function(row, log_scale) {
    W <- exp(-outer(exp(-log_scale), d))
    x <- drop(W %*% V[row, ]) / drop(W %*% C[row, ])
    x <- (x - input$reference[["mean"]]) / input$reference[["sd"]]
    a + b * x + future$offset[row]
  }
  known_scale <- intercept + trend * future$time[1L] +
    scale_sd * draw$scale_z[, st$z_index[iy]]
  known_eta <- fixed_eta(1L, known_scale) + response_sd * draw$response_z[, rt$z_index[is]]
  known <- stats::dpois(future$y[1L], exp(known_eta), log = TRUE) + log(weight)
  rule <- normal_rule(order)
  terms <- vapply(rule$x, function(z) {
    future_scale <- intercept + trend * future$time[2L] + scale_sd * z
    eta2 <- fixed_eta(2L, future_scale)
    eta3 <- fixed_eta(3L, future_scale)
    # Given both new effects, the two counts are independent Poisson draws.
    # The shared response shift adds to both log rates.
    base <- known + future$y[2L] * eta2 + future$y[3L] * eta3 -
      sum(lgamma(future$y[2:3] + 1))
    rate <- exp(eta2) + exp(eta3)
    k <- sum(future$y[2:3])
    inner <- vapply(rule$x, function(z2) {
      shift <- response_sd * z2
      log_sum(base + k * shift - rate * exp(shift))
    }, numeric(1))
    log_sum(inner + log(rule$w))
  }, numeric(1))
  -log_sum(terms + log(rule$w))
}

out <- file.path(output, "poisson-tail-diagnostic")
dir.create(out, showWarnings = FALSE)
for (method in c(methods, "nuts_independent")) {
  file <- file.path(out, paste0(method, ".csv"))
  if (file.exists(file)) next
  nuts <- method %in% c("nuts", "nuts_independent")
  independent <- method == "nuts_independent"
  fit <- if (nuts) multiScale_bayes(input$formula, input$inputs, train,
    scales = input$scales, family = input$family, row_id = "key", priors = input$priors,
    method = "sample", chains = 4L, iter = if (independent) 12000L else 3000L,
    warmup = if (independent) 6000L else 1500L,
    seed = case$seed + if (independent) 11L else 1L,
    control = list(adapt_delta = .99, max_treedepth = 12L))
  else multiScale_bayes(input$formula, input$inputs, train,
    scales = input$scales, family = input$family, row_id = "key", priors = input$priors,
    method = "marginal", seed = case$seed + 2L, laplace_draws = 256L,
    control = list(correction = "importance", conditional_beta = method == "conditional_beta",
      importance_min = 32768L, importance_max = 65536L, importance_df = 5,
      importance_refinements = 1L))
  eligible <- isTRUE(fit$diagnostics$passed) &&
    (!nuts || isTRUE(diagnose_bayes_scales(fit)$passed))
  orders <- c(24L, 48L, 72L)
  results <- data.frame(case = case$case, method = method, eligible = eligible,
    order = orders, neg_log_score = vapply(orders, function(n) score(fit, n), numeric(1)))
  utils::write.csv(results, file, row.names = FALSE)
  print(results, row.names = FALSE)
}

