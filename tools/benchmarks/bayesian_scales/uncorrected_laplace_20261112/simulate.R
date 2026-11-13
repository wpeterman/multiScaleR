# Design adapted from the frozen sibling qualification runner.
simulate_case <- function(case) {
  set.seed(case$seed)
  repeated <- case$design == "repeated"
  kernel <- if (repeated) "exp" else "gaussian"
  n_year <- if (repeated) 4L else 3L
  sites <- if (repeated) 9L else 6L
  train <- expand.grid(site = seq_len(sites), year = seq_len(n_year))
  if (repeated) train <- train[(train$site + train$year) %% 5L != 0, , drop = FALSE]
  future <- data.frame(site = c(1L, sites + 1L, sites + 1L),
    year = c(n_year, n_year + 1L, n_year + 1L))
  data <- rbind(train, future)
  data$time <- (data$year - 2) / 2
  data$key <- sprintf("row%02d", seq_len(nrow(data)))
  data$offset <- runif(nrow(data), -.2, .2)
  data$trials <- sample(4:9, nrow(data), replace = TRUE)
  d <- seq(.025, 1, length.out = 12L)
  n <- nrow(data)
  count <- matrix(sample(2:6, n * length(d), replace = TRUE), n)
  habitat <- outer(rnorm(n), sin(8 * d)) + outer(rnorm(n), cos(12 * d)) +
    matrix(rnorm(n * length(d), sd = .25), n)
  powers <- if (kernel == "gaussian") 2 else 1
  means <- function(log_scale) {
    w <- exp(-outer(exp(-powers * log_scale), d^powers - min(d)^powers) / powers)
    rowSums(w * habitat * count) / rowSums(w * count)
  }
  analysis <- seq_len(nrow(train))
  assessment <- setdiff(seq_len(n), analysis)
  p <- list(scale_median = 250, scale_intercept_sd = .45, scale_coef_sd = .25,
    scale_group_sd = .30, beta_sd = .70, response_group_sd = .35,
    dispersion_log_mean = if (case$family == "gaussian") log(.6) else log(2),
    dispersion_log_sd = .35)
  base <- rnorm(1, log(p$scale_median / 1000), p$scale_intercept_sd)
  trend <- if (repeated) rnorm(1, sd = p$scale_coef_sd) else 0
  scale_sd <- abs(rnorm(1, sd = p$scale_group_sd))
  log_scale <- base + trend * data$time + rnorm(n_year + 1L, sd = scale_sd)[data$year]
  beta <- rnorm(2, sd = p$beta_sd)
  response_sd <- if (repeated) abs(rnorm(1, sd = p$response_group_sd)) else 0
  re <- rnorm(sites + 1L, sd = response_sd)[data$site]
  dispersion <- exp(rnorm(1, p$dispersion_log_mean, p$dispersion_log_sd))
  ref <- means(rep(log(p$scale_median / 1000), n))[analysis]
  eta <- beta[1L] + beta[2L] * (means(log_scale) - mean(ref)) / sd(ref) + re + data$offset
  data$y <- switch(case$family, gaussian = rnorm(n, eta, dispersion),
    poisson = rpois(n, exp(eta)), binomial = rbinom(n, data$trials, plogis(eta)),
    nbinom2 = rnbinom(n, mu = exp(eta), size = dispersion))
  data$miss <- data$trials - data$y
  data$year <- factor(data$year)
  data$site <- factor(data$site)
  response <- if (case$family == "binomial") "cbind(y, miss)" else "y"
  formula <- as.formula(paste(response, "~ habitat + offset(offset)",
    if (repeated) "+ (1 | site)" else ""))
  scale <- list(habitat = if (repeated) ~ 1 + time + (1 | year) else ~ 1 + (1 | year))
  # Separate prepared inputs from training responses. Bins have fixed geometry.
  inputs <- list(kernel = kernel, min_D = 25, max_D = 1000, unit_conv = 1000,
    binned = list(bin_rep = d, vsum = list(habitat = habitat * count),
      csum = list(habitat = count), point_ids = data$key, n_points = n))
  list(data = data, analysis = analysis, assessment = assessment,
    formula = formula, scales = scale, inputs = inputs, priors = p, family = case$family,
    reference = c(mean = mean(ref), sd = sd(ref)), future_log_distance =
      log_scale[assessment[2L]] + log(1000),
    truth = c(beta = beta[2L], log_distance = log_scale[1L] + log(1000), scale_sd = scale_sd))
}
