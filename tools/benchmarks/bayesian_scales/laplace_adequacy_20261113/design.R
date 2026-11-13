# Frozen prior and geometry panel. The generating parameters and observation
# process are specified independently of the package likelihood and fit.
scenarios <- data.frame(scenario = c("gaussian_shared_narrow", "gaussian_shared_broad",
  "gaussian_annual_weak", "gaussian_annual_informative", "gaussian_repeated_funnel",
  "nb2_annual", "nb2_repeated_broad"),
  family = c(rep("gaussian", 5), rep("nbinom2", 2)),
  sites = c(24L, 8L, 8L, 24L, 10L, 16L, 10L),
  years = c(1L, 1L, 3L, 4L, 5L, 4L, 4L),
  kernel = c("gaussian", "exp", "gaussian", "gaussian", "exp", "gaussian", "exp"),
  beta = c(1.4, .3, .15, 1.5, 1.1, 1.1, .65),
  residual = c(.3, .8, .65, .25, .4, 6, 1.4),
  scale_sd = c(0, 0, .02, .35, .04, .30, .08),
  site_sd = c(0, 0, 0, 0, .45, 0, .55),
  scale_prior = c(.3, .9, 1, .5, .65, .5, .9),
  sd_prior = c(.3, .5, .7, .35, .9, .45, .85))
plan <- scenarios[rep(seq_len(nrow(scenarios)), each = 5L), ]
plan$replicate <- rep(seq_len(5L), nrow(scenarios))
plan$case <- sprintf("%s-%02d", plan$scenario, plan$replicate)
plan$seed <- 1130400L + seq_len(nrow(plan)) * 100L

simulate_case <- function(case) {
  set.seed(case$seed)
  d <- expand.grid(site = seq_len(case$sites), year = seq_len(case$years))
  if (case$site_sd > 0) d <- d[(d$site + d$year) %% 6L != 0L, , drop = FALSE]
  n <- nrow(d)
  d$key <- sprintf("point_%03d", seq_len(n))
  d$time <- (d$year - 1) / max(1L, case$years - 1L)
  d$offset <- runif(n, -.25, .25)
  distance <- seq(.025, 1, length.out = 14L)
  counts <- matrix(sample(2:6, n * length(distance), replace = TRUE), n)
  values <- outer(rnorm(n), sin(7 * distance)) + outer(rnorm(n), cos(15 * distance)) +
    matrix(rnorm(n * length(distance), sd = .18), n)
  power <- if (case$kernel == "gaussian") 2L else 1L
  kernel <- function(scale) {
    weights <- exp(-outer(1 / scale^power, distance^power - min(distance)^power) / power)
    rowSums(weights * values * counts) / rowSums(weights * counts)
  }
  true_scale <- exp(log(.25) + .2 * d$time +
    rnorm(case$years, sd = case$scale_sd)[d$year])
  reference <- kernel(rep(.25, n))
  x <- (kernel(true_scale) - mean(reference)) / sd(reference)
  eta <- (if (case$family == "gaussian") .4 else log(2)) + case$beta * x +
    rnorm(case$sites, sd = case$site_sd)[d$site] + d$offset
  d$y <- if (case$family == "gaussian") rnorm(n, eta, case$residual) else
    rnbinom(n, mu = exp(eta), size = case$residual)
  d$site <- factor(d$site)
  d$year <- factor(d$year)
  response <- if (case$site_sd > 0)
    y ~ habitat + offset(offset) + (1 | site) else y ~ habitat + offset(offset)
  scales <- list(habitat = if (case$years == 1L) ~ 1 else if (case$site_sd > 0)
    ~ 1 + time + (1 | year) else ~ 1 + (1 | year))
  inputs <- list(kernel = case$kernel, min_D = 25, max_D = 1000, unit_conv = 1000,
    binned = list(bin_rep = distance, vsum = list(habitat = values * counts),
      csum = list(habitat = counts), point_ids = d$key, n_points = n))
  priors <- list(scale_median = 250, scale_intercept_sd = case$scale_prior,
    scale_coef_sd = .4, scale_group_sd = case$sd_prior, beta_sd = 1.25,
    response_group_sd = .75, dispersion_log_mean = log(case$residual), dispersion_log_sd = .45)
  list(data = d, inputs = inputs, formula = response, scales = scales, priors = priors,
    reference = c(mean = mean(reference), sd = sd(reference)))
}
