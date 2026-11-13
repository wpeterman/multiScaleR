# Frozen design for actual posterior refits and matched predictive references.
plan <- expand.grid(replicate = seq_len(8L), family = c("gaussian", "nbinom2"),
  stringsAsFactors = FALSE)
plan$generator <- ifelse(plan$replicate %% 2L == 0L, "annual", "pooled")
plan$seed <- 863000L + 30L * seq_len(nrow(plan))
plan$case <- sprintf("%s-%02d", plan$family, plan$replicate)

simulate_cv_case <- function(case) {
  set.seed(case$seed)
  annual <- case$generator == "annual"
  gaussian <- case$family == "gaussian"
  point <- expand.grid(site = paste0("s", 1:7), year = paste0("y", 1:5),
    stringsAsFactors = FALSE)
  point$key <- paste(point$site, point$year, sep = "_")
  point <- point[point$year != "y5" | point$site %in% paste0("s", 1:4), ]
  data <- point[rep(seq_len(nrow(point)), each = 2L), ]
  rownames(data) <- paste0("visit", seq_len(nrow(data)))
  data$offset <- runif(nrow(data), -.15, .15)
  d <- seq(.035, .965, length.out = 12L)
  nbins <- 12L
  values <- outer(rnorm(nrow(point)), cos(9 * d)) +
    outer(rnorm(nrow(point)), sin(11 * d)) +
    outer(rnorm(nrow(point)), 1 + .2 * sin(4 * d)) +
    matrix(rnorm(nrow(point) * length(d), sd = .25), nrow(point))
  # One declared geometric distance in each bin. Training representative
  # distances therefore have an exact reference independent of fold size.
  raw <- setNames(lapply(seq_len(nrow(point)),
    function(i) cbind(cover = values[i, ])), point$key)
  dist <- setNames(rep(list(d), nrow(point)), point$key)
  inputs <- list(kernel = "gaussian", unit_conv = 1000, min_D = 35,
    max_D = 1000, kernel_dat = data.frame(habitat = 999, row.names = NULL),
    d_list = dist, raw_cov = raw, scale_vars = msr_vars(habitat = kernel_var("cover")),
    binned = list(point_ids = point$key, nbins = nbins, bin_rep = d,
      vsum = list(habitat = values), csum = list(habitat = matrix(1, nrow(point), nbins))))
  inputs$kernel_dat <- data.frame(habitat = rep(999, nrow(point)), row.names = point$key)
  rownames(values) <- point$key
  folds <- list(whole_year = list(analysis = rownames(data)[data$year != "y5"],
      assessment = rownames(data)[data$year == "y5"]),
    seen_year_new_site = list(analysis = rownames(data)[data$site != "s7"],
      assessment = rownames(data)[data$site == "s7"]))
  truth <- list(log_scale = log(.26), annual_sd = if (annual) .55 else 0,
    site_sd = if (gaussian) .55 else .65, dispersion = if (gaussian) .60 else 2.2,
    beta0 = if (gaussian) .1 else 1.0, beta1 = if (gaussian) 1.15 else .70)
  year <- setNames(rnorm(5L, truth$log_scale, truth$annual_sd), paste0("y", 1:5))
  site <- setNames(rnorm(7L, sd = truth$site_sd), paste0("s", 1:7))
  w <- exp(-.5 * (d / .26)^2)
  standard <- drop(values %*% (w / sum(w)))
  truth$center <- mean(standard)
  truth$spread <- sd(standard)
  weights <- exp(-.5 * outer(exp(-2 * year[point$year]), d^2))
  k <- setNames(rowSums(weights * values) / rowSums(weights), point$key)
  eta <- truth$beta0 + truth$beta1 * (k[data$key] - truth$center) / truth$spread +
    site[data$site] + data$offset
  data$y <- if (gaussian) rnorm(nrow(data), eta, truth$dispersion) else
    rnbinom(nrow(data), mu = exp(eta), size = truth$dispersion)
  data$site <- factor(data$site)
  data$year <- factor(data$year)
  priors <- list(scale_median = 260, scale_intercept_sd = .5, scale_group_sd = .55,
    beta_sd = 1.5, response_group_sd = .7,
    dispersion_log_mean = log(truth$dispersion), dispersion_log_sd = .5)
  list(data = data, point = point, values = values, distances = d, folds = folds,
    inputs = inputs, truth = truth, realized_year = year, realized_site = site,
    priors = priors, family = case$family)
}
