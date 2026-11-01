# No fitted/sampled outcomes are selected here. These are deterministic software
# fixtures with independent cell-weight response oracles and synthetic vectors.
spatial_response_fixture <- function(kernel = "gaussian", family = "gaussian",
    metric = "mean", groups = FALSE) {
  f <- bayesian_scale_fixture(n_site = 8, n_year = if (groups) 2 else 1,
    annual = FALSE, kernel = kernel)
  if (!groups) f$data$year <- NULL
  angle <- 2 * pi * (as.integer(f$data$site) - 1) / 8
  f$data$east <- cos(angle); f$data$north <- sin(angle)
  if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
  if (metric == "sq") {
    d <- f$inputs$binned$bin_rep; ids <- f$inputs$binned$point_ids
    f$inputs$d_list <- rep(list(d), length(ids))
    f$inputs$raw_cov <- lapply(seq_along(ids), function(i)
      data.frame(elev = sin(d * (i + 1)) + i * d^2))
    f$inputs$raw_cov[[1L]]$elev[] <- 2
    f$inputs$kernel_dat <- data.frame(row.names = ids)
    f$inputs$scale_vars <- msr_vars(habitat = surface_var("elev", "sq", weighted = TRUE))
    f$inputs$binned <- NULL
  }
  b <- .msr_spatial_scale_prepare(
    if (groups) count ~ habitat + temp + offset(log(effort)) + (1 | site) else
      count ~ habitat + temp + offset(log(effort)), f$inputs, f$data,
    scales = list(habitat = if (groups) ~1 + (1 | year) else ~1),
    family = family, row_id = "key", priors = list(scale_median = 200),
    covariate = "habitat", coordinates = c("east", "north"),
    knots = matrix(c(-.5, 0, .5, 0, .5, -.5), 3, 2), range = .7, field_sd_prior = .4)
  p <- b$parameters; p$beta <- c(.2, -.1, .5)
  p$scale_z[] <- seq(-.2, .4, length.out = length(p$scale_z))
  p$log_scale_sd[] <- log(.35)
  if (groups) p$response_z[] <- seq(-.3, .3, length.out = length(p$response_z))
  list(fixture = f, object = b, parameters = p)
}

# Calculate weights directly from the supplied cells/bins, independently of the
# new response context and evaluator. The existing scale-only adapter supplies
# the separately tested frozen spatial field, in physical distance units.
spatial_response_oracle <- function(x, data, parameters = x$parameters) {
  inputs <- x$fixture$inputs; spec <- x$object$specification
  sigma <- .msr_spatial_scale_predict(x$object, data, parameters = parameters) / inputs$unit_conv
  ids <- .msr_scale_point_ids(inputs); rows <- match(data$key, ids)
  raw <- vapply(seq_len(nrow(data)), function(i) {
    d <- if (is.null(inputs$binned)) inputs$d_list[[rows[i]]] else inputs$binned$bin_rep
    w <- if (inputs$kernel == "gaussian") exp(-.5 * (d / sigma[i])^2) else exp(-d / sigma[i])
    if (is.null(inputs$binned)) {
      v <- inputs$raw_cov[[rows[i]]]$elev
      sqrt(sum(w * (v - sum(w * v) / sum(w))^2) / sum(w))
    } else sum(w * inputs$binned$vsum$habitat[rows[i], ]) /
      sum(w * inputs$binned$csum$habitat[rows[i], ])
  }, 0)
  parameters$beta[1] + parameters$beta[2] * data$temp + log(data$effort) +
    parameters$beta[3] * (raw - spec$center[["habitat"]]) / spec$spread[["habitat"]]
}

