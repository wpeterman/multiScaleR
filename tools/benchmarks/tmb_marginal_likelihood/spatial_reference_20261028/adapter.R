# Research-only adapter for the frozen private spatial prototype. No fitting or
# producer side effects occur when this file is sourced.
sr_atomic <- function(value, path) {
  stopifnot(!file.exists(path))
  temporary <- paste0(path, ".tmp")
  writeLines(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", digits = 12), temporary)
  stopifnot(file.rename(temporary, path))
}

sr_case <- function(case_id) {
  stopifnot(case_id %in% c("null", "smooth"))
  data_seed <- if (case_id == "null") 202610281L else 202610282L
  f <- bayesian_scale_fixture(n_site = 24, n_year = 3, seed = data_seed)
  grid <- expand.grid(east = seq(0, 1000, length.out = 4), north = seq(0, 1000, length.out = 6))
  f$data$east <- grid$east[as.integer(f$data$site)]
  f$data$north <- grid$north[as.integer(f$data$site)]
  args <- list(formula = stats::as.formula("count ~ habitat + temp", env = .GlobalEnv),
    inputs = f$inputs, data = f$data,
    scales = list(habitat = stats::as.formula("~1 + (1 | year)", env = .GlobalEnv)),
    family = "gaussian", row_id = "key", priors = list(scale_median = 220),
    covariate = "habitat", coordinates = c("east", "north"),
    knots = rbind(c(200, 200), c(800, 200), c(500, 800)), range = 350, field_sd_prior = .4)
  generating <- do.call(multiScaleR:::.msr_spatial_scale_prepare, args)
  p <- generating$parameters
  p$beta <- c(.5, .2, 1.5); p$scale_beta[] <- log(.22)
  p$scale_z[] <- c(-.7, 0, .7, -.8, .6, .4)
  field_sd <- if (case_id == "null") 0 else .45
  p$log_scale_sd[] <- log(c(.15, field_sd)); p$log_dispersion[] <- log(.3)
  components <- multiScaleR:::.msr_bayes_components(generating$specification, p)
  set.seed(data_seed + 100L)
  f$data$count <- rnorm(nrow(f$data), components$eta, .3)
  args$data <- f$data
  object <- do.call(multiScaleR:::.msr_spatial_scale_prepare, args)
  log_scale <- log(components$scales$habitat)
  regions <- list(west = which(f$data$east <= 334), east = which(f$data$east >= 666),
    south = which(f$data$north <= 250), north = which(f$data$north >= 750))
  list(object = object, data = f$data, data_seed = data_seed,
    input_hash = digest::digest(list(data = f$data, inputs = f$inputs, args = args, truth = p), algo = "sha256"),
    regions = regions, variables = c(paste0("log_scale_", seq_along(log_scale)),
      "year_sd", "field_sd", "east_west", "north_south"),
    truth = c(log_scale, .15, field_sd,
      mean(log_scale[regions$east]) - mean(log_scale[regions$west]),
      mean(log_scale[regions$north]) - mean(log_scale[regions$south])))
}

sr_unpack <- function(raw, parameters) {
  if (!is.numeric(raw) || length(dim(raw)) != 3L || any(dim(raw) < 1L) ||
      is.null(dimnames(raw)[[3L]]) || anyNA(dimnames(raw)[[3L]]) ||
      any(!nzchar(dimnames(raw)[[3L]])) || anyDuplicated(dimnames(raw)[[3L]]))
    stop("Original chain array dimensions/names are invalid.")
  columns <- dimnames(raw)[[3L]]
  result <- lapply(names(parameters), function(nm) {
    n <- length(parameters[[nm]])
    keys <- paste0(nm, "[", seq_len(n), "]")
    if (n == 1L && nm %in% columns && keys %in% columns)
      stop("Ambiguous original scalar parameter layout.")
    if (n == 1L && nm %in% columns) keys <- nm
    if (n && !all(keys %in% columns)) stop("Missing original parameter layout.")
    if (!n) matrix(numeric(), prod(dim(raw)[1:2]), 0L) else
      matrix(raw[, , keys, drop = FALSE], nrow = prod(dim(raw)[1:2]))
  })
  names(result) <- names(parameters)
  result
}

sr_derive <- function(draws, case) {
  sm <- case$object$specification$scales$habitat
  n <- nrow(draws$scale_beta)
  x <- draws$scale_beta[, sm$beta_index, drop = FALSE] %*% t(sm$X)
  for (term in sm$random) {
    effect <- if (is.null(term$gaussian_block))
      draws$scale_z[, term$z_index[term$index], drop = FALSE] else
      draws$scale_z[, term$z_index, drop = FALSE] %*%
        t(term$gaussian_block$design %*% term$gaussian_block$chol)
    x <- x + effect * exp(draws$log_scale_sd[, term$sd_index])
  }
  x <- x + log(case$object$specification$unit_conv)
  regions <- case$regions
  output <- cbind(x,
    year_sd = exp(draws$log_scale_sd[, sm$random$year$sd_index]),
    field_sd = exp(draws$log_scale_sd[, sm$random$.spatial_field$sd_index]),
    east_west = rowMeans(x[, regions$east, drop = FALSE]) - rowMeans(x[, regions$west, drop = FALSE]),
    north_south = rowMeans(x[, regions$north, drop = FALSE]) - rowMeans(x[, regions$south, drop = FALSE]))
  colnames(output) <- case$variables
  stopifnot(nrow(output) == n)
  output
}

sr_derived_array <- function(raw, case) {
  values <- sr_derive(sr_unpack(raw, case$object$parameters), case)
  array(values, c(dim(raw)[1:2], ncol(values)), dimnames = list(NULL, NULL, case$variables))
}

sr_event <- function(chain, threshold) {
  finite <- all(is.finite(chain)) && is.finite(threshold)
  indicator <- if (finite) 1 * (chain <= threshold) else matrix(NA_real_, nrow(chain), ncol(chain))
  constant <- if (finite) length(unique(as.vector(indicator))) == 1L else NA
  reason <- if (!finite) "nonfinite values; no iterations dropped" else if (constant)
    "constant indicator; precision unavailable" else if (ncol(chain) < 2L || nrow(chain) < 4L)
      "insufficient original chains/iterations" else ""
  metrics <- c(mcse = NA_real_, ess_mean = NA_real_, rhat = NA_real_)
  if (!nzchar(reason)) {
    metrics <- c(mcse = posterior::mcse_mean(indicator),
      ess_mean = posterior::ess_mean(indicator), rhat = posterior::rhat(indicator))
    if (any(!is.finite(metrics)) || metrics[["mcse"]] <= 0 || metrics[["ess_mean"]] <= 0)
      reason <- "nonestimable indicator precision"
  }
  list(estimate = if (finite) mean(indicator) else NA_real_, mcse = unname(metrics[["mcse"]]),
    ess_mean = unname(metrics[["ess_mean"]]), rhat = unname(metrics[["rhat"]]),
    constant = constant, reason = reason,
    precision_pass = !nzchar(reason) && metrics[["mcse"]] <= .005,
    chain_probabilities = if (finite) unname(colMeans(indicator)) else rep(NA_real_, ncol(chain)))
}

sr_functional <- function(derived, thresholds) {
  stopifnot(all(c("variable", "quantile", "threshold") %in% names(thresholds)),
    all(thresholds$variable %in% dimnames(derived)[[3L]]))
  rows <- lapply(seq_len(nrow(thresholds)), function(i) {
    event <- sr_event(derived[, , thresholds$variable[i], drop = TRUE], thresholds$threshold[i])
    cbind(thresholds[i, , drop = FALSE], as.data.frame(event[setdiff(names(event), "chain_probabilities")]),
      as.data.frame(setNames(as.list(event$chain_probabilities), paste0("chain", seq_len(dim(derived)[2L])))))
  })
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

sr_quantiles <- function(derived) {
  do.call(rbind, lapply(dimnames(derived)[[3L]], function(variable) {
    chain <- derived[, , variable, drop = TRUE]
    finite <- all(is.finite(chain))
    estimate <- if (finite) unname(stats::quantile(chain, c(.025, .5, .975))) else rep(NA_real_, 3L)
    mcse <- if (finite) unname(posterior::mcse_quantile(chain, c(.025, .5, .975))) else rep(NA_real_, 3L)
    data.frame(variable = variable, quantile = c(.025, .5, .975), threshold = estimate, quantile_mcse = mcse)
  }))
}
