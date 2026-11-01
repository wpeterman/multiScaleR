# Existing Bayesian numerical inputs only; no recovery/sampling experiment.
prepared_formula_fixture <- function(kernel = "gaussian", surface = FALSE,
                                     two = TRUE) {
  f <- bayesian_scale_fixture(n_site = 8, n_year = 3, annual = FALSE, kernel = kernel)
  b <- f$inputs$binned
  if (surface) {
    f$inputs$d_list <- rep(list(b$bin_rep), length(b$point_ids))
    f$inputs$raw_cov <- lapply(seq_along(b$point_ids), function(i)
      data.frame(elev = sin(b$bin_rep * (i + 1)) + i * b$bin_rep^2))
    f$inputs$raw_cov[[1]]$elev[] <- 2
    f$inputs$raw_cov[[2]]$elev[2] <- NA_real_
    f$inputs$kernel_dat <- data.frame(row.names = b$point_ids)
    f$inputs$scale_vars <- msr_vars(habitat = kernel_var("habitat"),
      other = surface_var("elev", "sq", weighted = TRUE))
  } else if (two) {
    f$inputs$binned$vsum$other <- b$vsum$habitat[, rev(seq_along(b$bin_rep))]
    f$inputs$binned$csum$other <- b$csum$habitat[, rev(seq_along(b$bin_rep))]
  }
  f$covariates <- if (two) c("habitat", "other") else "habitat"
  f$formula <- if (two) count ~ habitat * other + habitat * year + poly(temp, 2) +
    offset(log(effort)) else count ~ habitat * year + poly(temp, 2) + offset(log(effort))
  f$scales <- setNames(rep(list(~1 + (1 | year)), length(f$covariates)), f$covariates)
  f
}

# Native terms carry the fitted polynomial/contrast metadata. Reconstructing
# from the original formula would relearn a basis and would be a different test.
prepared_formula_native_eta <- function(spec, parameters, data) {
  components <- multiScaleR:::.msr_bayes_components(spec, parameters, include_response = FALSE)
  for (nm in spec$covariates) data[[nm]] <- components$kernel[[nm]]
  tt <- spec$response$terms; attr(tt, ".Environment") <- asNamespace("stats")
  frame <- model.frame(tt, data, na.action = na.fail,
    xlev = spec$response$prediction$xlevels)
  X <- model.matrix(tt, frame, contrasts.arg = spec$response$prediction$contrasts)
  offset <- model.offset(frame)
  if (is.null(offset)) offset <- rep(0, nrow(X))
  list(eta = drop(X %*% parameters$beta) + offset, X = X)
}

prepared_formula_prior <- function(spec, p) {
  pr <- spec$priors
  value <- sum(dnorm(p$beta, 0, pr$beta_sd, log = TRUE)) +
    sum(dnorm(p$scale_z, log = TRUE)) + sum(dnorm(p$response_z, log = TRUE)) +
    sum(dnorm(exp(p$log_scale_sd), 0, pr$scale_group_sd, log = TRUE) + log(2) + p$log_scale_sd) +
    sum(dnorm(exp(p$log_response_sd), 0, pr$response_group_sd, log = TRUE) + log(2) + p$log_response_sd) +
    sum(dnorm(p$log_dispersion, pr$dispersion_log_mean, pr$dispersion_log_sd, log = TRUE))
  for (sm in spec$scales) value <- value + sum(dnorm(p$scale_beta[sm$beta_index],
    sm$prior_mean, sm$prior_sd, log = TRUE))
  value
}
