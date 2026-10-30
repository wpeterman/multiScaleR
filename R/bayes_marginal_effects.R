# Posterior regression curves at hypothetical covariate values. Spatial
# covariates use the retained fixed reference transformation, not a new scale
# or landscape extraction per draw. These are not spatial prediction maps.
.msr_bayes_marginal_effects <- function(object, ylab, length.out, draws) {
  x <- NULL # ggplot resolves the x column in its supplied data frame.
  validate_character_scalar(ylab, "ylab")
  validate_scalar_numeric(length.out, "length.out", integerish = TRUE, positive = TRUE)
  all_ids <- .msr_bayes_draw_ids(object, NULL)
  ids <- .msr_bayes_draw_ids(object, if (is.null(draws)) min(1000L, length(all_ids)) else draws)
  spec <- object$specification
  beta <- object$draws$beta[ids, , drop = FALSE]
  if (ncol(beta) != length(spec$beta_names) || any(!is.finite(beta)))
    stop("Stored finite response coefficients must match the fitted response design.", call. = FALSE)
  formula <- stats::formula(spec$response$terms)
  rhs <- formula[[length(formula)]]
  variables <- .msr_marginal_variables(rhs, offsets = TRUE)
  panels <- .msr_marginal_variables(rhs, offsets = FALSE)
  reference <- spec$data
  for (nm in spec$covariates) reference[[nm]] <- spec$response$frame[[nm]]
  if (!all(variables %in% names(reference)))
    stop("The fit must retain base covariates needed for effect curves.", call. = FALSE)
  training <- cbind(spec$X_fixed, as.matrix(reference[spec$covariates]))
  defaults <- lapply(reference[variables], function(value) {
    if (all(is.na(value))) stop("Effect covariates must contain observed values.", call. = FALSE)
    if (is.numeric(value)) mean(value, na.rm = TRUE) else value[which(!is.na(value))[1L]]
  })
  inverse <- switch(spec$family, gaussian = identity,
    poisson = exp, nbinom2 = exp, binomial = stats::plogis)
  caption <- paste0("Central 95% credible intervals from ", length(ids),
    " stored coefficient draws. Hypothetical covariate values; spatial scales are not reevaluated.\n",
    "Response random effects are zero, not integrated. Observation noise is excluded.")
  if (isFALSE(object$diagnostics$passed)) caption <- paste(caption, "Fit checks need attention.")
  else if (!isTRUE(object$diagnostics$passed)) caption <- paste(caption, "Stored fit checks unavailable.")
  plots <- stats::setNames(lapply(panels, function(nm) {
    value <- reference[[nm]]
    grid <- if (is.numeric(value)) seq(min(value, na.rm = TRUE), max(value, na.rm = TRUE),
      length.out = length.out) else if (is.factor(value))
      value[match(levels(droplevels(value)), as.character(value))] else unique(value[!is.na(value)])
    nd <- as.data.frame(lapply(defaults, function(x) rep(x, length(grid))),
                         stringsAsFactors = FALSE)
    nd[[nm]] <- grid
    design <- .msr_bayes_prediction_matrix(spec$response, nd, training, response = TRUE)
    offset <- stats::model.offset(design$frame)
    if (is.null(offset)) offset <- rep(0, length(grid))
    eta <- sweep(design$X %*% t(beta), 1L, offset, "+")
    value <- inverse(eta)
    if (any(!is.finite(value)))
      stop("Nonfinite posterior effect curves; inspect the retained fit and covariate range.", call. = FALSE)
    quantile <- t(apply(value, 1L, stats::quantile, probs = c(.025, .5, .975), names = FALSE))
    axis <- if (nm %in% spec$covariates) grid * spec$spread[[nm]] + spec$center[[nm]] else
      if (nm %in% spec$fixed_reference$covariates)
        grid * spec$fixed_reference$spread[[nm]] + spec$fixed_reference$center[[nm]] else grid
    df <- data.frame(x = axis, fit = quantile[, 2L], lwr = quantile[, 1L],
                      upr = quantile[, 3L], variable = nm)
    plot <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = fit)) +
      ggplot2::labs(x = nm, y = ylab, caption = caption) + cowplot::theme_cowplot()
    if (is.numeric(axis)) plot <- plot + ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = .25) else
      plot <- plot + ggplot2::geom_point() +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lwr, ymax = upr), width = .15)
    plot
  }), panels)
  attr(plots, "prediction_contract") <- .msr_prediction_contract("posterior",
    "posterior_covariate_slice", "expected_response", "response",
    population = "stored_marginal_response_coefficient_draws",
    kernel = "covariate_values_fixed_not_reevaluated",
    effects = list(response = "zero_deviation_reference"),
    conditioning = c("hypothetical_covariate_values", "retained_reference_normalization",
      "other_numeric_covariates_at_reference_means", "other_factors_at_first_observed_level"),
    intervals = "central_credible_coefficient_curve",
    omitted = c("landscape_reevaluation_under_kernel_draws", "response_random_effect_integration",
      "observation_noise"), coordinates = "stored_response_coefficients",
    diagnostics = object$diagnostics, fixed_inputs = TRUE)
  attr(plots, "draw_ids") <- ids
  lapply(plots, print)
  invisible(plots)
}
