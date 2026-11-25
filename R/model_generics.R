#' Standard model methods for multiScaleR fits
#'
#' @description
#' The usual model accessors work on a fit from [multiScale_optim()]. Each
#' reads the final model, the model refit at the estimated scales of effect
#' and stored in `object$opt_mod`. `coef(fit)`, `fitted(fit)`, and the other
#' methods therefore describe that model, never the starting model passed to
#' [multiScale_optim()]. Step-selection wrappers from `amt::fit_clogit()` and
#' `amt::fit_issf()` are first unwrapped to their native `survival` model.
#'
#' A fit from [multiScale_bayes()] has no single fitted model, so most of
#' these methods stop with a pointer to the Bayesian tools (see "Bayesian
#' fits").
#'
#' @param object,x A `multiScaleR` object from [multiScale_optim()], or for
#'   `coef()`, `fitted()`, `residuals()`, `vcov()`, `logLik()`, and `nobs()`
#'   a `multiScaleR_bayes` object from [multiScale_bayes()].
#' @param ... Passed to the native method of the final model. Ignored by the
#'   Bayesian methods.
#'
#' @return
#' * `coef()`, `fitted()`, `residuals()`, `nobs()`, `formula()`, `family()`,
#'   `deviance()`, and `df.residual()`: the value the same generic returns for
#'   the final model.
#' * `logLik()`: an object of class `"logLik"` holding the log-likelihood of
#'   the final model, with a `df` attribute that also counts the estimated
#'   kernel parameters and a `nobs` attribute (see "Log-likelihood and
#'   information criteria").
#' * `nobs()` for a Bayesian fit: the number of observations.
#'
#' @details
#' # Methods that delegate to the final model
#'
#' `coef()`, `fitted()`, `residuals()`, `nobs()`, `formula()`, `family()`,
#' `deviance()`, and `df.residual()` call the same generic on the final
#' model, so their results, and their errors, are the native ones. For
#' example, `unmarked` fits support `coef()`, `fitted()`, and `residuals()`
#' through the `unmarked` methods, while `family()` stops with the native
#' error. Coefficient standard errors from `summary(fit$opt_mod)` or
#' `vcov(fit$opt_mod)` condition on the estimated scales; see
#' [vcov.multiScaleR()] for the joint covariance of the analytic backends.
#'
#' # Log-likelihood and information criteria
#'
#' `logLik()` returns the maximized log-likelihood of the final model. Its
#' `df` attribute counts the native model parameters plus each estimated
#' scale (sigma) and each estimated kernel shape, exactly as `K` in
#' [aic_tab()]. `AIC(fit)` therefore equals the AIC reported by
#' `aic_tab(list(fit), AICc = FALSE)`, and `BIC(fit)` equals the BIC from
#' [bic_tab()]. Both are larger than `AIC(fit$opt_mod)` and
#' `BIC(fit$opt_mod)`, which leave the kernel parameters out. The `nobs`
#' attribute is the sample size those tables use, for example the event
#' count for Cox and conditional logistic models.
#'
#' # Bayesian fits
#'
#' `coef()`, `fitted()`, `residuals()`, `vcov()`, and `logLik()` stop for a
#' `multiScaleR_bayes` object, because a posterior has no single set of
#' estimates or maximized likelihood. Use `summary()` for posterior summaries
#' of the coefficients and scales, [posterior_scales()] for the scale draws,
#' and `predict()` for fitted values with their posterior uncertainty.
#' `nobs()` returns the number of observations in the fitted data.
#'
#' @seealso [multiScale_optim()], [summary.multiScaleR()], [aic_tab()],
#'   [bic_tab()], [vcov.multiScaleR()], [multiScaleR-class].
#'
#' @examples
#' data(pts)
#' data(count_data)
#' hab <- terra::rast(system.file("extdata", "hab.tif",
#'                                package = "multiScaleR"))
#' inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
#'                       kernel = "gaussian", verbose = FALSE)
#' fit <- multiScale_optim(glm(y ~ hab, family = poisson, data = count_data),
#'                         inputs, verbose = FALSE)
#' coef(fit)              # coefficients of the model at the estimated scale
#' head(fitted(fit))      # fitted values of that model
#' logLik(fit)            # df counts the estimated sigma
#' AIC(fit)               # equals aic_tab(list(fit), AICc = FALSE)$AIC
#' AIC(fit$opt_mod)       # native AIC, without the sigma parameter
#' @name multiScaleR-methods
NULL

# The final model of a multiScale_optim() fit, with amt wrappers unwrapped.
.msr_final_model <- function(object) {
  model <- object$opt_mod
  if (is.null(model))
    stop("`object$opt_mod` is missing; the fit holds no final model.", call. = FALSE)
  .analysis_model(model)
}

# The native generic for a final model. unmarked defines S4 methods for some
# stats generics; call those through the unmarked generic so S4 dispatch
# reaches them. Other classes use the stats generic and its native error.
.msr_native_generic <- function(name, model) {
  if (methods::is(model, "unmarkedFit") &&
      requireNamespace("unmarked", quietly = TRUE) &&
      name %in% getNamespaceExports("unmarked")) {
    return(getExportedValue("unmarked", name))
  }
  getExportedValue("stats", name)
}

# Call one generic on the final model of a fit.
.msr_delegate <- function(name, object, ...) {
  model <- .msr_final_model(object)
  .msr_native_generic(name, model)(model, ...)
}

#' @rdname multiScaleR-methods
#' @exportS3Method stats::coef
coef.multiScaleR <- function(object, ...) .msr_delegate("coef", object, ...)

#' @rdname multiScaleR-methods
#' @exportS3Method stats::fitted
fitted.multiScaleR <- function(object, ...) .msr_delegate("fitted", object, ...)

#' @rdname multiScaleR-methods
#' @exportS3Method stats::residuals
residuals.multiScaleR <- function(object, ...) .msr_delegate("residuals", object, ...)

#' @rdname multiScaleR-methods
#' @exportS3Method stats::nobs
nobs.multiScaleR <- function(object, ...) .msr_delegate("nobs", object, ...)

#' @rdname multiScaleR-methods
#' @exportS3Method stats::formula
formula.multiScaleR <- function(x, ...) .msr_delegate("formula", x, ...)

#' @rdname multiScaleR-methods
#' @exportS3Method stats::family
family.multiScaleR <- function(object, ...) .msr_delegate("family", object, ...)

#' @rdname multiScaleR-methods
#' @exportS3Method stats::deviance
deviance.multiScaleR <- function(object, ...) .msr_delegate("deviance", object, ...)

#' @rdname multiScaleR-methods
#' @exportS3Method stats::df.residual
df.residual.multiScaleR <- function(object, ...) .msr_delegate("df.residual", object, ...)

#' @rdname multiScaleR-methods
#' @exportS3Method stats::logLik
logLik.multiScaleR <- function(object, ...) {
  # Follow aic_tab() and bic_tab() step by step so that AIC(), BIC(), and the
  # tables agree: the same likelihood (including the checked negative
  # binomial BAM adjustment), the same parameter count with each estimated
  # scale and shape added, and the same sample size.
  model <- .msr_selection_model(object)
  checked <- .msr_bam_nb_selection(model)
  value <- if (is.null(checked)) .msr_model_loglik(model) else checked$loglik
  k <- .msr_selection_parameter_count(object, if (is.null(checked)) NULL else checked$K)
  structure(value, df = k, nobs = .msr_selection_nobs(model), class = "logLik")
}

# Bayesian fits ------------------------------------------------------------

# A posterior has no single set of estimates; name the tools that summarize it.
.msr_bayes_generic_stop <- function(what) {
  stop(what, " is not defined for a Bayesian multiScaleR fit, which has no ",
       "single set of estimates. Use summary() for posterior summaries, ",
       "posterior_scales() for the scale draws, and predict() for fitted ",
       "values with posterior uncertainty.", call. = FALSE)
}

#' @rdname multiScaleR-methods
#' @exportS3Method stats::coef
coef.multiScaleR_bayes <- function(object, ...) .msr_bayes_generic_stop("coef()")

#' @rdname multiScaleR-methods
#' @exportS3Method stats::fitted
fitted.multiScaleR_bayes <- function(object, ...) .msr_bayes_generic_stop("fitted()")

#' @rdname multiScaleR-methods
#' @exportS3Method stats::residuals
residuals.multiScaleR_bayes <- function(object, ...) .msr_bayes_generic_stop("residuals()")

#' @rdname multiScaleR-methods
#' @exportS3Method stats::vcov
vcov.multiScaleR_bayes <- function(object, ...) .msr_bayes_generic_stop("vcov()")

#' @rdname multiScaleR-methods
#' @exportS3Method stats::logLik
logLik.multiScaleR_bayes <- function(object, ...) {
  stop("logLik() is not defined for a Bayesian multiScaleR fit: a posterior ",
       "is not a maximized likelihood. Use summary(), posterior_scales(), ",
       "and predict() to describe the posterior.", call. = FALSE)
}

#' @rdname multiScaleR-methods
#' @exportS3Method stats::nobs
nobs.multiScaleR_bayes <- function(object, ...) {
  # summary() and print() of a Bayesian fit report the same count.
  length(object$specification$y)
}
