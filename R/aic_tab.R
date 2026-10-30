.msr_selection_model <- function(model) {
  if (inherits(model, "multiScaleR_bayes"))
    stop("AIC/BIC require maximized likelihoods; Bayesian posterior fits need an appropriate posterior predictive comparison, not a posterior objective substituted for logLik.", call. = FALSE)
  if (inherits(model, "multiScaleR")) {
    return(model$opt_mod)
  }

  model
}


.msr_kernel_label <- function(model) {
  if (inherits(model, "multiScaleR")) {
    return(model$kernel_inputs$kernel)
  }

  "NA"
}


.msr_model_nobs <- function(model) {
  analysis_mod <- .analysis_model(model)

  if (any(grepl("^unmarked", class(analysis_mod)))) {
    return(dim(analysis_mod@data@siteCovs)[1] - length(analysis_mod@sitesRemoved))
  }

  if (inherits(analysis_mod, "glm")) return(as.integer(stats::nobs(analysis_mod)))

  as.integer(n_obs(analysis_mod)[1])
}


.msr_model_loglik <- function(model) {
  analysis_mod <- .analysis_model(.msr_selection_model(model))
  value <- tryCatch(
    if (methods::is(analysis_mod, "unmarkedFit")) {
      unmarked::logLik(analysis_mod)
    } else {
      logLik(analysis_mod)
    },
    error = function(e) NULL
  )
  if (is.null(value) && !methods::is(analysis_mod, "unmarkedFit")) {
    value <- tryCatch(suppressWarnings(insight::get_loglikelihood(analysis_mod)),
                      error = function(e) NULL)
  }
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop("Could not extract a finite scalar log-likelihood from model class `",
         class(analysis_mod)[1], "`.", call. = FALSE)
  }
  as.numeric(value)
}


.msr_validate_parameter_count <- function(k, model, source) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k < 0 || k > .Machine$integer.max || abs(k - round(k)) > 1e-7) {
    stop("Could not determine fitted parameter count for model class `",
         class(model)[1], "`: invalid ", source, ".", call. = FALSE)
  }
  as.integer(round(k))
}


.msr_unmarked_parameter_count <- function(model) {
  ll <- .msr_model_loglik(model)
  aic <- model@AIC
  if (!is.numeric(aic) || length(aic) != 1L || !is.finite(aic)) {
    stop("Could not determine fitted parameter count for model class `",
         class(model)[1], "`: missing or malformed native AIC metadata.",
         call. = FALSE)
  }
  # Native AIC counts all modeled components. TMB's optimizer vector contains
  # fixed/hyperparameters, excluding integrated random-effect realizations.
  k <- .msr_validate_parameter_count((aic + 2 * ll) / 2, model,
                                      "count implied by native AIC and log-likelihood")
  par <- model@opt$par
  if (!is.null(par)) {
    if (!is.numeric(par) || !is.null(dim(par)) || any(!is.finite(par))) {
      stop("Could not determine fitted parameter count for model class `",
           class(model)[1], "`: malformed optimizer parameter metadata.",
           call. = FALSE)
    }
    if (length(par) != k) {
      stop("Inconsistent native AIC and optimizer parameter count for model class `",
           class(model)[1], "`; free fitted parameters cannot be inferred reliably.",
           call. = FALSE)
    }
  }
  if ("TMB" %in% methods::slotNames(model) && !is.null(model@TMB)) {
    tmb_par <- model@TMB$par
    if (!is.numeric(tmb_par) || !is.null(dim(tmb_par)) ||
        any(!is.finite(tmb_par)) || length(tmb_par) != k) {
      stop("Inconsistent native AIC and TMB free parameter metadata for model class `",
           class(model)[1], "` (including fixed or mapped parameters).", call. = FALSE)
    }
  }
  k
}


.msr_parameter_count <- function(model) {
  analysis_mod <- .analysis_model(model)

  if (methods::is(analysis_mod, "unmarkedFit")) {
    return(.msr_unmarked_parameter_count(analysis_mod))
  }

  # logLik degrees of freedom include ancillary likelihood parameters such as
  # negative-binomial theta, which coefficient tables may omit.
  loglik_df <- tryCatch(
    attr(stats::logLik(analysis_mod), "df"),
    error = function(e) NULL
  )
  if (!is.null(loglik_df)) {
    return(.msr_validate_parameter_count(loglik_df, analysis_mod, "logLik degrees of freedom"))
  }

  params <- tryCatch(
    suppressWarnings(get_parameters(analysis_mod)),
    error = function(e) NULL
  )
  if (is.data.frame(params) && nrow(params) > 0) {
    return(nrow(params))
  }

  coefs <- tryCatch(
    stats::coef(analysis_mod),
    error = function(e) NULL
  )
  if (length(coefs) > 0) {
    return(length(coefs))
  }

  stop("Could not determine the number of fitted parameters for model class `",
       class(analysis_mod)[1], "`: no reliable parameter metadata available.",
       call. = FALSE)
}


.msr_selection_formula_label <- function(model) {
  model <- .analysis_model(model)
  if (methods::is(model, "unmarkedFit")) {
    if ("formula" %in% methods::slotNames(model)) {
      return(paste(deparse(model@formula), collapse = " "))
    }
    return(paste(vapply(model@formlist, function(f) paste(deparse(f), collapse = " "),
                        character(1)), collapse = "; "))
  }
  paste((find_formula(model)$conditional)[-2], collapse = " ")
}


.msr_selection_parameter_count <- function(model) {
  if (!inherits(model, "multiScaleR")) {
    return(.msr_parameter_count(model))
  }

  extra_k <- 0L
  if (!is.null(model$scale_est)) {
    extra_k <- extra_k + nrow(model$scale_est)
  }
  if (!is.null(model$shape_est)) {
    hard <- .msr_hard_radius_covariates(model$kernel_inputs$scale_vars)
    extra_k <- extra_k + sum(!rownames(model$shape_est) %in% hard)
  }

  .msr_parameter_count(model$opt_mod) + extra_k
}


.msr_valid_observation_ids <- function(ids, n_obs) {
  length(ids) == n_obs &&
    !anyNA(ids) &&
    all(nzchar(ids)) &&
    !anyDuplicated(ids)
}


.msr_observation_ids <- function(model) {
  analysis_mod <- .analysis_model(model)
  mod_n <- .msr_model_nobs(model)

  if (any(grepl("^unmarked", class(analysis_mod)))) {
    dat <- analysis_mod@data@siteCovs
    if (length(analysis_mod@sitesRemoved))
      dat <- dat[-analysis_mod@sitesRemoved, , drop = FALSE]
    ids <- row.names(dat)
    if (.msr_valid_observation_ids(ids, mod_n)) {
      return(as.character(ids))
    }
    return(as.character(seq_len(mod_n)))
  }

  frame <- .msr_selection_frame(analysis_mod)
  if (!is.null(frame)) {
    if (inherits(analysis_mod, "glm") && !is.null(analysis_mod$prior.weights))
      frame <- frame[analysis_mod$prior.weights != 0, , drop = FALSE]
    ids <- rownames(frame)
    if (.msr_valid_observation_ids(ids, mod_n)) return(as.character(ids))
  }

  dat <- .model_data(analysis_mod, effects = "all")
  if (is.null(dat)) {
    dat <- extract_model_data(model)
  }

  if (is.data.frame(dat)) {
    ids <- row.names(dat)
    if (.msr_valid_observation_ids(ids, mod_n)) {
      return(as.character(ids))
    }
  }

  residual_ids <- tryCatch(
    names(stats::residuals(analysis_mod)),
    error = function(e) NULL
  )
  if (.msr_valid_observation_ids(residual_ids, mod_n)) {
    return(as.character(residual_ids))
  }

  fitted_ids <- tryCatch(
    names(stats::fitted(analysis_mod)),
    error = function(e) NULL
  )
  if (.msr_valid_observation_ids(fitted_ids, mod_n)) {
    return(as.character(fitted_ids))
  }

  NULL
}


.msr_validate_comparable_models <- function(mod_list) {
  mod_dims <- vapply(mod_list, .msr_model_nobs, integer(1))
  if (length(unique.default(mod_dims)) != 1L) {
    stop("\nYou are attempting to compare models with different number of sample locations. These are not valid comparisons.\n")
  }

  obs_ids <- lapply(mod_list, .msr_observation_ids)
  if (all(vapply(obs_ids, Negate(is.null), logical(1)))) {
    ref_ids <- sort(obs_ids[[1]])
    same_obs <- vapply(
      obs_ids[-1],
      function(x) identical(sort(x), ref_ids),
      logical(1)
    )
    if (length(same_obs) > 0 && !all(same_obs)) {
      stop("\nYou are attempting to compare models fit to different observation sets. These are not valid comparisons.\n")
    }
  }

  .msr_assert_selection_targets(mod_list)
  mod_dims
}

# AIC table ---------------------------------------------------------------
#' @author Bill Peterman
#' @title multiScaleR model selection
#'
#' @description Creates a model selection table ranked by AIC or AICc for a
#' list of fitted models. Mixed lists containing both \code{multiScaleR}
#' objects and plain model objects (e.g., \code{glm}) are supported, allowing
#' comparison of scale-optimized models against fixed-scale alternatives.
#'
#' @param mod_list List containing fitted \code{multiScaleR} objects (from
#'   \code{\link{multiScale_optim}}) and/or plain fitted model objects (e.g.,
#'   \code{glm}, \code{lm}). All models must have been fit to the same set of
#'   observations; an error is raised if sample sizes differ.
#' @param AICc Logical. If \code{TRUE} (default), the second-order corrected
#'   AIC (AICc) is used. Set to \code{FALSE} for standard AIC. AICc is
#'   recommended when the ratio of observations to parameters is small
#'   (n / K < 40 is a common rule of thumb).
#' @param mod_names Optional character vector of model names. Length must equal
#'   \code{length(mod_list)}. By default, names are constructed as
#'   \code{[kernel]right-hand-side-formula} for \code{multiScaleR} objects.
#' @param verbose Logical. If \code{TRUE}, the table is also printed to the
#'   console. Default: \code{FALSE}.
#' @param ... Additional arguments (not used).
#'
#' @return A data frame of class \code{"aictab"} (from the \code{AICcmodavg}
#' package) sorted by ascending AIC(c), containing:
#' \describe{
#'   \item{\code{Modnames}}{Model name (from \code{mod_names} or auto-generated).}
#'   \item{\code{K}}{Number of fitted likelihood parameters, including
#'     ancillary parameters such as negative-binomial dispersion, plus optimized
#'     scale and shape parameters for \code{multiScaleR} models.}
#'   \item{\code{AICc} or \code{AIC}}{Information criterion value.}
#'   \item{\code{Delta_AICc} or \code{Delta_AIC}}{Difference from the top model.}
#'   \item{\code{ModelLik}}{Relative likelihood (\code{exp(-0.5 * Delta)}).}
#'   \item{\code{AICcWt} or \code{AICWt}}{Akaike weight.}
#'   \item{\code{LL}}{Log-likelihood.}
#'   \item{\code{Cum.Wt}}{Cumulative Akaike weight.}
#' }
#'
#' @export
#'
#' @details
#' The table is constructed using \code{\link[AICcmodavg]{aictabCustom}} from
#' the \code{AICcmodavg} package. For \code{multiScaleR} objects, the optimized
#' (refitted) model stored in \code{opt_mod} is used for log-likelihood and
#' parameter count. The sigma parameter (and shape for \code{"expow"}) is added
#' to K automatically. For plain model objects in a mixed list, sigma is not
#' added to K. The fitted model's log-likelihood degrees of freedom are used
#' when available, so negative-binomial dispersion contributes one parameter.
#' For \code{unmarkedFit} models, native AIC and log-likelihood metadata
#' determine the count, including detection coefficients, factor contrasts,
#' mixture parameters, and estimated random-effect variances. Optimizer and
#' TMB free-parameter metadata are checked for agreement when present. Latent
#' random effects and the integration limit \code{K} are not counted.
#' Missing or inconsistent parameter metadata and non-finite log-likelihoods
#' raise an error instead of guessing from predictor names.
#' Hard-radius metrics in native expow mixtures contribute their estimated
#' radius but no inapplicable shape parameter. Fixed-radius metrics contribute
#' neither an estimated radius nor shape.
#'
#' When native metadata can be recovered, comparison also checks used row IDs,
#' response values, binomial trials, likelihood weights, and the likelihood's
#' response measure. Binomial counts and proportions with their trial weights
#' are compared on the same count representation. Zero-weight GLM rows and
#' removed unmarked sites are excluded from the used-observation count.
#' ML and REML fits cannot be mixed; REML comparisons require the same
#' recoverable fixed-effect design. Cox partial-likelihood comparisons require
#' matching methods and equivalent stratum partitions. These checks cannot
#' certify metadata unavailable from a native provider, and the user must
#' still establish scientific comparability of the likelihoods.
#'
#' Bayesian posterior fits are rejected. A posterior mode or penalized target
#' is not a maximized sampling likelihood for AIC or BIC. Posterior predictive
#' comparison requires a separate workflow with its own target and diagnostics.
#' AICc tables raise an error if any model has \code{n <= K + 1}; use
#' \code{AICc = FALSE} for an AIC table in that case.
#'
#' @usage
#' aic_tab(mod_list,
#'         AICc = TRUE,
#'         mod_names = NULL,
#'         verbose = FALSE,
#'         ...)
#' @rdname aic_tab
#' @importFrom insight find_formula get_loglikelihood get_parameters n_obs
#' @importFrom AICcmodavg aictabCustom bictabCustom
#'
#' @examples
#' ## Bundled counts simulated with a habitat scale of 75 map units
#' data("pts")
#' data("count_data")
#' hab <- terra::rast(system.file("extdata", "hab.tif", package = "multiScaleR"))
#' prepared <- kernel_prep(pts, hab, max_D = 300, sigma = 75, verbose = FALSE)
#' dat <- data.frame(y = count_data$y, prepared$kernel_dat)
#' fixed_model <- glm(y ~ hab, family = poisson(), data = dat)
#' null_model <- glm(y ~ 1, family = poisson(), data = dat)
#' spatial_model <- multiScale_optim(fixed_model, prepared, par = 75 / 300,
#'                                   n_cores = 1, verbose = FALSE)
#' models <- list(spatial = spatial_model, fixed = fixed_model, null = null_model)
#'
#' ## The spatial model pays for one extra estimated sigma parameter.
#' aic_tab(models, AICc = FALSE, mod_names = names(models))
#' aic_tab(models, AICc = TRUE, mod_names = names(models))
aic_tab <- function(mod_list,
                    AICc = TRUE,
                    mod_names = NULL,
                    verbose = FALSE,
                    ...){

  p <- list(...)

  opt_list <- lapply(mod_list, .msr_selection_model)
  mod_df <- .msr_validate_comparable_models(opt_list)
  mod_kernel <- vapply(mod_list, .msr_kernel_label, character(1))
  k <- vapply(mod_list, .msr_selection_parameter_count, integer(1))

  if(is.null(mod_names)){
    mod_eq <- vapply(opt_list, .msr_selection_formula_label, character(1))
    mod_names <- paste0("[" ,mod_kernel, "]",mod_eq)
  }

  mod_loglik <- vapply(opt_list, .msr_model_loglik, numeric(1))
  if (isTRUE(AICc) && any(mod_df <= k + 1L)) {
    stop("AICc requires n > K + 1 for every model; use `AICc = FALSE` for AIC.",
         call. = FALSE)
  }

  # AIC ---------------------------------------------------------------------
  if(AICc == FALSE) {

    tab <- aictabCustom(logL = mod_loglik,
                        K = k,
                        modnames = mod_names,
                        second.ord = FALSE,
                        nobs = mod_df,
                        sort = TRUE)

    if(isTRUE(verbose)){
      return(print(tab))
    } else {
      return(tab)
    }

    # AICc --------------------------------------------------------------------

  } else {

    tab <- aictabCustom(logL = mod_loglik,
                        K = k,
                        modnames = mod_names,
                        second.ord = TRUE,
                        nobs = mod_df,
                        sort = TRUE)

    if(isTRUE(verbose)){
      return(print(tab))
    } else {
      return(tab)
    }
  }
} ## End function


# BIC table ---------------------------------------------------------------
#' @author Bill Peterman
#' @title multiScaleR model selection
#'
#' @description Creates a model selection table ranked by BIC for a list of
#' fitted models. Mixed lists containing both \code{multiScaleR} objects and
#' plain model objects are supported.
#'
#' @param mod_list List containing fitted \code{multiScaleR} objects (from
#'   \code{\link{multiScale_optim}}) and/or plain fitted model objects. All
#'   models must have been fit to the same observations and comparable
#'   likelihoods, including matching responses, weights, and likelihood measure.
#' @param mod_names Optional character vector of model names. Length must equal
#'   \code{length(mod_list)}. By default, names are constructed as
#'   \code{[kernel]right-hand-side-formula} for \code{multiScaleR} objects.
#' @param verbose Logical. If \code{TRUE}, the table is also printed to the
#'   console. Default: \code{FALSE}.
#' @param ... Additional arguments (not used).
#'
#' @return A data frame of class \code{"bictab"} (from the \code{AICcmodavg}
#' package) sorted by ascending BIC, containing:
#' \describe{
#'   \item{\code{Modnames}}{Model name.}
#'   \item{\code{K}}{Number of fitted likelihood parameters, including
#'     ancillary parameters such as negative-binomial dispersion, plus optimized
#'     scale and shape parameters for \code{multiScaleR} models.}
#'   \item{\code{BIC}}{Bayesian information criterion value.}
#'   \item{\code{Delta_BIC}}{Difference from the top model.}
#'   \item{\code{BICWt}}{BIC weight.}
#'   \item{\code{Cum.Wt}}{Cumulative BIC weight.}
#'   \item{\code{LL}}{Log-likelihood.}
#' }
#'
#' @details
#' The table is constructed using \code{\link[AICcmodavg]{bictabCustom}} from
#' the \code{AICcmodavg} package. BIC penalizes model complexity more heavily
#' than AIC as sample size grows, and is often preferred when the goal is
#' identifying the single best-supported model rather than model averaging.
#' Sigma (and shape for \code{"expow"}) is counted in K.
#' Fitted parameter counts and log-likelihood validation follow
#' \code{\link{aic_tab}}, including all estimated components of unmarked models.
#' The same recoverable response, row, weight, ML/REML, and partial-likelihood
#' checks apply. Inapplicable hard-radius shape positions are excluded.
#' Bayesian posterior fits are rejected; the name BIC does not make a posterior
#' target a sampling likelihood.
#' @export
#'
#' @usage
#' bic_tab(mod_list,
#'         mod_names = NULL,
#'         verbose = FALSE,
#'         ...)
#' @rdname bic_tab
#'
#' @examples
#' ## Simulate data
#' set.seed(555)
#'
#' points <- vect(cbind(c(5,7,9,11,13),
#'                      c(13,11,9,7,5)))
#'
#' mat_list <- list(r1 = rast(matrix(rnorm(20^2),
#'                                   nrow = 20)),
#'                  r2 = rast(matrix(rnorm(20^2),
#'                                   nrow = 20)))
#' rast_stack <- rast(mat_list)
#' kernel_inputs <- kernel_prep(pts = points,
#'                              raster_stack = rast_stack,
#'                              max_D = 5,
#'                              kernel = 'gaussian',
#'                              sigma = NULL)
#' ## Example response data
#' y <- rnorm(5)
#'
#' ## Create data frame with raster variables
#' dat <- data.frame(y = y,
#'                   kernel_inputs$kernel_dat)
#' mod1 <- glm(y ~ r1,
#'             data = dat)
#' mod2 <- glm(y ~ r2,
#'             data = dat)
#' mod3 <- glm(y ~ r1 + r2,
#'             data = dat)
#'
#' ## NOTE: This code is only for demonstration
#' ## Optimization results will have no meaning
#'
#' opt_mod1 <- multiScale_optim(fitted_mod = mod1,
#'                              kernel_inputs = kernel_inputs,
#'                              par = NULL,
#'                              n_cores = NULL)
#'
#' opt_mod2 <- multiScale_optim(fitted_mod = mod2,
#'                              kernel_inputs = kernel_inputs,
#'                              par = NULL,
#'                              n_cores = NULL)
#' opt_mod3 <- multiScale_optim(fitted_mod = mod3,
#'                              kernel_inputs = kernel_inputs,
#'                              par = NULL,
#'                              n_cores = NULL)
#'
#'  ## BIC table
#'  mod_list <- list(opt_mod1, opt_mod2, opt_mod3)
#'
#'  bic_tab(mod_list = mod_list)
#'
#'  ## BIC table with specified names
#'  bic_tab(mod_list = mod_list,
#'          mod_names = c('mod1', 'mod2', 'mod3'))

bic_tab <- function(mod_list,
                    mod_names = NULL,
                    verbose = FALSE,
                    ...){

  p <- list(...)

  opt_list <- lapply(mod_list, .msr_selection_model)
  mod_df <- .msr_validate_comparable_models(opt_list)
  mod_kernel <- vapply(mod_list, .msr_kernel_label, character(1))
  k <- vapply(mod_list, .msr_selection_parameter_count, integer(1))

  if(is.null(mod_names)){
    mod_eq <- vapply(opt_list, .msr_selection_formula_label, character(1))
    mod_names <- paste0("[" ,mod_kernel, "]",mod_eq)
  }

  mod_loglik <- vapply(opt_list, .msr_model_loglik, numeric(1))
  # * BIC ---------------------------------------------------------------------

  tab <- bictabCustom(logL = mod_loglik,
                      K = k,
                      modnames = mod_names,
                      nobs = mod_df,
                      sort = TRUE)

  if(isTRUE(verbose)){
    return(print(tab))
  } else {
    return(tab)
  }
}  ## End function
