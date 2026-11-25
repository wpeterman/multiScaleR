.msr_selection_family_label <- function(model) {
  model <- .analysis_model(model)
  family <- tryCatch(stats::family(model)$family, error = function(e) NULL)
  if (is.character(family) && length(family) == 1L) return(family)
  if (inherits(model, c("lm", "gls", "lme"))) return("gaussian")
  if (inherits(model, "coxph")) return("partial likelihood")
  if (methods::is(model, "unmarkedFit")) {
    label <- class(model)[1L]
    if ("mixture" %in% methods::slotNames(model)) {
      mixture <- methods::slot(model, "mixture")
      if (is.character(mixture) && length(mixture) == 1L && !is.na(mixture))
        label <- paste0(label, " (mixture=", mixture, ")")
    }
    return(label)
  }
  "unknown"
}

.msr_selection_model <- function(model) {
  retained <- if (inherits(model, "multiScaleR")) model$opt_mod else model
  if (methods::is(retained, "unmarkedFitGDRField")) .msr_gdr_field_guard(retained)
  if (inherits(model, "multiScaleR_cv"))
    stop("CV results contain held-out scores, not a maximized likelihood; use cv_tab() or select_cv() with retained training fits.", call. = FALSE)
  if (inherits(model, "multiScaleR_bayes"))
    stop("AIC/BIC require maximized likelihoods; Bayesian posterior fits need an appropriate posterior predictive comparison, not a posterior objective substituted for logLik.", call. = FALSE)
  if (inherits(model, "msr_spatial_scale_prototype"))
    stop("Private spatial prototypes have no qualified maximized marginal likelihood for AIC/BIC; conditional response densities and posterior objectives are not substitutes.", call. = FALSE)
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


# Recover original site identities before dropping sites, including null models.
.msr_unmarked_observations <- function(model) {
  n <- unmarked::numSites(model@data)
  used <- unmarked::sampleSize(model)
  valid_n <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= 0 && x <= .Machine$integer.max && x == floor(x)
  if (!valid_n(n) || !valid_n(used) || used > n)
    stop("Invalid native unmarked site-count metadata.", call. = FALSE)
  omitted <- model@sitesRemoved
  if (!is.numeric(omitted) || anyNA(omitted) || any(!is.finite(omitted)) ||
      any(omitted != floor(omitted) | omitted < 1 | omitted > n) || anyDuplicated(omitted))
    stop("Invalid native unmarked omitted-site metadata.", call. = FALSE)
  rows <- setdiff(seq_len(n), omitted)
  if (length(rows) != used)
    stop("Native unmarked retained-site count disagrees with omitted-site metadata.", call. = FALSE)
  ids <- rownames(model@data@siteCovs)
  if (!.msr_valid_observation_ids(ids, n)) ids <- rownames(unmarked::getY(model@data))
  if (!.msr_valid_observation_ids(ids, n)) ids <- as.character(seq_len(n))
  list(rows = rows, ids = as.character(ids[rows]))
}


.msr_model_nobs <- function(model) {
  analysis_mod <- .analysis_model(model)

  if (methods::is(analysis_mod, "unmarkedFit"))
    return(as.integer(length(.msr_unmarked_observations(analysis_mod)$rows)))

  if (inherits(analysis_mod, "glm")) return(as.integer(stats::nobs(analysis_mod)))

  as.integer(n_obs(analysis_mod)[1])
}


.msr_selection_nobs <- function(model) {
  analysis_mod <- .analysis_model(.msr_selection_model(model))
  if (!inherits(analysis_mod, "coxph")) return(.msr_model_nobs(analysis_mod))

  # survival's likelihood metadata counts events, while row matching must
  # retain every alternative in each conditional choice set.
  n <- tryCatch(attr(stats::logLik(analysis_mod), "nobs"),
                error = function(e) NULL)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 1 || n > .Machine$integer.max || n != floor(n)) {
    stop("Cox information criteria require a positive integer native logLik event count; alternative rows cannot substitute.",
         call. = FALSE)
  }
  as.integer(n)
}


.msr_model_loglik <- function(model) {
  analysis_mod <- .analysis_model(.msr_selection_model(model))
  .msr_gdr_native_guard(analysis_mod)
  checked <- .msr_bam_nb_selection(analysis_mod)
  if (!is.null(checked)) return(checked$loglik)
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

  if (inherits(analysis_mod, "HLfit")) {
    if (!requireNamespace("spaMM", quietly = TRUE))
      stop("spaMM is required to recover HLfit likelihood parameter counts.", call. = FALSE)
    # DoF is the provider's marginal-IC count, including fitted dispersion,
    # covariance and correlation terms. extractAIC's edf counts fixed terms only.
    k <- sum(spaMM::DoF(analysis_mod))
    return(.msr_validate_parameter_count(k, analysis_mod, "spaMM marginal parameter metadata"))
  }

  checked <- .msr_bam_nb_selection(analysis_mod)
  if (!is.null(checked)) return(checked$K)

  # logLik degrees of freedom include ancillary likelihood parameters such as
  # negative-binomial theta, which coefficient tables may omit.
  loglik_df <- tryCatch(
    attr(stats::logLik(analysis_mod), "df"),
    error = function(e) NULL
  )
  if (!is.null(loglik_df)) {
    # GAM likelihoods use effective degrees of freedom, which need not be integer.
    if (inherits(analysis_mod, "gam") && is.numeric(loglik_df) &&
        length(loglik_df) == 1L && is.finite(loglik_df) && loglik_df >= 0)
      return(as.numeric(loglik_df))
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


.msr_selection_parameter_count <- function(model, native_k = NULL) {
  if (!inherits(model, "multiScaleR")) {
    return(if (is.null(native_k)) .msr_parameter_count(model) else native_k)
  }

  extra_k <- 0L
  if (!is.null(model$scale_est)) {
    extra_k <- extra_k + nrow(model$scale_est)
  }
  if (!is.null(model$shape_est)) {
    hard <- .msr_hard_radius_covariates(model$kernel_inputs$scale_vars,
      model$kernel_inputs$kernel, rownames(model$shape_est))
    extra_k <- extra_k + sum(!rownames(model$shape_est) %in% hard)
  }

  if (is.null(native_k)) native_k <- .msr_parameter_count(model$opt_mod)
  native_k + extra_k
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

  if (methods::is(analysis_mod, "unmarkedFit"))
    return(.msr_unmarked_observations(analysis_mod)$ids)

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


# Check the shape of `mod_list` before any model is inspected, so an empty
# list, a bare fit, or a NULL element gets a message about `mod_list` itself.
# Bayesian fits and CV results pass here and keep their own refusals in
# .msr_selection_model().
.msr_validate_mod_list <- function(mod_list, fn) {
  if (!is.list(mod_list) || is.object(mod_list))
    stop(sprintf(paste0("`mod_list` must be a plain list of fitted models, for ",
      "example %s(list(fit1, fit2)); wrap a single fit in list()."), fn),
      call. = FALSE)
  if (!length(mod_list))
    stop(sprintf(paste0("`mod_list` is empty. Supply a list of one or more ",
      "fitted models, for example %s(list(fit1, fit2))."), fn), call. = FALSE)
  not_model <- vapply(mod_list, function(m) is.null(m) || is.atomic(m) ||
    is.function(m) || inherits(m, "formula"), logical(1))
  if (any(not_model)) {
    labels <- names(mod_list)
    if (is.null(labels)) labels <- rep("", length(mod_list))
    labels <- ifelse(nzchar(labels), paste0("`", labels, "`"),
                     paste0("element ", seq_along(mod_list)))
    what <- vapply(mod_list[not_model], function(m)
      if (is.null(m)) "NULL" else paste0("a `", class(m)[1L], "`"), character(1))
    stop(paste0("Every `mod_list` element must be a fitted model. Not a model: ",
      paste0(labels[not_model], " (", what, ")", collapse = ", "),
      ". Remove these elements or replace them with fitted models."),
      call. = FALSE)
  }
  invisible(TRUE)
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
#' @title Rank models by AIC or AICc
#'
#' @description Build a model-selection table ranked by AICc (the default) or
#' AIC for a list of fitted models. Each estimated scale of effect (sigma), and
#' each estimated kernel shape, counts as a parameter, so a scale-optimized
#' model competes fairly with fixed-scale or non-spatial alternatives in the
#' same list. Lists can mix `multiScaleR` objects with plain fitted models
#' such as `glm` objects.
#'
#' `vignette("v02-user-guide", package = "multiScaleR")` compares kernels and
#' model structures with these tables, including an N-mixture model.
#'
#' @param mod_list A list of fitted models: `multiScaleR` objects from
#'   [multiScale_optim()], plain fitted models (for example `glm` or `lm`), or
#'   both. Every model must be fitted to the same observations with a
#'   comparable likelihood (see "Check that models are comparable"). Bayesian
#'   fits, cross-validation results, and private spatial prototypes are
#'   rejected.
#' @param AICc Logical. If `TRUE` (the default), rank by the second-order
#'   corrected AIC (AICc). Set `FALSE` for standard AIC. AICc is recommended
#'   when the number of observations is small relative to the number of
#'   parameters (n / K < 40 is a common rule of thumb). AICc requires
#'   n > K + 1 for every model and raises an error otherwise. Use
#'   `AICc = FALSE` in that case.
#' @param mod_names Optional character vector of model names, one per model
#'   in `mod_list`. By default, `multiScaleR` objects are labeled
#'   `[kernel=...; family=...]` followed by the right-hand side of the model
#'   formula. `unmarked` labels name the model class and, where present, the
#'   mixture, so Poisson, negative binomial, and ZIP alternatives get distinct
#'   labels.
#' @param verbose Logical. If `TRUE`, the table is also printed to the
#'   console. Default `FALSE`.
#' @param ... Not used.
#'
#' @return A data frame of class `"aictab"` (from the AICcmodavg package),
#' sorted from the best (smallest) criterion value, with one row per model:
#' * `Modnames`: model name, from `mod_names` or the default label.
#' * `K`: number of estimated parameters. It counts the fitted likelihood
#'   parameters, including ancillary ones such as negative binomial
#'   dispersion (or native effective degrees of freedom for GAMs), plus each
#'   estimated scale and shape of a `multiScaleR` model.
#' * `AICc` or `AIC`: the information criterion. Smaller is better.
#' * `Delta_AICc` or `Delta_AIC`: difference from the best model. Models
#'   within about 2 units of the best are generally considered to have
#'   comparable support.
#' * `ModelLik`: relative likelihood of the model, `exp(-0.5 * Delta)`.
#' * `AICcWt` or `AICWt`: Akaike weight, the relative support for each model
#'   within this set. Weights sum to 1.
#' * `LL`: maximized log-likelihood.
#' * `Cum.Wt`: cumulative Akaike weight.
#'
#' When a likelihood was adjusted (see "Count parameters"), the table carries
#' a `likelihood_adjustments` attribute.
#'
#' @export
#'
#' @details
#' # Count parameters
#'
#' For a `multiScaleR` object, the log-likelihood and parameter count come
#' from the refit model in `opt_mod`, and each estimated sigma (and shape, for
#' `kernel = "expow"`) is added to `K`. The AIC printed by the native model in
#' `opt_mod` leaves those parameters out. Hard-radius metrics in an `"expow"`
#' model add their estimated radius but no shape, and fixed-radius metrics add
#' neither. Plain models in the list get no sigma. Otherwise the native
#' log-likelihood degrees of freedom are used when available, so negative
#' binomial dispersion counts as one parameter. Missing or inconsistent
#' parameter metadata and non-finite log-likelihoods raise an error rather
#' than guessing from predictor names.
#'
#' * **`unmarked` models.** The native AIC and log-likelihood determine the
#'   count. It includes every abundance or occupancy coefficient, every
#'   detection coefficient, factor contrasts, mixture parameters, and
#'   estimated random-effect standard deviations. Optimizer and TMB
#'   free-parameter metadata are checked for agreement when present. Latent
#'   abundances and realized random effects add nothing. For example, four
#'   abundance coefficients, three detection coefficients, and one
#'   random-effect standard deviation give eight native parameters, or nine
#'   with one estimated Gaussian scale. A table of fixed coefficients can
#'   therefore show fewer parameters than the fitted count.
#' * **The `K` argument of `unmarked::pcount()`** is an upper bound for the
#'   abundance sum in the likelihood, not a parameter count. Check that the
#'   log-likelihood does not change when you refit with a larger bound.
#' * **spaMM `HLfit` models.** `spaMM::DoF()` supplies the full count,
#'   including estimated residual, covariance, and correlation parameters.
#'   Realized random effects are not counted. The native default
#'   log-likelihood target is kept (`p_v` for ordinary maximum likelihood,
#'   `p_bv` for REML), so absolute REML criteria can differ from spaMM's
#'   marginal AIC, which uses `p_v`.
#' * **mgcv negative binomial BAMs.** The table checks the native fit against
#'   the normalized count likelihood. A retained outer scale flag can
#'   otherwise add one to both the reported log-likelihood and the parameter
#'   count, even though the family scale is one. Verified fits use the
#'   normalized likelihood and native effective degrees of freedom plus an
#'   estimated theta, with no extra scale parameter. A fixed theta adds no
#'   parameter. The `likelihood_adjustments` attribute records the native and
#'   used likelihoods and regression-model counts for the affected models,
#'   before the scale and shape parameters are added. With large case
#'   weights, verification checks the native AIC arithmetic and a rescaled
#'   response-density comparison at the original tolerance, and records that
#'   method. Reported likelihoods keep the original weights. Unknown
#'   inconsistencies raise an error, and fitted coefficients, covariance, and
#'   diagnostics are unchanged.
#'
#' # Check that models are comparable
#'
#' Information criteria compare models only when they describe the same
#' observations with comparable likelihoods. The table checks the number of
#' observations and, when native metadata can be recovered, the used row IDs,
#' response values, binomial trials, likelihood weights, and the likelihood's
#' response measure. Zero-weight GLM rows and removed `unmarked` sites are
#' excluded from the used observations. Binomial counts and proportions with
#' their trial weights are compared on the same count scale, and logical
#' `FALSE`/`TRUE` outcomes, including success and failure columns, compare as
#' 0 and 1 when trial totals and weights agree. Fitted responses, covariance,
#' and parameter counts are not changed.
#'
#' * **N-mixture models.** Mixed `unmarked::pcount()` models keep their
#'   native marginal likelihood. Two models can share sites and the original
#'   count matrix yet use different visits when a detection covariate has
#'   missing values. The table recovers the visits each model used from its
#'   component formulas, including formula-specific missing observations and
#'   omitted sites, and stops when they differ. Equal count matrices and site
#'   totals alone do not establish comparability, and missing or inconsistent
#'   used-response metadata raises an error without falling back on raw
#'   data. Null abundance models and detection-only models without site
#'   covariates keep the native site count and their original response-row
#'   identities.
#' * **Distance-removal models.** `gdistremoval()` comparisons check both
#'   distance and removal histories and their stored survey protocol,
#'   including the primary-period layout, distance bins and units, and
#'   removal durations. Equal paired totals are not enough. Distances in
#'   meters and kilometers are converted to meters in the comparison copy,
#'   and the fitted objects keep their native units. Harmonize time units and
#'   event definitions before comparing, because a native frame does not
#'   store a removal-time unit. Mixture and distance-key alternatives may
#'   share the same observed target. Check sensitivity to the integration
#'   cutoff separately.
#' * **Maximum likelihood and REML.** ML and REML fits cannot be mixed, and
#'   REML comparisons require the same recoverable fixed-effect design.
#' * **Cox and conditional logistic models.** Partial-likelihood comparisons
#'   require matching methods and equivalent stratum partitions.
#'
#' These checks cannot certify metadata that a native provider does not
#' store. You still need to establish that the likelihoods are
#' scientifically comparable.
#'
#' # Sample size for AICc and BIC
#'
#' * **`unmarked` models.** AICc and BIC use the retained site rows, not visits or
#'   random-effect realizations. When groups repeat across rows, this
#'   convention does not establish a valid small-sample correction or an
#'   effective independent sample size.
#' * **Cox and conditional logistic models.** AICc and BIC use the event count from
#'   the native log-likelihood metadata, so each supported choice set with
#'   one used step contributes one event, while every alternative row is
#'   still checked for comparability. This follows the native survival
#'   BIC convention. It does not adjust for dependence among steps or animals
#'   or establish that AICc is valid for a particular sampling design. Missing
#'   event-count metadata raises an error.
#'
#' A finite criterion and an optimizer code of zero do not establish
#' convergence, identifiability, or interval calibration. Read the native
#' fit diagnostics, including the gradient and curvature when available.
#'
#' # Models the table rejects
#'
#' * **Bayesian fits** from [multiScale_bayes()]. A posterior mode or
#'   penalized target is not a maximized sampling likelihood for AIC or BIC.
#'   Posterior predictive comparison needs a separate workflow with its own
#'   target and diagnostics.
#' * **Cross-validation results**, which hold held-out scores rather than a
#'   maximized likelihood. Use [cv_tab()] or [select_cv()].
#' * **Private spatial prototypes**, which have no qualified maximized
#'   marginal likelihood. Conditional response densities and posterior
#'   objectives are not substitutes.
#' * **`gdistremoval()` models with active formula offsets or TMB zero
#'   inflation.** These are rejected before fitting, model selection, or
#'   prediction. In the audited `unmarked` 1.5.2 provider, offsets did not
#'   enter the fitted likelihood, and the TMB ZIP likelihood ignored its
#'   zero-inflation parameter. These combinations need separate provider
#'   qualification. multiScaleR does not remove offsets, switch engines,
#'   alter parameter counts, or substitute uncertainty estimates. Offset-free
#'   C-engine ZIP models keep their native likelihood and parameter count.
#'
#' The table itself is built with [AICcmodavg::aictabCustom()].
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
#' @seealso [bic_tab()], [multiScale_optim()], [cv_tab()] for comparisons by
#'   held-out prediction, and
#'   `vignette("v08-cross-validation", package = "multiScaleR")`.
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
#' ## The fixed model already uses the true 75-unit scale, so the spatial model
#' ## reaches the same log-likelihood and pays for one extra parameter (sigma).
#' aic_tab(models, AICc = FALSE, mod_names = names(models))
#' aic_tab(models, AICc = TRUE, mod_names = names(models))
#'
#' ## K for a stored N-mixture fit: three abundance coefficients, one
#' ## detection coefficient, and two estimated scales
#' nmix <- readRDS(system.file("extdata", "cache-user-guide-nmix.rds",
#'                             package = "multiScaleR"))
#' aic_tab(list(nmix), mod_names = "bin1 + cont2")
aic_tab <- function(mod_list,
                    AICc = TRUE,
                    mod_names = NULL,
                    verbose = FALSE,
                    ...){

  p <- list(...)

  .msr_validate_mod_list(mod_list, "aic_tab")
  opt_list <- lapply(mod_list, .msr_selection_model)
  .msr_validate_comparable_models(opt_list)
  mod_n <- vapply(opt_list, .msr_selection_nobs, integer(1))
  mod_kernel <- vapply(mod_list, .msr_kernel_label, character(1))
  likelihood <- lapply(opt_list, .msr_bam_nb_selection)
  k <- vapply(seq_along(mod_list), function(i)
    .msr_selection_parameter_count(mod_list[[i]], if (is.null(likelihood[[i]])) NULL else likelihood[[i]]$K), numeric(1))

  if(is.null(mod_names)){
    mod_eq <- vapply(opt_list, .msr_selection_formula_label, character(1))
    mod_names <- paste0("[kernel=", mod_kernel, "; family=",
      vapply(opt_list, .msr_selection_family_label, character(1)), "]", mod_eq)
  }

  mod_loglik <- vapply(seq_along(opt_list), function(i) if (is.null(likelihood[[i]]))
    .msr_model_loglik(opt_list[[i]]) else likelihood[[i]]$loglik, numeric(1))
  if (isTRUE(AICc) && any(mod_n <= k + 1L)) {
    stop("AICc requires n > K + 1 for every model; use `AICc = FALSE` for AIC.",
         call. = FALSE)
  }

  # AIC ---------------------------------------------------------------------
  if(AICc == FALSE) {

    tab <- aictabCustom(logL = mod_loglik,
                        K = k,
                        modnames = mod_names,
                        second.ord = FALSE,
                        nobs = mod_n,
                        sort = TRUE)

    adjustments <- .msr_selection_adjustments(likelihood)
    if (length(adjustments)) attr(tab, "likelihood_adjustments") <- adjustments
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
                        nobs = mod_n,
                        sort = TRUE)

    adjustments <- .msr_selection_adjustments(likelihood)
    if (length(adjustments)) attr(tab, "likelihood_adjustments") <- adjustments
    if(isTRUE(verbose)){
      return(print(tab))
    } else {
      return(tab)
    }
  }
} ## End function


# BIC table ---------------------------------------------------------------
#' @author Bill Peterman
#' @title Rank models by BIC
#'
#' @description Build a model-selection table ranked by BIC (the Bayesian
#' information criterion) for a list of fitted models. As in [aic_tab()], each
#' estimated scale of effect (sigma) and kernel shape counts as a parameter,
#' and lists can mix `multiScaleR` objects with plain fitted models.
#'
#' @param mod_list A list of fitted models: `multiScaleR` objects from
#'   [multiScale_optim()], plain fitted models, or both. Every model must be
#'   fitted to the same observations with a comparable likelihood, including
#'   matching responses, weights, and likelihood measure. Bayesian fits,
#'   cross-validation results, and private spatial prototypes are rejected.
#' @param mod_names Optional character vector of model names, one per model
#'   in `mod_list`. By default, `multiScaleR` objects are labeled
#'   `[kernel=...; family=...]` followed by the right-hand side of the model
#'   formula.
#' @param verbose Logical. If `TRUE`, the table is also printed to the
#'   console. Default `FALSE`.
#' @param ... Not used.
#'
#' @return A data frame of class `"bictab"` (from the AICcmodavg package),
#' sorted from the best (smallest) BIC, with one row per model:
#' * `Modnames`: model name.
#' * `K`: number of estimated parameters, counted as in [aic_tab()]: the
#'   fitted likelihood parameters (including ancillary ones such as negative
#'   binomial dispersion, or native effective degrees of freedom for GAMs)
#'   plus each estimated scale and shape of a `multiScaleR` model.
#' * `BIC`: the Bayesian information criterion. Smaller is better.
#' * `Delta_BIC`: difference from the best model.
#' * `BICWt`: BIC weight, the relative support for each model within this
#'   set. Weights sum to 1.
#' * `Cum.Wt`: cumulative BIC weight.
#' * `LL`: maximized log-likelihood.
#'
#' A `likelihood_adjustments` attribute is added when a likelihood was
#' adjusted, as described in [aic_tab()].
#'
#' @details
#' BIC penalizes each parameter by the log of the sample size rather than by
#' 2, so it penalizes complexity more heavily than AIC once there are more
#' than about 7 observations. It is often preferred when the goal is the
#' single best-supported model rather than model averaging.
#'
#' Every rule in [aic_tab()] applies here: the parameter count (including all
#' estimated components of `unmarked` models, spaMM degrees of freedom, and
#' the checked negative binomial BAM adjustment), the comparability checks
#' (responses, rows, weights, binomial outcomes, ML and REML, partial
#' likelihoods, used `pcount()` visits, and distance-removal protocols), the
#' sample-size conventions (retained site rows for `unmarked` models and
#' native event counts for Cox and conditional logistic models), and the
#' rejected models. Inapplicable hard-radius shape positions are excluded
#' from `K`. Calling the criterion BIC does not make a posterior target a
#' sampling likelihood, so Bayesian fits are rejected here too.
#'
#' The table itself is built with [AICcmodavg::bictabCustom()].
#'
#' @export
#'
#' @usage
#' bic_tab(mod_list,
#'         mod_names = NULL,
#'         verbose = FALSE,
#'         ...)
#' @rdname bic_tab
#'
#' @seealso [aic_tab()], [multiScale_optim()], and
#'   `vignette("v02-user-guide", package = "multiScaleR")`.
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
#'
#' ## BIC charges log(75), about 4.3, for the estimated sigma, versus 2 in AIC
#' bic_tab(list(spatial_model, fixed_model, null_model),
#'         mod_names = c("spatial", "fixed", "null"))
bic_tab <- function(mod_list,
                    mod_names = NULL,
                    verbose = FALSE,
                    ...){

  p <- list(...)

  .msr_validate_mod_list(mod_list, "bic_tab")
  opt_list <- lapply(mod_list, .msr_selection_model)
  .msr_validate_comparable_models(opt_list)
  mod_n <- vapply(opt_list, .msr_selection_nobs, integer(1))
  mod_kernel <- vapply(mod_list, .msr_kernel_label, character(1))
  likelihood <- lapply(opt_list, .msr_bam_nb_selection)
  k <- vapply(seq_along(mod_list), function(i)
    .msr_selection_parameter_count(mod_list[[i]], if (is.null(likelihood[[i]])) NULL else likelihood[[i]]$K), numeric(1))

  if(is.null(mod_names)){
    mod_eq <- vapply(opt_list, .msr_selection_formula_label, character(1))
    mod_names <- paste0("[kernel=", mod_kernel, "; family=",
      vapply(opt_list, .msr_selection_family_label, character(1)), "]", mod_eq)
  }

  mod_loglik <- vapply(seq_along(opt_list), function(i) if (is.null(likelihood[[i]]))
    .msr_model_loglik(opt_list[[i]]) else likelihood[[i]]$loglik, numeric(1))
  # * BIC ---------------------------------------------------------------------

  tab <- bictabCustom(logL = mod_loglik,
                      K = k,
                      modnames = mod_names,
                      nobs = mod_n,
                      sort = TRUE)

  adjustments <- .msr_selection_adjustments(likelihood)
  if (length(adjustments)) attr(tab, "likelihood_adjustments") <- adjustments
  if(isTRUE(verbose)){
    return(print(tab))
  } else {
    return(tab)
  }
}  ## End function
