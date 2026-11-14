# Resolve the explicit spatial refit transform before any fold is fitted.
.msr_cv_covariate_transform <- function(optim_args, backend, spatial = TRUE) {
  selected <- optim_args$covariate_transform
  if (is.null(selected)) return("candidate")
  selected <- match.arg(selected, c("candidate", "raw"))
  if (selected == "raw" && (!isTRUE(spatial) || backend != "refit"))
    stop("Raw optimized covariates require spatial CV with backend = 'refit'.", call. = FALSE)
  selected
}

.msr_cv_raw_identity <- function(center, spread, transform, scale_vars) {
  if (identical(transform, "raw")) {
    optimized <- if (is.null(scale_vars)) names(center) else
      scale_vars$covariate[scale_vars$optimize]
    optimized <- intersect(names(center), optimized)
    center[optimized] <- 0
    spread[optimized] <- 1
  }
  list(mean = center, sd = spread)
}

# Validate the whole split plan before starting any fold fits. IDs, rather than
# changing row positions, bind raw data, spatial buffers and assessment output.
.msr_cv_folds <- function(folds, ids) {
  valid_ids <- function(x) is.character(x) && length(x) > 0L &&
    !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x) && all(x %in% ids)
  if (!is.list(folds) || !length(folds) || is.null(names(folds)) ||
      anyNA(names(folds)) || any(!nzchar(names(folds))) || anyDuplicated(names(folds)))
    stop("folds must be a nonempty, uniquely named list of explicit splits.", call. = FALSE)
  for (name in names(folds)) {
    split <- folds[[name]]
    if (!is.list(split) || !setequal(names(split), c("analysis", "assessment")) ||
        length(split) != 2L || !valid_ids(split$analysis) || !valid_ids(split$assessment))
      stop("Fold `", name, "` must contain analysis and assessment character row IDs, without missing, duplicate or unknown IDs.", call. = FALSE)
    if (length(intersect(split$analysis, split$assessment)))
      stop("Fold `", name, "` has overlapping analysis and assessment rows.", call. = FALSE)
  }
  folds
}

# Binning is itself learned preparation. Its physical domain is declared by
# max_D; representative distances use training buffers only and stay frozen.
.msr_cv_bins <- function(inputs, reference = NULL) {
  specs <- inputs$scale_vars[inputs$scale_vars$type == "kernel", , drop = FALSE]
  if (is.null(inputs$binned) || !nrow(specs)) return(NULL)
  n <- length(inputs$d_list)
  nb <- inputs$binned$nbins
  if (!is.numeric(nb) || length(nb) != 1L || !is.finite(nb) || nb < 2 || nb != floor(nb))
    stop("CV binning requires a valid prepared nbins value.", call. = FALSE)
  domain <- inputs$max_D / inputs$unit_conv
  # Extraction selects intersecting cells, whose centers may exceed max_D.
  # A fixed cell-diagonal envelope retains them without using assessment data.
  cell_size <- inputs$cell_size
  if (!is.null(cell_size)) {
    if (!is.numeric(cell_size) || length(cell_size) != 2L ||
        any(!is.finite(cell_size) | cell_size <= 0))
      stop("CV cell_size must declare two finite positive grid dimensions.", call. = FALSE)
    cell_size <- unname(as.numeric(cell_size))
    half_diagonal <- (max(cell_size) / 2) * sqrt(sum((cell_size / max(cell_size))^2))
    domain <- (inputs$max_D + half_diagonal) / inputs$unit_conv
    if (!is.finite(domain) || domain <= 0)
      stop("Invalid declared CV extraction geometry.", call. = FALSE)
  }
  if (is.null(reference)) {
    breaks <- seq(0, domain * (1 + 1e-9), length.out = nb + 1L)
  } else {
    breaks <- reference$breaks
    if (length(breaks) != nb + 1L || any(!is.finite(breaks)) || any(diff(breaks) <= 0) ||
        length(reference$bin_rep) != nb || any(!is.finite(reference$bin_rep)) ||
        !identical(reference$domain, domain) || !identical(reference$cell_size, cell_size))
      stop("Invalid frozen training bin reference.", call. = FALSE)
  }
  codes <- lapply(inputs$d_list, function(d) {
    if (any(!is.finite(d) | d < 0 | d > domain * (1 + 1e-9)))
      stop("Prepared cell distances must lie inside the declared extraction envelope; rebuild legacy inputs with kernel_prep() to record grid geometry.", call. = FALSE)
    .bincode(d, breaks, include.lowest = TRUE)
  })
  if (is.null(reference)) {
    total <- count <- numeric(nb)
    for (i in seq_len(n)) {
      total <- total + .msr_bin_accumulate(inputs$d_list[[i]], codes[[i]], nb)
      count <- count + tabulate(codes[[i]], nb)
    }
    rep <- total / count
    empty <- !is.finite(rep)
    rep[empty] <- ((breaks[-1L] + breaks[-length(breaks)]) / 2)[empty]
    reference <- list(breaks = breaks, bin_rep = rep, domain = domain)
    if (!is.null(cell_size)) reference$cell_size <- cell_size
  }
  vsum <- csum <- stats::setNames(vector("list", nrow(specs)), specs$covariate)
  for (j in seq_len(nrow(specs))) {
    V <- C <- matrix(0, n, nb)
    for (i in seq_len(n)) {
      value <- as.numeric(inputs$raw_cov[[i]][, specs$source[j]])
      if (length(value) != length(codes[[i]]))
        stop("Prepared distances and source-cell values must align.", call. = FALSE)
      keep <- is.finite(value)
      V[i, ] <- .msr_bin_accumulate(value[keep], codes[[i]][keep], nb)
      C[i, ] <- tabulate(codes[[i]][keep], nb)
    }
    vsum[[j]] <- V; csum[[j]] <- C
  }
  list(bin_rep = reference$bin_rep, vsum = vsum, csum = csum,
    covariates = specs$covariate, sources = stats::setNames(specs$source, specs$covariate),
    nbins = nb, n_points = n, point_ids = rownames(inputs$kernel_dat),
    cv_reference = reference)
}

.msr_cv_inputs <- function(inputs, ids, covariates, reference = NULL) {
  rows <- match(ids, rownames(inputs$kernel_dat))
  out <- inputs
  out$kernel_dat <- inputs$kernel_dat[rows, covariates, drop = FALSE]
  out$d_list <- inputs$d_list[rows]
  out$raw_cov <- inputs$raw_cov[rows]
  if (!is.null(inputs$raster_group)) out$raster_group <- inputs$raster_group[ids]
  names(out$d_list) <- names(out$raw_cov) <- ids
  specs <- inputs$scale_vars
  old_opt <- specs$covariate[specs$optimize]
  out$scale_vars <- specs[match(covariates, specs$covariate), , drop = FALSE]
  selected <- match(out$scale_vars$covariate[out$scale_vars$optimize], old_opt)
  out$sigma <- inputs$sigma[selected]
  out$shape <- if (is.null(inputs$shape)) NULL else inputs$shape[selected]
  out$n_covs <- length(selected)
  out$binned <- .msr_cv_bins(out, reference)
  out
}

.msr_cv_initial_values <- function(inputs, fit = NULL) {
  specs <- inputs$scale_vars
  sigma <- inputs$sigma
  shape <- inputs$shape
  if (!is.null(fit$diagnostics$radius_profile)) {
    # CV preparation can order predictors differently from the native formula.
    # Use the selected coordinates by name and preserve fitted hard cutoffs.
    optimized <- specs$covariate[specs$optimize]
    fitted_names <- rownames(fit$scale_est)
    idx <- match(optimized, fitted_names)
    if (anyNA(idx))
      stop("The assessment predictors do not match the fitted radius grid.", call. = FALSE)
    sigma <- fit$optim_results$par[idx] * (fit$kernel_inputs$unit_conv / inputs$unit_conv)
    if (inputs$kernel == "expow") shape <- fit$optim_results$par[length(fitted_names) + idx]
  }
  value <- matrix(vapply(seq_along(inputs$d_list), function(i)
    .msr_eval_scale_vars(inputs$d_list[[i]], inputs$raw_cov[[i]], specs,
      sigma, shape, inputs$kernel, inputs$unit_conv,
      inputs$resolution, inputs$n_cols), numeric(nrow(specs))),
    nrow = length(inputs$d_list), ncol = nrow(specs), byrow = TRUE)
  dimnames(value) <- list(rownames(inputs$kernel_dat), specs$covariate)
  if (any(!is.finite(value))) stop("Initial spatial covariates must be finite within each fold.", call. = FALSE)
  value
}

# Restrict outcomes to an explicit recorded observation definition. Predictor
# equations remain native GLM formulas, including learned poly/spline metadata.
.msr_cv_outcome <- function(formula, data, family) {
  if (.msr_cv_distribution_family(family)) return(.msr_cv_distribution_outcome(formula, data, family))
  expr <- formula[[2L]]
  while (is.call(expr) && identical(expr[[1L]], as.name("I"))) expr <- expr[[2L]]
  simple <- is.symbol(expr) || (is.call(expr) && identical(expr[[1L]], as.name("cbind")) &&
    length(expr) == 3L && all(vapply(as.list(expr)[-1L], is.symbol, logical(1))))
  if (!simple) stop("CV responses must be a data column or cbind(successes, failures); predefine other observation transformations explicitly.", call. = FALSE)
  y <- eval(expr, envir = data, enclos = environment(formula))
  if (is.matrix(y) && ncol(y) == 1L) y <- y[, 1L]
  trials <- rep(1, nrow(data))
  if (family == "binomial") {
    if (is.factor(y)) {
      if (nlevels(y) != 2L) stop("A binomial response factor must declare two levels.", call. = FALSE)
      y <- as.integer(y != levels(y)[1L])
    }
    if (is.matrix(y) && ncol(y) == 2L) {
      if ((!is.numeric(y) && !is.logical(y)) || any(!is.finite(y) | y < 0 | y != floor(y)))
        stop("Binomial successes and failures must be finite nonnegative integer counts.", call. = FALSE)
      trials <- rowSums(y)
      y <- y[, 1L]
    } else if ((!is.numeric(y) && !is.logical(y)) || !is.null(dim(y)) || anyNA(y) || any(!y %in% c(0, 1))) {
      stop("Use binary binomial outcomes or cbind(successes, failures), not weighted proportions.", call. = FALSE)
    }
  } else {
    if ((!is.numeric(y) && !is.logical(y)) || !is.null(dim(y)) || any(!is.finite(y)))
      stop("Gaussian and Poisson outcomes must be finite numeric vectors.", call. = FALSE)
    if (family == "poisson" && any(y < 0 | y != floor(y)))
      stop("Poisson outcomes must be nonnegative integer counts.", call. = FALSE)
  }
  if (length(y) != nrow(data)) stop("Response rows must match the explicit data rows.", call. = FALSE)
  list(observed = as.numeric(y), trials = as.numeric(trials))
}

.msr_cv_score <- function(predictions, family, metrics) {
  if (identical(family, "conditional_choice")) return(.msr_choice_score(predictions, metrics))
  if (.msr_cv_distribution_family(family) || "distribution_family" %in% names(predictions))
    return(.msr_cv_distribution_score(predictions, metrics))
  y <- predictions$observed; mu <- predictions$prediction
  n <- predictions$trials; w <- predictions$weight
  if (any(!is.finite(mu)) || (family == "binomial" && any(mu < 0 | mu > 1)) ||
      (family == "poisson" && any(mu < 0)))
    stop("Held-out predictions must be finite and inside the family's response domain.", call. = FALSE)
  eligible <- w > 0 & (family != "binomial" | n > 0)
  if (!any(eligible)) stop("The assessment set has no positive-weight observations or trials.", call. = FALSE)
  predictions$eligible <- eligible
  for (metric in metrics) {
    value <- switch(metric, rmse = (y - mu)^2, mae = abs(y - mu),
      neg_log_score = if (family == "poisson") -stats::dpois(y, mu, log = TRUE) else
        -stats::dbinom(y, n, mu, log = TRUE),
      deviance = stats::poisson()$dev.resids(y, mu, rep(1, length(y))),
      brier = (y * (1 - mu)^2 + (n - y) * mu^2) / pmax(n, 1))
    value[!eligible] <- NA_real_
    predictions[[metric]] <- value
  }
  predictions
}

.msr_cv_metric <- function(predictions, metric, family, complete, expected) {
  if (identical(family, "conditional_choice")) return(.msr_choice_metric(predictions, complete))
  keep <- predictions$eligible & !is.na(predictions[[metric]])
  w <- predictions$weight[keep]
  if (metric == "brier") w <- w * predictions$trials[keep]
  value <- if (complete && length(w) && sum(w) > 0)
    sum(w * predictions[[metric]][keep]) / sum(w) else NA_real_
  if (metric == "rmse") value <- sqrt(value)
  data.frame(metric = metric, estimate = value,
    unit = if (metric == "brier") "trial" else "observation",
    n_assessment = expected, n_scored = sum(keep), weight_sum = sum(w), complete = complete,
    stringsAsFactors = FALSE)
}

#' Cross-Validate Estimated Scales with Explicit Folds
#'
#' Estimate the scales of effect and the GLM coefficients again inside every
#' fold, using only that fold's training rows, then score the predictions for
#' the rows it holds out. A *fold* is one split of the data into *analysis*
#' rows, which fit the model, and *assessment* rows, which are held out and
#' predicted. Use this function to ask how well a scale-of-effect model
#' predicts observations it never saw.
#'
#' @param formula Two-sided model formula, written as for [stats::glm()]. The
#'   response must be a column of `data`, optionally wrapped in `I()`, or
#'   `cbind(successes, failures)` for grouped binomial data. Create any other
#'   transformed response as a column before the call. Predictors can include
#'   factors, interactions, offsets, and the data-dependent terms described in
#'   "Formula terms that learn from the data". Every variable must be a column
#'   of `data` or a covariate prepared in `kernel_inputs`, and at least one
#'   predictor must be a covariate whose scale is optimized.
#' @param data Data frame with one row per observation, holding the raw
#'   responses and every nonspatial predictor. Its row names are the row IDs.
#'   They must be unique and must be the same set of IDs as the row names of
#'   `kernel_inputs$kernel_dat`. The folds refer to rows by these IDs. Spatial
#'   covariate columns are rebuilt inside every fold, so any values `data`
#'   holds for them are ignored. Rows with missing values are never dropped
#'   silently, so a fold that contains one fails.
#' @param kernel_inputs A `multiScaleR_data` object from [kernel_prep()] or
#'   [kernel_prep_by_group()] for the same observations. It must keep the
#'   cell-level distances and values (`store_cell_data = TRUE`, the default),
#'   because every fold rebuilds its covariates from the raw cells. Lean
#'   objects are rejected. Choose `max_D`, the starting scales, and any metric
#'   definitions without looking at the assessment outcomes.
#' @param folds Named list of explicit splits, for example from
#'   [grouped_cv_folds()] or [spatial_cv_folds()]. Each split is a list of two
#'   character vectors of row IDs, `analysis` and `assessment`, which must not
#'   overlap. A row can be left out of a split entirely (for example, by a
#'   spatial buffer) and can be assessed in more than one split. Keep rows that
#'   share a site, animal, or other sampling unit on the same side of every
#'   split.
#' @param family Character response family. Choose `"gaussian"` (identity
#'   link, the default), `"poisson"` (log), `"binomial"` (logit), `"nbinom2"`
#'   (negative binomial, log), or `"Gamma"` (log). Lowercase `"gamma"` is
#'   accepted. For `"nbinom2"` the variance is mu + mu^2 / theta, where mu is
#'   the expected count and theta the size parameter, so a small theta means
#'   strong extra variation. For `"Gamma"` the variance is dispersion x mu^2.
#'   Other links and quasi-likelihood families are not supported.
#' @param backend Fitting route for the scales, passed to
#'   [multiScale_optim()]. The default `"profile_glm"` is the fast analytic
#'   route. `"refit"` refits the GLM at every trial scale and accepts formula
#'   terms the default rejects, such as a spline of an optimized covariate
#'   without fixed knots. `"joint_nb2"` estimates the scales, coefficients,
#'   and theta together, and applies only to `family = "nbinom2"` with `theta`
#'   estimated. If you leave `backend` unset, `"Gamma"` and an NB2 model with
#'   a fixed `theta` use `"refit"`, the only route they support. Naming another
#'   backend for them is an error. The usual starting-value, convergence, and
#'   kernel or metric checks of [multiScale_optim()] apply in every fold. See
#'   [multiScaleR-fitting-methods].
#' @param weights Optional case frequencies, one nonnegative whole number per
#'   row of `data` in the same order. If the vector has names, they must equal
#'   the row names. A weight of 3 counts the row three times. A weight of 0
#'   keeps the row in the records but out of the fit and the scores. The
#'   default `NULL` gives every row a weight of one. Gaussian models accept
#'   only weights of one. For binomial data, put the trial counts in a
#'   `cbind(successes, failures)` response, not in `weights`. Weighted
#'   proportions are rejected. Some scale-fitting backends add their own
#'   weight requirements.
#' @param optim_args Named list of extra [multiScale_optim()] arguments used
#'   in every fold, such as `par` or `profile_control`. Supply `radius_grid`
#'   with `backend = "refit"` to profile hard radii, together with any smooth
#'   scales, in each training split. Supply `covariate_transform = "raw"` with
#'   `backend = "refit"` to keep optimized covariates in their own units (see
#'   "Keep optimized covariates in their own units"). Arguments that supply
#'   data, the starting model, a custom refit function, parallel workers,
#'   `backend`, or `verbose` are reserved. Every fold runs on one worker. The
#'   default `list()` keeps the backend's own defaults. Choose starting values
#'   without using the assessment outcomes.
#' @param glm_control Named list of [stats::glm.control()] settings for the
#'   initial GLM fitted to the analysis rows of each fold. The default keeps
#'   the settings of `glm()`.
#' @param metrics Optional character vector of scores to record. The defaults
#'   depend on `family`. They are `c("rmse", "mae")` for Gaussian,
#'   `c("rmse", "mae", "neg_log_score", "deviance")` for Poisson,
#'   `c("neg_log_score", "brier")` for binomial, and
#'   `c("rmse", "mae", "neg_log_score")` for NB2 and Gamma. You can request a
#'   subset of the defaults, but not other scores. NB2 and Gamma leave out the
#'   deviance because its value depends on theta or the dispersion, which
#'   makes it unsuitable for comparing models. Lower is better for every
#'   score, and the [multiScaleR_cv-class] page defines each one.
#' @param keep_fits Logical. The default `FALSE` keeps the predictions, scores,
#'   fold estimates, and compact diagnostics. `TRUE` also keeps the full fitted
#'   model of every fold in `fits`, which can use much more memory. You need
#'   full fits only to inspect the fold models or to have [select_cv()] return
#'   the selected fits.
#' @param keep_selection Logical. `TRUE` keeps a compact record of each fold's
#'   training likelihood and parameter count in `selection_records`, which is
#'   what [select_cv()] needs when `keep_fits = FALSE`. The records come from
#'   the fits already run, so they add no fitting time. The default is
#'   `FALSE`.
#' @param theta Optional NB2 size parameter fixed in advance, one positive
#'   number. The default `NULL` estimates theta separately from the analysis
#'   rows of each fold. A fixed theta adds no parameter to K (the parameter
#'   count used by AIC) and requires the `"refit"` backend. Fix it only from
#'   information independent of the assessment outcomes. Other families
#'   reject it.
#' @param dispersion Optional Gamma dispersion fixed in advance, one positive
#'   number. The default `NULL` estimates it in each fold as the training
#'   deviance divided by the sum of the frequency weights, the value native
#'   Gamma AIC uses. A fixed dispersion adds no parameter to K and is built
#'   into the Gamma density before fitting. Fix it only from information
#'   independent of the assessment outcomes. Other families reject it.
#' @return `cross_validate()` returns a list of class `multiScaleR_cv`, with
#'   class `multiScaleR_distribution_cv` added for `"nbinom2"` and `"Gamma"`.
#'   Read `complete` first, then `summary` (pooled scores), `metrics` (scores
#'   per fold), `predictions` (one row per held-out observation),
#'   `fold_results`, and `training` (the scales and coefficients estimated in
#'   each fold). See [multiScaleR_cv-class] for every element and how to read
#'   each score. `print()` shows how many folds were scored and the pooled
#'   scores, and returns `x` invisibly.
#' @details
#' # How each fold is fitted
#'
#' Suppose you estimated sigma, the kernel scale, once from all of the data
#' and then cross-validated only the regression. Every held-out row would
#' already have helped choose the scale that defines its own predictors, and
#' the scores would flatter the model. `cross_validate()` avoids this by
#' repeating the whole fit inside each fold, using the analysis rows alone.
#'
#' 1. Rebuild each spatial covariate from the raw cells at the starting scales
#'    in `kernel_inputs`, then center and scale it with the analysis-row mean
#'    and standard deviation. Each covariate must have a positive, finite
#'    standard deviation in the analysis rows.
#' 2. Fit the initial GLM to the analysis rows.
#' 3. Estimate the scales and coefficients with [multiScale_optim()]. No fit
#'    to the full data is used as a starting point.
#' 4. Rebuild the assessment covariates with the fitted scales, the training
#'    standardization, and the training distance bins, then predict the
#'    expected response and score it.
#'
#' The fold fails if the initial GLM, the final GLM, or the scale optimizer
#' does not converge. Folds run one after another, and a failed fold is
#' recorded and never retried.
#'
#' # Build the folds
#'
#' Row names tie the pieces together. Before any fitting, the function checks
#' that every split names known, unique IDs and that its analysis and
#' assessment rows do not overlap. It cannot check ecological independence
#' from the IDs you supply. Keep observations that share a site, animal, year,
#' or other sampling unit on the same side of each split, and choose the
#' design before you look at any scores.
#'
#' Rows left out of a split, such as those removed by a buffer, are neither
#' fitted nor scored in that split. A row assessed in several splits is
#' scored each time, every occurrence enters the pooled scores, and
#' `assessment_counts` records how often each row was assessed. A factor level
#' present in the assessment rows but absent from the analysis rows makes that
#' fold fail at prediction. No zero effect is substituted.
#' `vignette("v08-cross-validation", package = "multiScaleR")` compares
#' random, grouped, spatial, and buffered fold designs on the same model.
#'
#' # Choose a family and scores
#'
#' | `family` | Link | Default `metrics` | Backend when unset |
#' |:---|:---|:---|:---|
#' | `"gaussian"` | identity | `rmse`, `mae` | `"profile_glm"` |
#' | `"poisson"` | log | `rmse`, `mae`, `neg_log_score`, `deviance` | `"profile_glm"` |
#' | `"binomial"` | logit | `neg_log_score`, `brier` | `"profile_glm"` |
#' | `"nbinom2"` | log | `rmse`, `mae`, `neg_log_score` | `"profile_glm"`, or `"refit"` with a fixed `theta` |
#' | `"Gamma"` | log | `rmse`, `mae`, `neg_log_score` | `"refit"` |
#'
#' Lower is better for every score. RMSE and MAE judge only the predicted
#' mean, in the units of the response. The negative log score judges the
#' whole predicted distribution, so it also penalizes a model that is too
#' confident. Binomial responses can be 0/1 numbers, logical values, a
#' two-level factor, or `cbind(successes, failures)`. [cross_validate_glm()]
#' explains how binomial log and Brier scores are weighted, and
#' the [multiScaleR_cv-class] page defines every score.
#'
#' ## Negative binomial and Gamma responses
#'
#' Both families carry a nuisance parameter, a parameter of the response
#' distribution that is not part of the mean model. It is theta for NB2 and
#' the dispersion for Gamma. Unless you fix it with `theta` or `dispersion`,
#' each fold estimates it from the analysis rows and freezes it to score the
#' assessment rows. The negative log score therefore uses a complete,
#' normalized density that the held-out outcomes did not help shape. The
#' frozen value of each fold is stored in `distribution` and repeated in the
#' `theta` and `dispersion` columns of `predictions`.
#'
#' * An estimated nuisance parameter adds one to the native K, and a fixed one
#'   adds zero. The native `nobs` counts rows with positive weight, not the
#'   sum of the weights.
#' * Estimated Gamma dispersion is the training deviance divided by the sum of
#'   the frequency weights, the value native Gamma AIC uses. It differs from
#'   the Pearson dispersion that `summary()` uses for coefficient standard
#'   errors, and it does not replace that standard-error convention. A fixed
#'   Gamma dispersion is built into the Gamma density before fitting, so the
#'   likelihood, AIC, and K agree. Setting only the `dispersion` field of a
#'   native Gamma family would leave its AIC estimating the dispersion from
#'   the deviance.
#' * NB2 responses must be nonnegative whole numbers and Gamma responses
#'   strictly positive, including rows with zero weight.
#' * A theta warning from the native NB2 fit, such as a reached iteration
#'   limit, makes the fold fail. The warning is kept in `warnings`.
#' * `backend = "joint_nb2"` starts the coefficients from a GLM with
#'   theta = 1 and then estimates theta jointly with the scales. That
#'   starting value is neither a fixed theta nor a prior. A fold reports an
#'   infinite theta, the exact Poisson limit, only when the joint fit verifies
#'   and flags that limit. A large finite theta never stands in for it.
#'
#' A Gamma model for a strictly positive response uses the same call, and the
#' backend switches to `"refit"` on its own, for example
#' `cross_validate(mass ~ habitat, dat, inputs, folds, family = "Gamma")`.
#' Compare NB2 and Poisson candidates with the negative log score on the same
#' counts. [cv_tab()] lists the comparisons each family allows.
#'
#' # Formula terms that learn from the data
#'
#' Some formula terms learn constants from the data. `scale()` learns a mean
#' and standard deviation, `poly()` learns polynomial coefficients, and
#' `splines::ns()` and `splines::bs()` learn knots. In each fold these
#' constants come from the analysis rows and are reused unchanged for the
#' assessment rows. This holds with or without a `base::`, `stats::`, or
#' `splines::` prefix, and inside arithmetic, `I()`, pointwise functions such
#' as `log()`, and offsets. For spatial covariates the constants come from the
#' final fitted training values, not the starting-scale values. Arguments you
#' fix in the formula, such as `knots`, stay fixed. Coefficients, standard
#' errors, likelihoods, and K are those of the native fit, and fits kept with
#' `keep_fits = TRUE` carry the same training constants.
#'
#' Only these four functions are recognized, and only with numeric inputs and
#' supported arguments. An unsupported argument makes the fold fail. Any other
#' operation that learns from the data, such as `mean(x)` inside a formula or
#' your own transformation function, is not detected. Compute it within each
#' analysis split and apply it unchanged to the assessment rows, or the scores
#' will leak assessment information. For the same reason, do not supply
#' predictors standardized with all rows. Every variable must still be a
#' column of `data` or a prepared covariate. Variables are never taken from
#' the formula's environment.
#'
#' With the default `"profile_glm"` backend, a spline of an optimized
#' covariate needs fixed `knots` and `Boundary.knots` and no `df`. Otherwise
#' the fold fails with a message that names the fix. `backend = "refit"`
#' accepts adaptive splines.
#'
#' # Keep optimized covariates in their own units
#'
#' By default each fold centers and scales every optimized covariate with its
#' analysis-row mean and standard deviation, so a coefficient is the effect of
#' a one standard deviation change. With `backend = "refit"`, the setting
#' `optim_args = list(covariate_transform = "raw")` skips that final step for
#' the optimized covariates. They keep the units of their kernel or metric
#' calculation, such as a kernel-weighted proportion of habitat. Use this when
#' a fixed spline knot or another formula constant is meant in those units.
#'
#' The initial fit and every scale refit within a fold use the same units,
#' and prediction applies a center of zero and a scale of one. Covariates
#' whose scales are not optimized are still standardized with analysis-row
#' values, and normalization applied upstream, such as the reference
#' normalization of surface metrics, is unchanged. Assessment rows never
#' choose the transformation, the knots, or the scale. Decide on raw units,
#' and on any knots in those units, before you look at the results. Because a
#' change of units can alter a constrained or nonlinear model, compare the
#' two versions by their held-out scores rather than assume they are
#' equivalent.
#'
#' # Use separate maps for groups of surveys
#'
#' Inputs from [kernel_prep_by_group()], which take each row's covariates
#' from its own map (for example, the map of its survey year), are supported.
#' Before fitting, the function checks that the named `raster_group` labels
#' align with the prepared row IDs. It records the map labels of each fold's
#' analysis and assessment rows in `raster_groups`. Standardization and
#' distance bins still come from the pooled analysis rows. The map labels are
#' fixed inputs. They add no random effects, no extra parameters, and no
#' scales that vary among maps. Group effects you write in the formula, such
#' as a year factor, keep their native fit, parameter count, and new-level
#' check. See `vignette("v06-year-specific-rasters", package = "multiScaleR")`.
#'
#' # Distance bins and the extraction envelope
#'
#' When `kernel_inputs` holds distance-binned summaries ([kernel_prep()] with
#' `bin = TRUE`, the default), each fold rebuilds them. The bin edges span
#' `max_D` plus half a cell diagonal, so cells that intersect the buffer but
#' whose centers lie just beyond `max_D` keep their values. These edges depend
#' on the declared geometry only, not on any data. The representative
#' distance of each bin is the mean cell distance in that bin over the
#' analysis rows, and the assessment rows reuse it. Inputs saved by older
#' versions of [kernel_prep()] lack the cell geometry and keep the stricter
#' rule that every cell center lies within `max_D`. If such inputs are
#' rejected, prepare them again with the current [kernel_prep()].
#'
#' # Check every fold before comparing models
#'
#' Read `complete` first. `TRUE` means every planned fold was fitted and
#' scored. When a fold fails, `complete` is `FALSE` and the pooled scores in
#' `summary` are `NA`. The scores of the folds that finished, the failed
#' assessment rows, the diagnostics, and the warnings all stay available, and
#' `errors` names the stage and message of each failure. [cv_tab()] and
#' [select_cv()] refuse an incomplete result rather than rank a model on its
#' easier folds. A complete result says only that the calculations finished,
#' so inspect `diagnostics` and `warnings` before you rank models.
#'
#' # What the scores do and do not include
#'
#' The predictions and scores treat each fold's fitted parameters as known.
#' They do not average over uncertainty in the coefficients or scales, give
#' prediction intervals, or establish calibration. Compare candidates only on
#' identical folds, assessment rows, responses, and weights, which [cv_tab()]
#' checks for you. The scores complement [aic_tab()] but are not AIC values
#' or Akaike weights. If you use the scores to choose a model, kernel, or
#' metric, the winning score is optimistic. Assess the chosen procedure with
#' an outer validation layer or an independent test set. [select_cv()] scores
#' a predeclared training-AIC rule without this bias, because the assessment
#' outcomes never choose the winner within a fold.
#'
#' `cross_validate()` fits GLMs. For other model structures use
#' [cross_validate_mixed()] (a random intercept for repeated surveys),
#' [cross_validate_clogit()] (step selection), [cross_validate_field()] (a
#' spatial random field), or [cross_validate_bayes()] (Bayesian models). GAMs
#' and `unmarked` models have no cross-validation function.
#' [cross_validate_glm()] scores reference models without scales on the same
#' folds.
#'
#' @examples
#' data("landscape_counts")   # counts and one site-level covariate
#' data("surv_pts")           # 100 survey locations
#' pts <- vect(surv_pts)
#'
#' # Folds are lists of row IDs, so build and check them before any fitting.
#' # The simulated coordinates are in meters, and crs = 32617 declares that.
#' folds <- spatial_cv_folds(crds(pts), block_size = 1000, k = 3, crs = 32617)
#' folds   # rows fitted (analysis) and held out (assessment) in each fold
#'
#' \donttest{
#' land <- rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))
#' inputs <- kernel_prep(pts = pts, raster_stack = land[[c("land1", "land2")]],
#'                       max_D = 1000, kernel = "gaussian", verbose = FALSE)
#' # Row names of dat match those of inputs$kernel_dat ("1" to "100")
#' dat <- data.frame(landscape_counts, inputs$kernel_dat)
#'
#' # Scales and coefficients are estimated from the analysis rows of each fold
#' cv <- cross_validate(counts ~ site + land1 + land2, data = dat,
#'                      kernel_inputs = inputs, folds = folds,
#'                      family = "poisson", keep_selection = TRUE)
#' cv            # pooled held-out scores, lower is better
#' cv$complete   # TRUE: every fold was fitted and scored
#'
#' # Sigma (map units) estimated in each fold. Similar values across folds
#' # suggest that no single block of surveys drives the estimated scales.
#' round(sapply(cv$training, function(fold) fold$scale_est$Mean))
#'
#' # Compare with a reference model that has no landscape covariates
#' site_cv <- cross_validate_glm(counts ~ site, dat, folds, family = "poisson",
#'                               keep_selection = TRUE)
#' cv_tab(list(landscape = cv, site_only = site_cv), metric = "rmse")
#' }
#' @seealso [multiScaleR_cv-class] for the returned object,
#'   [cross_validate_glm()] for reference models, [cv_tab()] and [select_cv()]
#'   to compare candidates, [grouped_cv_folds()] and [spatial_cv_folds()] to
#'   build folds, and `vignette("v08-cross-validation", package = "multiScaleR")`
#'   for the full workflow.
#' @export
cross_validate <- function(formula, data, kernel_inputs, folds,
    family = c("gaussian", "poisson", "binomial", "nbinom2", "Gamma"),
    backend = c("profile_glm", "refit", "joint_nb2"), weights = NULL,
    optim_args = list(), glm_control = stats::glm.control(), metrics = NULL,
    keep_fits = FALSE, keep_selection = FALSE, theta = NULL, dispersion = NULL) {
  if (identical(family, "gamma")) family <- "Gamma"
  family <- match.arg(family)
  if (.msr_cv_distribution_family(family)) {
    selected_backend <- if (missing(backend) && (family == "Gamma" || !is.null(theta))) "refit" else match.arg(backend)
    out <- .msr_cv_distribution(formula, data, folds, family, weights, glm_control, metrics,
      keep_fits, keep_selection, theta, dispersion, kernel_inputs, selected_backend, optim_args)
    out$call <- match.call(); return(out)
  }
  if (!is.null(theta) || !is.null(dispersion)) stop("theta/dispersion controls apply only to NB2/Gamma CV.", call. = FALSE)
  backend <- match.arg(backend)
  if (backend == "joint_nb2") stop("joint_nb2 CV requires family = 'nbinom2'.", call. = FALSE)
  validate_scalar_logical(keep_fits, "keep_fits")
  validate_scalar_logical(keep_selection, "keep_selection")
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("formula must be a two-sided GLM formula.", call. = FALSE)
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)))
    stop("data must be a nonempty data frame with unique columns.", call. = FALSE)
  ids <- rownames(data)
  if (!.msr_valid_observation_ids(ids, nrow(data)))
    stop("data requires unique, nonmissing row IDs.", call. = FALSE)
  inputs <- kernel_inputs
  if (!inherits(inputs, "multiScaleR_data") || is.null(inputs$scale_vars) ||
      !is.list(inputs$d_list) || !is.list(inputs$raw_cov) ||
      length(inputs$d_list) != nrow(inputs$kernel_dat) ||
      length(inputs$raw_cov) != nrow(inputs$kernel_dat) ||
      !.msr_valid_observation_ids(rownames(inputs$kernel_dat), nrow(inputs$kernel_dat)) ||
      !setequal(ids, rownames(inputs$kernel_dat)))
    stop("CV requires matching kernel_prep() row IDs and stored cell-level inputs; use store_cell_data = TRUE.", call. = FALSE)
  if (!is.null(inputs$raster_group)) {
    groups <- inputs$raster_group
    if (!is.character(groups) || length(groups) != nrow(inputs$kernel_dat) ||
        !identical(names(groups), rownames(inputs$kernel_dat)))
      stop("raster_group must align exactly with prepared row IDs.", call. = FALSE)
    if (anyNA(groups) || any(!nzchar(groups)))
      stop("raster_group must contain nonmissing, nonempty fixed map labels.", call. = FALSE)
  }
  variables <- all.vars(formula)
  spatial <- inputs$scale_vars$covariate
  covariates <- if ("." %in% variables) spatial else intersect(spatial, variables)
  if (!any(inputs$scale_vars$optimize[match(covariates, spatial)]))
    stop("The formula must use an optimized spatial predictor.", call. = FALSE)
  if (length(setdiff(variables, c(names(data), spatial, "."))))
    stop("Supply every formula variable explicitly in data or spatial inputs.", call. = FALSE)
  if (length(intersect(all.vars(formula[[2L]]), spatial)))
    stop("A spatial predictor cannot be the CV response.", call. = FALSE)
  folds <- .msr_cv_folds(folds, ids)
  defaults <- list(gaussian = c("rmse", "mae"),
    poisson = c("rmse", "mae", "neg_log_score", "deviance"),
    binomial = c("neg_log_score", "brier"))
  if (is.null(metrics)) metrics <- defaults[[family]]
  if (!is.character(metrics) || !length(metrics) || anyNA(metrics) ||
      anyDuplicated(metrics) || any(!metrics %in% defaults[[family]]))
    stop("metrics must uniquely name scores supported for this family.", call. = FALSE)
  if (is.null(weights)) weights <- rep(1, nrow(data))
  if (!is.numeric(weights) || !is.null(dim(weights)) || length(weights) != nrow(data) ||
      any(!is.finite(weights) | weights < 0 | weights != floor(weights)) ||
      (!is.null(names(weights)) && !identical(names(weights), ids)) ||
      (family == "gaussian" && any(weights != 1)))
    stop("weights must align with data and be nonnegative integer frequencies; Gaussian weights must equal one.", call. = FALSE)
  allowed <- setdiff(names(formals(multiScale_optim)), c("fitted_mod", "kernel_inputs",
    "join_by", "n_cores", "PSOCK", "refit_fn", "backend", "verbose"))
  if (!is.list(optim_args) || (length(optim_args) &&
      (is.null(names(optim_args)) || anyNA(names(optim_args)) || anyDuplicated(names(optim_args)) ||
       any(!names(optim_args) %in% allowed))) ||
      (!is.null(optim_args$profile_control) && !is.list(optim_args$profile_control)) ||
      (is.list(optim_args$profile_control) && !is.null(optim_args$profile_control$worker_pool)))
    stop("optim_args must name supported serial scale controls; object/data/worker overrides are not allowed.", call. = FALSE)
  if (!is.list(glm_control) || is.null(names(glm_control)) || anyDuplicated(names(glm_control)) ||
      any(!names(glm_control) %in% names(formals(stats::glm.control))))
    stop("glm_control must be a named glm.control() list.", call. = FALSE)
  covariate_transform <- .msr_cv_covariate_transform(optim_args, backend)
  # Validate controls before any fit, without changing native defaults.
  glm_control <- do.call(stats::glm.control, glm_control)
  family_object <- switch(family, gaussian = stats::gaussian(), poisson = stats::poisson(), binomial = stats::binomial())
  output <- list(predictions = list(), metrics = list(), fold_results = list(),
    training = list(), diagnostics = list(), warnings = list(), errors = list(),
    prediction_contracts = list(), fits = if (keep_fits) list() else NULL,
    selection_records = if (keep_selection) list() else NULL)
  output$raster_groups <- if (is.null(inputs$raster_group)) NULL else list()
  for (name in names(folds)) {
    split <- folds[[name]]
    if (!is.null(inputs$raster_group))
      output$raster_groups[[name]] <- lapply(split, function(rows) inputs$raster_group[rows])
    ai <- match(split$analysis, ids); ti <- match(split$assessment, ids)
    pred <- data.frame(fold = name, row_id = split$assessment, observed = NA_real_, trials = NA_real_,
      prediction = NA_real_, weight = unname(weights[ti]), eligible = FALSE, stringsAsFactors = FALSE)
    for (metric in metrics) pred[[metric]] <- NA_real_
    fit <- initial <- NULL; messages <- character(); stage <- "preparation"
    error <- tryCatch(withCallingHandlers({
      outcome <- .msr_cv_outcome(formula, data[ti, , drop = FALSE], family)
      pred$observed <- outcome$observed; pred$trials <- outcome$trials
      .msr_cv_outcome(formula, data[ai, , drop = FALSE], family)
      train <- .msr_cv_inputs(inputs, split$analysis, covariates)
      raw <- .msr_cv_initial_values(train)
      center <- colMeans(raw); spread <- apply(raw, 2L, stats::sd)
      if (any(!is.finite(spread) | spread <= 0))
        stop("Training spatial predictors must have positive finite standard deviations.", call. = FALSE)
      normalization <- .msr_cv_raw_identity(center, spread, covariate_transform, train$scale_vars)
      center <- normalization$mean; spread <- normalization$sd
      train$kernel_dat <- as.data.frame(sweep(sweep(raw, 2L, center), 2L, spread, "/"), optional = TRUE)
      train$scl_params <- list(mean = center, sd = spread)
      training_data <- data[ai, , drop = FALSE]
      for (nm in covariates) training_data[[nm]] <- train$kernel_dat[[nm]]
      stage <- "initial_glm"
      initial <- do.call(stats::glm, list(formula = formula, data = training_data,
        family = family_object, weights = unname(weights[ai]), control = glm_control,
        na.action = stats::na.fail, model = TRUE))
      output$training[[name]] <- list(initial_normalization = train$scl_params,
        initial_coefficients = stats::coef(initial),
        bin_reference = if (is.null(train$binned)) NULL else train$binned$cv_reference)
      if (!isTRUE(initial$converged)) stop("Initial training GLM did not converge.", call. = FALSE)
      if (is.factor(stats::model.response(initial$model)) &&
          length(unique(stats::model.response(initial$model))) != 2L)
        stop("A training binomial factor must contain both declared response levels.", call. = FALSE)
      stage <- "scale_fit"
      fit <- do.call(multiScale_optim, c(list(fitted_mod = initial, kernel_inputs = train,
        n_cores = 1L, backend = backend, verbose = FALSE), optim_args))
      output$training[[name]] <- list(initial_normalization = train$scl_params,
        initial_coefficients = stats::coef(initial),
        fitted_normalization = fit$scl_params, coefficients = stats::coef(fit$opt_mod),
        scale_est = fit$scale_est, shape_est = fit$shape_est,
        bin_reference = if (is.null(train$binned)) NULL else train$binned$cv_reference,
        optimizer = fit$optim_results[c("par", "value", "convergence", "message")],
        start_results = fit$profile_fit$start_results)
      if (!isTRUE(fit$opt_mod$converged) || !identical(as.integer(fit$optim_results$convergence), 0L))
        stop("Training GLM or scale optimizer did not converge.", call. = FALSE)
      stage <- "prediction"
      fit <- .msr_cv_prediction_metadata(fit, training_data, train)
      assess <- .msr_cv_inputs(inputs, split$assessment, covariates,
        if (is.null(train$binned)) NULL else train$binned$cv_reference)
      raw_assess <- .msr_cv_initial_values(assess, fit)
      assess$kernel_dat <- as.data.frame(sweep(sweep(raw_assess, 2L, center), 2L, spread, "/"), optional = TRUE)
      assess$scl_params <- train$scl_params
      assessment_data <- data[ti, , drop = FALSE]
      for (nm in covariates) assessment_data[[nm]] <- assess$kernel_dat[[nm]]
      controls <- list(object = fit, newdata = assessment_data, kernel_inputs = assess,
        type = "response", se.fit = FALSE, uncertainty = "conditional", interval = "none")
      if (backend == "refit") controls$include_fixed <- TRUE
      mu <- do.call(stats::predict, controls)
      if (!is.numeric(mu) || length(mu) != nrow(pred)) stop("Prediction rows do not match the assessment set.", call. = FALSE)
      output$prediction_contracts[[name]] <- attr(mu, "prediction_contract", exact = TRUE)
      pred$prediction <- as.numeric(mu)
      stage <- "scoring"
      pred <- .msr_cv_score(pred, family, metrics)
      NULL
    }, warning = function(w) {
      messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) {
      result <- list(stage = stage, message = conditionMessage(e))
      if (!is.null(e$radius_profile)) result$radius_profile <- e$radius_profile
      result
    })
    complete <- is.null(error)
    if (!complete) {
      pred$prediction <- NA_real_; pred$eligible <- FALSE
      for (metric in metrics) pred[[metric]] <- NA_real_
    }
    output$predictions[[name]] <- pred
    output$metrics[[name]] <- do.call(rbind, lapply(metrics, function(metric)
      cbind(fold = name, .msr_cv_metric(pred, metric, family, complete, nrow(pred)))))
    output$fold_results[[name]] <- data.frame(fold = name, n_analysis = length(ai),
      n_assessment = length(ti), status = if (complete) "scored" else "failed", stringsAsFactors = FALSE)
    output$warnings[name] <- list(messages); output$errors[name] <- list(error)
    output$diagnostics[name] <- list(list(
      initial_glm = if (is.null(initial)) NULL else initial[c("converged", "boundary", "iter")],
      fitted_glm = if (is.null(fit)) NULL else fit$opt_mod[c("converged", "boundary", "iter")],
      scale = if (is.null(fit)) NULL else fit$diagnostics))
    if (keep_selection) output$selection_records[name] <- list(.msr_cv_selection_record(fit,
      output$diagnostics[[name]], .msr_cv_selection_controls(list(backend = backend,
        optim_args = optim_args, glm_control = glm_control))))
    if (keep_fits) output$fits[name] <- list(fit)
  }
  for (item in c("predictions", "metrics", "fold_results")) {
    output[[item]] <- do.call(rbind, output[[item]])
    rownames(output[[item]]) <- NULL
  }
  output$complete <- all(output$fold_results$status == "scored")
  output$summary <- do.call(rbind, lapply(metrics, function(metric)
    .msr_cv_metric(output$predictions, metric, family, output$complete, nrow(output$predictions))))
  output$folds <- folds
  output$assessment_counts <- stats::setNames(tabulate(match(output$predictions$row_id, ids), length(ids)), ids)
  output$targets <- .msr_cv_targets(formula, data, family, weights, folds)
  output$settings <- list(family = family, link = family_object$link, backend = backend,
    metrics = metrics, optim_args = optim_args, glm_control = glm_control,
    normalization = "training_only", bin_reference = "training_only_declared_buffer_domain",
    prediction_target = "held_out_GLM_expected_response", parameter_uncertainty = "not_integrated",
    observation_intervals = "not_provided", calibration = "not_established", workers = 1L)
  output$call <- match.call()
  class(output) <- "multiScaleR_cv"
  output
}

#' @rdname cross_validate
#' @param x A `multiScaleR_cv` result, for `print()`.
#' @param ... Not used. Reserved for future print controls.
#' @export
print.multiScaleR_cv <- function(x, ...) {
  cat("Cross-validation:", sum(x$fold_results$status == "scored"), "of",
    nrow(x$fold_results), "folds scored;", sum(x$assessment_counts > 0), "of",
    length(x$assessment_counts), "input rows assessed.\n")
  if (!x$complete) cat("Incomplete: inspect errors and diagnostics before comparing models.\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}
