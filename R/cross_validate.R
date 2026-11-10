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

.msr_cv_initial_values <- function(inputs) {
  specs <- inputs$scale_vars
  value <- matrix(vapply(seq_along(inputs$d_list), function(i)
    .msr_eval_scale_vars(inputs$d_list[[i]], inputs$raw_cov[[i]], specs,
      inputs$sigma, inputs$shape, inputs$kernel, inputs$unit_conv,
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
#' Refit ecological scales and GLM coefficients using only each training set,
#' then score predictions for the explicitly held-out rows. A fold is one
#' training and assessment split. This function supports ordinary Gaussian,
#' Poisson, binomial, NB2 and Gamma GLMs and runs folds sequentially.
#'
#' @param formula Two-sided GLM formula. The response is a data column or
#'   \code{cbind(successes, failures)}. An \code{I()} wrapper is also accepted.
#'   Predictor expressions, interactions, factors, offsets and supported learned
#'   transformations use the existing backend's native formula machinery.
#'   Supply every variable in \code{data} or the prepared spatial inputs.
#' @param data Data frame containing raw responses and nonspatial formula
#'   variables. Unique row names must match the prepared spatial row IDs.
#'   Prepared spatial columns are rebuilt inside each fold, replacing any
#'   supplied values. No fitting rows are silently omitted for missing values.
#' @param kernel_inputs A \code{multiScaleR_data} object from
#'   \code{\link{kernel_prep}} with cell-level distances and source values
#'   retained. Use \code{store_cell_data = TRUE}. Its buffer, initial scales,
#'   metric definitions and other controls must be chosen independently of the
#'   assessment outcomes. Modern \code{\link{kernel_prep_by_group}} inputs
#'   are supported when their named \code{raster_group} map aligns with prepared
#'   row IDs. Lean-only preparation is not supported.
#' @param folds Uniquely named list of explicit splits. Each split contains
#'   nonempty character vectors \code{analysis} and \code{assessment} of row
#'   IDs. They must be disjoint within a split. Splits can omit rows for spatial
#'   buffers and can repeat assessment rows across splits. Keep dependent sites,
#'   animals or other sampling groups together when defining the splits.
#' @param family Character response family: \code{"gaussian"} uses the identity
#'   link, \code{"poisson"}, \code{"nbinom2"} and \code{"Gamma"} use log links,
#'   and \code{"binomial"} uses the logit link. Lowercase \code{"gamma"} is an
#'   alias. NB2 means variance equals mean plus mean squared divided by theta;
#'   Gamma variance equals dispersion times mean squared. The default is
#'   \code{"gaussian"}. Other links and quasi-likelihood targets are unqualified.
#' @param backend Scale-fitting backend passed to \code{\link{multiScale_optim}}.
#'   \code{"profile_glm"} remains the default, except that Gamma and fixed-theta
#'   NB2 default to their required native \code{"refit"} route. Estimated NB2
#'   also supports \code{"joint_nb2"}. Explicit incompatible backends are
#'   rejected. Existing starts, convergence, cold/limiting, curvature and
#'   kernel/metric guards remain in force.
#' @param weights Optional nonnegative integer frequency weights, one per data
#'   row in its supplied order. Names, when present, must match that order.
#'   Gaussian weights must all equal one. Binomial trials belong in the
#'   two-column response, separately from frequency weights; weighted proportions
#'   are not accepted. A backend may impose stricter native weight requirements.
#' @param optim_args Named list of additional \code{multiScale_optim()} controls,
#'   such as \code{par} or \code{profile_control}. With spatial native refits,
#'   \code{covariate_transform = "raw"} keeps optimized predictors in their
#'   kernel/metric units for both initial training fits and later prediction.
#'   Fixed-metric transforms still use training rows only. Defaults preserve that backend's
#'   controls. Starts must be selected without assessment outcomes. Object/data
#'   arguments, custom refit functions, worker pools and parallel controls cannot
#'   be replaced. Each fold uses one fit worker.
#' @param glm_control Named list from \code{stats::glm.control()}, used for each
#'   initial training-only GLM. The default keeps native GLM controls.
#' @param metrics Optional character vector of score names. Defaults are
#'   \code{c("rmse", "mae")} for Gaussian, \code{c("rmse", "mae",
#'   "neg_log_score", "deviance")} for Poisson, and
#'   \code{c("neg_log_score", "brier")} for binomial. NB2/Gamma default to
#'   \code{c("rmse", "mae", "neg_log_score")}. Their nuisance-dependent
#'   deviances are not comparison metrics. Lower values are better.
#' @param keep_fits Logical. Default \code{FALSE} keeps compact predictions,
#'   scales and diagnostics. Set \code{TRUE} to retain full fold models in memory
#'   for inspection; this can use substantially more memory.
#' @param keep_selection Logical. Default \code{FALSE} does not retain native
#'   training-selection records. Set \code{TRUE} to use \code{\link{select_cv}}
#'   without keeping full fits. Records are captured from each actual training
#'   fit and add no model fitting or diagnostic eligibility.
#' @param theta Optional finite positive NB2 size fixed before CV, or \code{NULL}
#'   (default) to estimate it separately in each training set. Fixed theta adds
#'   no nuisance parameter to K and requires \code{backend = "refit"} for scales.
#'   Do not choose theta from assessment outcomes. Other families reject it.
#' @param dispersion Optional finite positive Gamma dispersion fixed before CV,
#'   or \code{NULL} (default) for the native training Gamma AIC convention,
#'   deviance divided by the sum of frequency weights. Fixed dispersion adds
#'   no nuisance parameter to K. It is supplied to a scoped known-dispersion
#'   family before fitting, with its normalized density and native coefficient
#'   SE convention. Other families reject it.
#' @return An object of class \code{multiScaleR_cv}, additionally inheriting
#'   \code{multiScaleR_distribution_cv} for NB2/Gamma, a list containing:
#'   \describe{
#'     \item{\code{summary}}{Pooled scores with columns \code{metric},
#'       \code{estimate}, scoring \code{unit}, planned \code{n_assessment},
#'       \code{n_scored}, \code{weight_sum} and \code{complete}. Estimates are
#'       \code{NA} if any planned fold fails.}
#'     \item{\code{metrics}}{The same score columns per \code{fold}. RMSE is the
#'       square root of pooled squared errors, not the average of fold RMSEs.}
#'     \item{\code{predictions}}{Every planned assessment occurrence, including
#'       failed folds: \code{fold}, \code{row_id}, \code{observed} response/count,
#'       \code{trials}, response-scale \code{prediction}, frequency \code{weight},
#'       \code{eligible} scoring flag and requested score columns. The \code{rmse}
#'       column stores squared errors. Failed predictions/scores are \code{NA}.
#'       NB2/Gamma add \code{distribution_family}, frozen training \code{theta}
#'       or \code{dispersion}, and the existing joint-NB2 \code{poisson_limit}
#'       flag. Inapplicable nuisance columns are \code{NA}.}
#'     \item{\code{distribution}}{NB2/Gamma training distribution metadata per
#'       fold: family, theta/dispersion, estimated/fixed convention, Poisson-limit
#'       flag, native rank/K, log likelihood, effective nobs, deviance and weight
#'       total. A failed training distribution may leave an element \code{NULL}.}
#'     \item{\code{fold_results}}{Per-fold \code{fold}, \code{n_analysis},
#'       \code{n_assessment} and \code{status} (\code{"scored"} or \code{"failed"}).}
#'     \item{\code{folds}, \code{assessment_counts}}{The explicit split plan and
#'       named counts of assessment occurrences for every input row, including
#'       unassessed rows. Repeated assessment occurrences each enter pooled scores.}
#'     \item{\code{training}}{Per-fold initial/fitted normalization,
#'       coefficients on the link scale, scales in map units, unitless shapes,
#'       frozen bin reference and optimizer/start
#'       records when available. Failed preparation can leave an element empty.}
#'     \item{\code{diagnostics}, \code{warnings}, \code{errors}}{Original scale-fit
#'       diagnostics, initial/final GLM convergence metadata and captured messages
#'       per fold. Failures are retained and never retried automatically.}
#'     \item{\code{prediction_contracts}}{Existing prediction contracts per fold,
#'       retaining the native conditioning and uncertainty distinctions.}
#'     \item{\code{raster_groups}}{For grouped raster inputs, named fixed-map
#'       labels for analysis/assessment IDs per fold, retained even on failure.
#'       Otherwise \code{NULL}. These are input maps, not fitted random effects.}
#'     \item{\code{targets}}{Compact used-row responses, trials and weights for
#'       guarded comparison with \code{\link{cv_tab}}. Never used in fitting.}
#'     \item{\code{selection_records}}{Native training likelihood, complete
#'       parameter counts, effective row IDs, outcomes/trials/frequencies, formula
#'       text, controls, diagnostics and metadata failures when
#'       \code{keep_selection = TRUE}; otherwise \code{NULL}. These compact records
#'       support \code{\link{select_cv}} without full fold models. Each record has
#'       \code{version}, \code{controls}, \code{diagnostics}, metadata \code{warnings},
#'       \code{data} and, on extraction failure, an \code{error} message.
#'       Successful \code{data} contains \code{logLik} on the natural-log scale;
#'       \code{native_K}, \code{ecological_K} and their total \code{K}; effective
#'       \code{nobs} and \code{observation_ids}; analysis-row \code{targets};
#'       \code{response_layout}, native likelihood \code{signature}, \code{formula}
#'       text, \code{family} and \code{link}. These are provenance from the fit;
#'       preserve them unchanged. Records cannot reconstruct native models.}
#'     \item{\code{fits}}{Full fold fits when \code{keep_fits = TRUE}; otherwise
#'       \code{NULL}. These fits use only their own training data.}
#'     \item{\code{complete}, \code{settings}, \code{call}}{Whether all planned
#'       folds were scored, family/backend/metric settings and the original call.
#'       Completion does not certify diagnostic eligibility or calibration.}
#'   }
#' @details
#' Raw spatial buffers are reused without changing them. Initial spatial
#' normalization and pooled representative bin distances are learned from the
#' training rows only. The bin domain is fixed from the declared buffer and,
#' for modern preparation, the raster cell-size envelope. Intersecting cells
#' can have centers beyond the nominal buffer; their values are retained.
#' The geometric domain is independent of assessment data. Legacy stored inputs
#' without cell geometry retain the original strict buffer-center gate; rebuild
#' them with the current preparation function if it rejects their distances.
#' Assessment values use the frozen training bins and normalization. Existing
#' coefficient, scale, shape and numerical-diagnostic calculations are unchanged.
#' Grouped raster maps select fixed source covariates by row. Their group map
#' is validated before fitting and subset with the buffers. Scales and
#' coefficients use the pooled training rows; assessment preprocessing remains
#' frozen. Map labels do not introduce annual/random effects, extra fitted
#' parameters or spatially varying scales. Explicit formula group effects keep
#' their native fitting, parameter counts and new-level prediction guards.
#'
#' Supported learned formula transformations keep analysis-row constants for
#' assessment prediction. These include bare or namespaced \code{scale()},
#' \code{poly()}, \code{splines::ns()} and \code{splines::bs()}, including calls
#' inside arithmetic, \code{I()}, pointwise functions and offsets. Spatial fits
#' use their final fitted training predictors, not initial-scale values. Original
#' formulas, coefficients, coefficient SEs, likelihoods, K and numerical guards
#' are unchanged. Retained fits carry the same repaired prediction terms.
#' Predeclared numeric transform arguments remain fixed, and native abbreviated
#' scale arguments keep their meanings. Formula variables must still be supplied
#' explicitly; this does not permit formula-environment predictor fallback.
#' Recovery evaluates only these known learners with numeric inputs and supported
#' pointwise arguments; arbitrary user functions are not re-executed to infer
#' preprocessing. Other learned operations, such as \code{mean(x)} inside a
#' formula or a custom transformation, need an explicit recipe trained within
#' each analysis split and applied unchanged to its assessment rows. Unsupported
#' arguments to a known learner produce a retained fold failure. Existing backend
#' restrictions on expressions involving optimized spatial predictors still apply.
#'
#' Negative log scores retain Poisson/binomial likelihood constants and average
#' over positive-weight observations using frequency weights. Binomial Brier
#' scores average squared prediction errors across trials, weighted by frequency
#' times trial count. Zero-trial rows contribute no score. Infinite log scores
#' correctly record an observed outcome assigned zero probability. Point scores
#' use fitted parameters; this function adds no parameter uncertainty, observation
#' intervals, calibration claim or new eligibility for guarded methods.
#'
#' Use identical folds and assessment targets when comparing candidate models.
#' If CV scores also choose a kernel, metric or model, assess the selected
#' procedure with outer validation or an independent test set. These scores
#' complement \code{\link{aic_tab}}; they are not AIC values or Akaike weights.
#'
#' Read \code{complete} first. Inspect failed folds and original diagnostics
#' before ranking models. Partial fold scores remain available for diagnosis,
#' while an incomplete pooled summary stays \code{NA}. This stage does not
#' score GAMs, mixed, unmarked or Bayesian fits. Use \code{\link{cross_validate_field}}
#' for the separately qualified Gaussian spatial-field workflow. See
#' \code{\link{grouped_cv_folds}} and \code{\link{spatial_cv_folds}} for
#' deterministic whole-group and spatial-grid fold construction.
#' NB2 and Gamma held-out negative log scores use complete normalized densities
#' with nuisance parameters learned only from analysis rows. Gamma uses its
#' native AIC deviance scale when dispersion is estimated; this differs from
#' the Pearson scale used for native coefficient SEs. Fixed Gamma constructs
#' its declared density before fitting, rather than setting only a native
#' dispersion field that would leave Gamma's AIC estimator inconsistent.
#' Estimated nuisance parameters each add one to native K; known ones add zero.
#' Native nobs counts positive-weight rows, not the sum of frequency weights.
#' NB2 counts must be nonnegative integers and Gamma outcomes strictly positive,
#' including zero-weight rows. Native theta/convergence warnings retain failed
#' folds. No score integrates parameter uncertainty or supplies observation
#' intervals. The exact Poisson limit is available only when the original
#' joint-NB2 fit verifies and flags it; no large finite theta substitutes for it.
#' Joint-NB2 initializes coefficients with a native fixed-size GLM (theta = 1),
#' then estimates dispersion through the original joint engine and boundary
#' checks. The initializer's fixed size is not a fixed fitted theta or a prior.
#'
#' @examples
#' # Explicit folds can be checked before any fitting is attempted.
#' ids <- paste0("site", 1:6)
#' folds <- list(first = list(analysis = ids[1:4], assessment = ids[5:6]),
#'               second = list(analysis = ids[3:6], assessment = ids[1:2]))
#' intersect(folds$first$analysis, folds$first$assessment) # no overlap
#' @seealso \code{\link{kernel_prep}}, \code{\link{multiScale_optim}},
#'   \code{\link{prepare_prediction}}, \code{\link{aic_tab}}
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
      raw_assess <- .msr_cv_initial_values(assess)
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
    }), error = function(e) list(stage = stage, message = conditionMessage(e)))
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
#' @param x A \code{multiScaleR_cv} result.
#' @param ... Reserved for future print controls.
#' @export
print.multiScaleR_cv <- function(x, ...) {
  cat("Cross-validation:", sum(x$fold_results$status == "scored"), "of",
    nrow(x$fold_results), "folds scored;", sum(x$assessment_counts > 0), "of",
    length(x$assessment_counts), "input rows assessed.\n")
  if (!x$complete) cat("Incomplete: inspect errors and diagnostics before comparing models.\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}
