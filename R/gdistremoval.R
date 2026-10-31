# GDR keeps inactive formulas in its call. Use its estimated components to
# distinguish singing availability from additional between-visit availability.
.msr_gdr_components <- function(model) {
  c("lambda", if (model@data@numPrimary > 1L) "phi",
    if (model@keyfun != "uniform") "dist", "rem")
}

.msr_gdr_prediction_type <- function(model, type = NULL) {
  if (is.null(type)) return("lambda")
  allowed <- .msr_gdr_components(model)
  if (!is.character(type) || length(type) != 1L || is.na(type) || !type %in% allowed)
    stop("GDR prediction type must be an active native component: ",
      paste(allowed, collapse = ", "),
      ". Use predict_gdistremoval() for combined detection and expected counts.", call. = FALSE)
  type
}

.msr_gdr_component_info <- function(model, type) {
  quantity <- switch(type,
    lambda = if (model@mixture == "ZIP") "poisson_component_intensity" else
      if (model@output == "density") "expected_density" else "expected_latent_abundance",
    phi = "between_visit_availability", rem = "per_time_unit_removal_probability",
    dist = switch(model@keyfun, halfnorm = "distance_standard_deviation",
      exp = "distance_decay_length", hazard = "hazard_distance_scale"))
  list(quantity = quantity, units = if (type %in% c("phi", "rem")) "probability" else
    if (type == "dist") model@data@unitsIn else
      if (model@output == "density") paste0("individuals_per_", model@unitsOut) else "individuals_per_survey_circle",
    mixture = model@mixture, native_output = model@output,
    unconditional_mean = type == "lambda" && model@mixture != "ZIP")
}

.msr_gdr_formulas <- function(model) {
  formulas <- model@formlist
  aliases <- c(lambda = "lambdaformula", phi = "phiformula",
               dist = "distanceformula", rem = "removalformula")
  out <- lapply(.msr_gdr_components(model), function(nm) {
    f <- formulas[[nm]]
    if (is.null(f)) f <- formulas[[aliases[[nm]]]]
    if (!inherits(f, "formula")) stop("Could not recover the active GDR ", nm,
      " formula. Refit the model with your installed unmarked version.", call. = FALSE)
    f
  })
  stats::setNames(out, .msr_gdr_components(model))
}

.msr_gdr_validate_spatial <- function(model, covs, frame = model@data) {
  formulas <- .msr_gdr_formulas(model)
  for (component in intersect(c("phi", "dist", "rem"), names(formulas))) {
    vars <- intersect(covs, all.vars(formulas[[component]]))
    higher <- names(frame@yearlySiteCovs)
    if (component == "rem") higher <- union(higher, names(frame@obsCovs))
    shadowed <- intersect(vars, higher)
    if (length(shadowed)) stop("Optimized GDR site covariates are shadowed by yearly or observation covariates: ",
      paste(shadowed, collapse = ", "), ". Use distinct covariate names.", call. = FALSE)
  }
  invisible(NULL)
}

# Match native covariate precedence, but use data-frame prediction separately
# for each component. This preserves missing prediction rows without calling
# getP(), whose historical GDR implementation reuses the first visit's sigma.
.msr_gdr_covariates <- function(frame) {
  m <- nrow(frame@yDistance)
  tt <- frame@numPrimary
  j <- ncol(frame@yRemoval) / tt
  bind <- function(a, b) {
    out <- cbind(a, b)
    out[, !duplicated(names(out)), drop = FALSE]
  }
  sc <- frame@siteCovs
  if (is.null(sc)) sc <- data.frame(.dummy = rep(0, m))
  annual <- frame@yearlySiteCovs
  if (is.null(annual)) annual <- data.frame(.dummy = rep(0, m * tt))
  annual <- bind(annual, sc[rep(seq_len(m), each = tt), , drop = FALSE])
  obs <- frame@obsCovs
  if (is.null(obs)) obs <- data.frame(.dummy = rep(0, m * tt * j))
  obs <- bind(obs, annual[rep(seq_len(m * tt), each = j), , drop = FALSE])
  list(lambda = sc, phi = annual, dist = annual, rem = obs)
}

.msr_gdr_distance <- function(key, parameter, breaks, n, hazard = NULL) {
  radius <- max(breaks)
  lower <- utils::head(breaks, -1L) / radius
  upper <- utils::tail(breaks, -1L) / radius
  if (key == "uniform") return(matrix(rep(upper^2 - lower^2, each = n), nrow = n))
  if (length(parameter) != n || any(parameter < 0, na.rm = TRUE))
    stop("GDR distance predictions must be nonnegative and follow site/period order.", call. = FALSE)
  if (key == "hazard" && (length(hazard) != 1L || !is.finite(hazard) || hazard <= 0))
    stop("The hazard-rate exponent must be finite and positive.", call. = FALSE)
  # Evaluate repeated parameter values once. Closed forms operate on all
  # unique scales together; only the hazard key needs numerical quadrature.
  # All integrals use r/radius and the native point-transect area measure.
  unique_par <- unique(parameter[!is.na(parameter)])
  probabilities <- matrix(0, length(unique_par), length(lower))
  ratio <- unique_par / radius
  flat <- is.infinite(ratio)
  if (key == "halfnorm") flat <- flat | (ratio > 0 & 0.5 / ratio^2 == 0)
  if (any(flat)) probabilities[flat, ] <- matrix(rep(upper^2 - lower^2,
    each = sum(flat)), nrow = sum(flat))
  active <- which(ratio > 0 & !flat)
  if (key == "halfnorm") {
    a <- 0.5 / ratio[active]^2
    # Preserve the zero-probability limit at scales below floating-point range.
    active <- active[is.finite(a)]
    a <- a[is.finite(a)]
    for (j in seq_along(lower)) probabilities[active, j] <-
      exp(-a * lower[j]^2) * (-expm1(-a * (upper[j]^2 - lower[j]^2))) / a
  } else if (key == "exp") {
    # Integral of 2*x*exp(-x/ratio) is a shape-two gamma probability times
    # 2*ratio^2. Choose CDF or survival differences to avoid cancellation.
    scales <- ratio[active]
    for (j in seq_along(lower)) {
      tail <- lower[j] >= scales
      for (use_tail in c(FALSE, TRUE)) {
        rows <- which(tail == use_tail)
        if (!length(rows)) next
        hi <- stats::pgamma(upper[j] / scales[rows], shape = 2,
          lower.tail = !use_tail, log.p = TRUE)
        lo <- stats::pgamma(lower[j] / scales[rows], shape = 2,
          lower.tail = !use_tail, log.p = TRUE)
        larger <- if (use_tail) lo else hi
        smaller <- if (use_tail) hi else lo
        log_mass <- rep(-Inf, length(rows))
        finite <- is.finite(larger)
        log_mass[finite] <- larger[finite] + log(-expm1(smaller[finite] - larger[finite]))
        probabilities[active[rows], j] <- exp(log(2) + 2 * log(scales[rows]) + log_mass)
      }
    }
  } else for (i in active) {
    scale <- ratio[i]
    integrand <-
      function(x) 2 * x * (-expm1(-exp(-hazard * log(x / scale))))
    probabilities[i, ] <- vapply(seq_along(lower), function(j) stats::integrate(integrand,
      lower[j], upper[j], rel.tol = 1e-8, subdivisions = 1000L)$value, numeric(1))
  }
  out <- matrix(NA_real_, nrow = n, ncol = length(lower))
  keep <- !is.na(parameter)
  if (any(keep)) out[keep, ] <- probabilities[match(parameter[keep], unique_par), , drop = FALSE]
  out
}

.msr_gdr_removal <- function(p, lengths) {
  out <- p
  survival <- rep(1, nrow(p))
  for (j in seq_along(lengths)) {
    interval <- -expm1(lengths[j] * log1p(-p[, j]))
    out[, j] <- survival * interval
    survival <- survival * (1 - interval)
  }
  out
}

.msr_gdr_probability <- function(x, label) {
  if (any(!is.na(x) & (!is.finite(x) | x < -1e-10 | x > 1 + 1e-10)))
    stop(label, " must be a probability. Check the fitted detection model.", call. = FALSE)
  # Only integration/roundoff excursions are corrected, never invalid fits.
  pmin(1, pmax(0, x))
}

.msr_gdr_prepare <- function(object, newdata, kernel_inputs, join_by, include_fixed) {
  validate_scalar_logical(include_fixed, "include_fixed")
  spatial <- inherits(object, "multiScaleR")
  model <- if (spatial) object$opt_mod else object
  if (!methods::is(model, "unmarkedFitGDR"))
    stop("object must be an unmarkedFitGDR or a multiScaleR fit containing one.", call. = FALSE)
  frame <- if (is.null(newdata)) model@data else newdata
  if (!methods::is(frame, "unmarkedFrameGDR"))
    stop("newdata must be an explicit unmarkedFrameGDR with the prediction covariates and survey protocol.", call. = FALSE)
  if (!is.null(kernel_inputs)) {
    if (!spatial || is.null(newdata))
      stop("kernel_inputs require a multiScaleR fit and an explicit prediction frame.", call. = FALSE)
    covs <- intersect(rownames(object$scale_est), .unmarked_model_predictors(model))
    .msr_gdr_validate_spatial(model, covs, frame)
    frame@siteCovs <- prepare_unmarked_data(object, frame@siteCovs,
      kernel_inputs, join_by, include_fixed)
  } else if (!is.null(join_by) || include_fixed)
    stop("join_by and include_fixed require kernel_inputs.", call. = FALSE)
  m <- nrow(frame@yDistance)
  tt <- frame@numPrimary
  breaks <- frame@dist.breaks
  lengths <- frame@period.lengths
  if (!m || length(tt) != 1L || !is.finite(tt) || tt < 1 || tt != round(tt) ||
      nrow(frame@yRemoval) != m || ncol(frame@yDistance) %% tt || ncol(frame@yRemoval) %% tt ||
      length(breaks) != ncol(frame@yDistance) / tt + 1L ||
      length(breaks) < 2L || any(!is.finite(breaks)) || breaks[1] != 0 || any(diff(breaks) <= 0) ||
      length(lengths) != ncol(frame@yRemoval) / tt || !length(lengths) ||
      any(!is.finite(lengths) | lengths <= 0 | lengths != round(lengths)) ||
      !frame@unitsIn %in% c("m", "km") || !model@unitsOut %in% c("ha", "kmsq"))
    stop("GDR prediction requires a valid site/period layout, increasing distance breaks starting at zero, and positive integer removal durations.", call. = FALSE)
  for (nm in c("siteCovs", "yearlySiteCovs", "obsCovs")) {
    data <- methods::slot(frame, nm)
    expected <- switch(nm, siteCovs = m, yearlySiteCovs = m * tt,
      obsCovs = m * tt * length(lengths))
    if (!is.null(data) && (nrow(data) != expected || anyDuplicated(names(data))))
      stop("GDR ", nm, " must have unique names and ", expected,
        " rows in native site/period order.", call. = FALSE)
  }
  # State abundance from an abundance-output fit refers to its fitted circle.
  # Do not silently reinterpret it as a count in a different survey footprint.
  if (!isTRUE(all.equal(breaks, model@data@dist.breaks, tolerance = 0, check.attributes = FALSE)) ||
      !identical(frame@unitsIn, model@data@unitsIn))
    stop("Prediction distance breaks and input units must match the fitted GDR survey footprint.", call. = FALSE)
  list(model = model, frame = frame, m = m, tt = tt, spatial = spatial)
}

.msr_gdr_evaluate <- function(object, newdata, kernel_inputs, join_by, include_fixed, re.form) {
  prep <- .msr_gdr_prepare(object, newdata, kernel_inputs, join_by, include_fixed)
  model <- prep$model
  frame <- prep$frame
  m <- prep$m; tt <- prep$tt
  covs <- .msr_gdr_covariates(frame)
  native <- function(type, n) {
    value <- unmarked::predict(model, type = type, newdata = covs[[type]],
      re.form = re.form, na.rm = FALSE, level = NULL)$Predicted
    if (length(value) != n) stop("Native GDR predictions changed the site/period row layout.", call. = FALSE)
    value
  }
  lambda <- native("lambda", m)
  if (any(!is.na(lambda) & (!is.finite(lambda) | lambda < 0)))
    stop("Native GDR intensities must be finite and nonnegative; check model fit and extrapolation.", call. = FALSE)
  phi <- if ("phi" %in% .msr_gdr_components(model)) native("phi", m * tt) else rep(1, m * tt)
  rem <- matrix(native("rem", m * tt * length(frame@period.lengths)),
                nrow = m * tt, byrow = TRUE)
  if (any(rem < 0 | rem > 1, na.rm = TRUE) || any(phi < 0 | phi > 1, na.rm = TRUE))
    stop("Native GDR removal and availability predictions must be probabilities.", call. = FALSE)
  rem <- .msr_gdr_removal(rem, frame@period.lengths)
  distance <- .msr_gdr_distance(model@keyfun,
    if (model@keyfun == "uniform") NULL else native("dist", m * tt),
    frame@dist.breaks, m * tt,
    if (model@keyfun == "hazard") exp(as.numeric(unmarked::coef(model, type = "scale"))) else NULL)
  availability <- .msr_gdr_probability(rowSums(rem), "Whole-survey availability")
  perceptibility <- .msr_gdr_probability(rowSums(distance), "Integrated perceptibility")
  area <- pi * max(frame@dist.breaks)^2 * if (frame@unitsIn == "m") 1e-6 else 1
  if (model@unitsOut == "ha") area <- area * 100
  base_abundance <- if (model@output == "density") lambda * area else lambda
  alpha <- if (model@mixture == "P") NULL else as.numeric(unmarked::coef(model, type = "alpha"))
  if (!model@mixture %in% c("P", "NB", "ZIP") ||
      (model@mixture != "P" && (length(alpha) != 1L || !is.finite(alpha))))
    stop("A supported GDR mixture and a finite mixture parameter are required.", call. = FALSE)
  if (model@mixture == "NB" && (!is.finite(exp(alpha)) || exp(alpha) <= 0))
    stop("The GDR negative-binomial size must be finite and positive.", call. = FALSE)
  zero <- if (model@mixture == "ZIP") stats::plogis(alpha) else 0
  abundance <- base_abundance * (1 - zero)
  detection <- phi * availability * perceptibility
  out <- data.frame(site = rep(seq_len(m), each = tt), period = rep(seq_len(tt), m),
    lambda = rep(lambda, each = tt), zero_probability = zero,
    abundance = rep(abundance, each = tt), density = rep(abundance / area, each = tt),
    phi = phi, availability = availability, perceptibility = perceptibility,
    detection = detection, expected_count = rep(abundance, each = tt) * detection)
  attr(out, "density_units") <- model@unitsOut
  attr(out, "survey_area") <- area
  attr(out, "lambda_output") <- model@output
  attr(out, "prediction_contract") <- .msr_prediction_contract("native", "gdistremoval",
    "expected_response", "component_specific_response_scales",
    population = "fixed_at_estimate", kernel = if (prep$spatial) "fixed_at_estimate" else "not_applicable",
    effects = list(response = if (.msr_native_zero_random_effects(model, list(re.form = re.form)))
      "zero_deviation_reference" else "native_conditional_effects"),
    conditioning = c("survey_protocol", "native_parameter_estimates"),
    omitted = c("population_parameter_uncertainty", "kernel_parameter_uncertainty", "observation_noise"))
  list(table = out, distance = distance, removal = rem, base_abundance = base_abundance,
    alpha = alpha, prep = prep)
}

#' Refit a Distance-Removal Model During Spatial Optimization
#'
#' Supply this function as \code{refit_fn} in \code{multiScale_optim()} to skip
#' intermediate coefficient Hessians. Returned fitted models retain the original
#' call's standard-error setting. The likelihood and coefficient starts are
#' unchanged.
#'
#' @param model A fitted \code{unmarkedFitGDR} to update.
#' @param data An \code{unmarkedFrameGDR} with the updated spatial site covariates.
#' @param context A list supplied by \code{multiScale_optim()}. Only an explicit
#'   \code{return_model = FALSE} disables coefficient standard errors. The default
#'   \code{list(return_model = TRUE)} performs the complete native refit.
#' @return An updated \code{unmarkedFitGDR}. For scalar objective evaluations,
#'   coefficient covariance is omitted. When a fitted model is requested,
#'   native coefficient covariance follows the original model's \code{se}
#'   setting. Use \code{unmarked::coef()} and \code{unmarked::logLik()} as usual.
#' @details
#' This is an optional callback, with no change to the default refit policy.
#' It keeps the native engine, formulas, mixture, detection key, survey protocol,
#' abundance cutoff, optimizer controls, and starting values in the model call.
#' An original \code{se = FALSE} stays false for the returned model. Omitting
#' intermediate Hessians does not quantify uncertainty in spatial scales.
#'
#' The function contains no references to local data objects and can be
#' serialized for PSOCK workers. See
#' \code{vignette("distance-removal", package = "multiScaleR")} for its use in
#' the spatial workflow.
#' @examples
#' y <- matrix(c(2, 1, 0, 1, 1, 0, 3, 1, 1, 2, 0, 0), ncol = 3, byrow = TRUE)
#' frame <- unmarked::unmarkedFrameGDR(y, y, dist.breaks = c(0, 25, 50, 100),
#'   unitsIn = "m")
#' model <- unmarked::gdistremoval(data = frame, K = 25, se = FALSE)
#' objective_fit <- refit_gdistremoval(model, frame, list(return_model = FALSE))
#' as.numeric(unmarked::logLik(objective_fit))  # Native likelihood at the updated data
#' @seealso \code{\link{multiScale_optim}}, \code{\link{predict_gdistremoval}}
#' @export
refit_gdistremoval <- function(model, data, context = list(return_model = TRUE)) {
  if (!methods::is(model, "unmarkedFitGDR") || !methods::is(data, "unmarkedFrameGDR"))
    stop("Supply an unmarkedFitGDR model and unmarkedFrameGDR data.", call. = FALSE)
  if (!is.list(context)) stop("context must be the refit context list.", call. = FALSE)
  if (identical(context$return_model, FALSE))
    unmarked::update(model, data = data, se = FALSE)
  else
    unmarked::update(model, data = data)
}

#' Predict Combined Distance and Removal Processes
#'
#' Calculate detection-adjusted abundance and density, whole-survey detection,
#' and expected observed counts from a fitted \code{unmarked::gdistremoval()}
#' model. Native coefficients and fitted spatial scales remain fixed. The
#' function integrates the distance key separately for every site and visit.
#'
#' @param object A fitted \code{unmarkedFitGDR}, or a \code{multiScaleR} object
#'   containing that native model. Use the native-refit backend for fitting.
#' @param newdata An \code{unmarkedFrameGDR} describing prediction sites,
#'   visits, covariates, and removal durations. The default \code{NULL} uses
#'   the fitted frame. For new sites, supply zero-filled response matrices with
#'   equal totals for distance and removal. Distance breaks must start at zero
#'   and, with input units, match the fitted survey footprint. Removal durations
#'   must be positive integers and may change.
#' @param kernel_inputs Optional \code{kernel_prep()} inputs for the prediction
#'   sites. With a \code{multiScaleR} fit and explicit \code{newdata}, these
#'   recalculate spatial site covariates using the fitted transforms and training
#'   normalization. Without these inputs, frame covariates are passed unchanged.
#' @param join_by Optional data frame of unique prediction-site keys, required
#'   with \code{kernel_inputs}. Keys follow prepared spatial order and also
#'   occur in the prediction frame's site covariates.
#' @param include_fixed Logical. Default \code{FALSE}. With spatial inputs,
#'   set \code{TRUE} to reconstruct and verify modeled fixed-radius metrics.
#' @param re.form Native unmarked random-effect policy. The default \code{NULL}
#'   conditions on fitted group effects. Set \code{NA} for zero group deviations
#'   at new sites. This does not integrate over the random-effect distribution.
#'
#' @return A data frame in site-major, then visit order, with columns:
#'   \describe{
#'     \item{site, period}{One-based frame row and primary-period indices.}
#'     \item{lambda}{Native abundance or density intensity, before ZIP adjustment.}
#'     \item{zero_probability}{Structural-zero probability, zero for P and NB.}
#'     \item{abundance}{Expected latent abundance within the survey circle,
#'       including ZIP adjustment. It is a superpopulation mean when phi is modeled.}
#'     \item{density}{Abundance divided by surveyed area in the fitted output units.}
#'     \item{phi}{Between-visit availability, or one for a single-visit fit.}
#'     \item{availability}{Probability of giving a detectable cue during at least
#'       one removal interval, conditional on presence.}
#'     \item{perceptibility}{Distance-integrated detection given cue availability.}
#'     \item{detection}{Product of phi, availability, and perceptibility.}
#'     \item{expected_count}{Expected observed individuals per survey.}
#'   }
#'   Attributes \code{density_units}, \code{survey_area}, and \code{lambda_output}
#'   identify the units and native intensity definition. The
#'   \code{prediction_contract} attribute records conditioning and omitted
#'   uncertainty. These are point estimates, without SEs or confidence limits.
#' @details
#' Time and distance are independent conditional on modeled covariates. Removal
#' probabilities refer to one time unit; unequal integer durations are applied
#' before combining intervals. Uniform distance keys need no distance parameter.
#' For ZIP, the structural-zero probability uses the fitted parameter's logit
#' link, rather than treating the Poisson-component intensity as the full mean.
#'
#' This helper does not change unmarked's likelihood or methods. It calculates
#' its own integrated probabilities, including visit-specific distance scales,
#' and avoids division by a zero interval-removal probability. Static spatial
#' predictors may enter any active component. Their environmental kernel widths
#' describe landscape effects and are distinct from the distance detection key.
#' Neither this function nor \code{simulate_gdistremoval()} propagates coefficient,
#' model-selection, spatial-scale, or random-effect variance uncertainty.
#'
#' For a statewide density map, native \code{predict(fit, ..., type = "lambda")}
#' and raster projection remain available. Native lambda is a count-component
#' intensity for ZIP. Use this helper on grid-site frames when an unconditional
#' density or a complete survey scenario is required. See
#' \code{vignette("distance-removal", package = "multiScaleR")}.
#' @examples
#' y <- matrix(c(2, 1, 0, 1, 1, 0, 3, 1, 1, 2, 0, 0), ncol = 3, byrow = TRUE)
#' frame <- unmarked::unmarkedFrameGDR(y, y, dist.breaks = c(0, 25, 50, 100),
#'   unitsIn = "m", period.lengths = c(1, 2, 2))
#' model <- unmarked::gdistremoval(data = frame, K = 25, se = FALSE)
#' prediction <- predict_gdistremoval(model)
#' prediction[, c("density", "detection", "expected_count")]
#' # Density corrects for detection; expected_count describes the survey.
#' @seealso \code{\link{simulate_gdistremoval}},
#'   \code{\link{refit_gdistremoval}}, \code{\link{prepare_unmarked_data}},
#'   \code{\link{multiScale_optim}}
#' @export
predict_gdistremoval <- function(object, newdata = NULL, kernel_inputs = NULL,
                                 join_by = NULL, include_fixed = FALSE, re.form = NULL) {
  .msr_gdr_evaluate(object, newdata, kernel_inputs, join_by, include_fixed, re.form)$table
}

#' Simulate Combined Distance and Removal Counts
#'
#' Draw paired distance and removal histories conditional on fitted parameters,
#' using visit-specific integrated detection and the fitted abundance mixture.
#' Each visit has identical distance and removal totals. Use these replicates
#' to investigate model behavior or develop conditional predictive checks.
#'
#' @inheritParams predict_gdistremoval
#' @param nsim Positive integer number of replicate frames. Default \code{1}.
#' @return A list of \code{nsim} \code{unmarkedFrameGDR} objects. Covariates and
#'   protocol are retained from the prepared prediction frame. Counts are random
#'   draws, and visits with an incomplete distance or removal history remain
#'   missing in both histories. A \code{prediction_contract} attribute states
#'   that fitted parameters and spatial scales are fixed, with observation and
#'   latent-abundance variation simulated.
#' @details
#' The latent abundance is drawn once per site and shared across visits. P uses
#'   Poisson abundance, NB uses the native size parameter, and ZIP draws structural
#'   absence with its fitted logit probability before Poisson abundance. Simulated
#'   totals are allocated separately to distance and first-detection-time classes
#'   under the conditional-independence assumption.
#'
#' Set the random seed before calling this function for reproducibility. Group
#'   effects follow \code{re.form}; new random effects are not sampled. These
#'   conditional simulations do not by themselves provide a calibrated bootstrap
#'   for random-effect variances, selected scales, or model selection. Such an
#'   analysis must reproduce the entire fitting and sampling design.
#' @examples
#' y <- matrix(c(2, 1, 0, 1, 1, 0, 3, 1, 1, 2, 0, 0), ncol = 3, byrow = TRUE)
#' frame <- unmarked::unmarkedFrameGDR(y, y, dist.breaks = c(0, 25, 50, 100), unitsIn = "m")
#' model <- unmarked::gdistremoval(data = frame, K = 25, se = FALSE)
#' set.seed(14)
#' replicate <- simulate_gdistremoval(model)[[1]]
#' # These two matrices describe the same simulated individuals.
#' rowSums(replicate@yDistance) == rowSums(replicate@yRemoval)
#' @seealso \code{\link{predict_gdistremoval}}
#' @export
simulate_gdistremoval <- function(object, nsim = 1, newdata = NULL, kernel_inputs = NULL,
                                  join_by = NULL, include_fixed = FALSE, re.form = NULL) {
  validate_scalar_numeric(nsim, "nsim", positive = TRUE, integerish = TRUE)
  evaluated <- .msr_gdr_evaluate(object, newdata, kernel_inputs, join_by, include_fixed, re.form)
  prep <- evaluated$prep
  frame <- prep$frame
  m <- prep$m; tt <- prep$tt
  jd <- ncol(frame@yDistance) / tt
  jr <- ncol(frame@yRemoval) / tt
  observed <- matrix(TRUE, nrow = m, ncol = tt)
  for (t in seq_len(tt)) observed[, t] <-
    rowSums(is.na(frame@yDistance[, (t - 1L) * jd + seq_len(jd), drop = FALSE])) == 0L &
    rowSums(is.na(frame@yRemoval[, (t - 1L) * jr + seq_len(jr), drop = FALSE])) == 0L
  valid <- !is.na(evaluated$base_abundance)
  out <- lapply(seq_len(nsim), function(i) {
    abundance <- rep(NA_real_, m)
    abundance[valid] <- if (prep$model@mixture == "NB") stats::rnbinom(sum(valid),
      size = exp(evaluated$alpha), mu = evaluated$base_abundance[valid]) else
        stats::rpois(sum(valid), evaluated$base_abundance[valid])
    if (prep$model@mixture == "ZIP") abundance[valid] <- abundance[valid] *
      stats::rbinom(sum(valid), 1, 1 - stats::plogis(evaluated$alpha))
    yd <- matrix(NA_real_, m, tt * jd)
    yr <- matrix(NA_real_, m, tt * jr)
    for (site in seq_len(m)) for (t in seq_len(tt)) {
      row <- (site - 1L) * tt + t
      p <- evaluated$table$detection[row]
      if (!observed[site, t] || is.na(abundance[site]) || is.na(p)) next
      count <- stats::rbinom(1L, abundance[site], p)
      yd[site, (t - 1L) * jd + seq_len(jd)] <- if (count == 0) rep(0, jd) else
        as.vector(stats::rmultinom(1L, count, evaluated$distance[row, ]))
      yr[site, (t - 1L) * jr + seq_len(jr)] <- if (count == 0) rep(0, jr) else
        as.vector(stats::rmultinom(1L, count, evaluated$removal[row, ]))
    }
    replicate <- frame
    dimnames(yd) <- dimnames(frame@yDistance)
    dimnames(yr) <- dimnames(frame@yRemoval)
    replicate@yDistance <- yd
    replicate@yRemoval <- replicate@y <- yr
    replicate
  })
  contract <- attr(evaluated$table, "prediction_contract")
  contract$quantity <- "replicate_response"
  contract$uncertainty$observation_noise <- TRUE
  contract$uncertainty$omitted <- setdiff(contract$uncertainty$omitted, "observation_noise")
  attr(out, "prediction_contract") <- contract
  out
}
