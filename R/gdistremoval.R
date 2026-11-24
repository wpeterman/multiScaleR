# Whole-survey predictions and simulations for unmarked::gdistremoval() fits.
#
# These helpers combine the fitted components the way the native likelihood
# does (unmarked 1.5.2, src/nll_gdistremoval.cpp): abundance or density times
# the survey area, phi for between-visit availability, the removal cell
# probabilities for each interval, and the distance-bin probabilities of the
# point-transect key. They do not call unmarked::getP(), whose GDR method in
# unmarked 1.5.2 integrates every primary period with the first period's
# distance parameter.

# Active native components, in the order the native likelihood uses them.
.msr_gdr_components <- function(model) names(.msr_gdr_active_formulas(model))

# Covariate tables for each component, with the native precedence: removal
# rows see observation, yearly and site covariates; phi and distance rows see
# yearly and site covariates; lambda rows see site covariates. Each component
# is predicted from its own table, which keeps missing rows in place.
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

# Probability that an individual in the survey circle is in each distance bin
# and is detected there, given that it is available. Rows follow the
# site-by-period order of `parameter` (the native distance parameter: sigma
# for "halfnorm", the decay length for "exp", the distance scale for "hazard").
# All integrals use x = r / radius and the point-transect measure 2 x dx, so a
# uniform key gives each bin its share of the circle's area.
.msr_gdr_distance <- function(key, parameter, breaks, n, hazard = NULL) {
  radius <- max(breaks)
  lower <- utils::head(breaks, -1L) / radius
  upper <- utils::tail(breaks, -1L) / radius
  if (key == "uniform") return(matrix(rep(upper^2 - lower^2, each = n), nrow = n))
  if (length(parameter) != n || any(parameter < 0, na.rm = TRUE))
    stop("GDR distance predictions must be nonnegative and follow site/period order.", call. = FALSE)
  if (key == "hazard" && (length(hazard) != 1L || !is.finite(hazard) || hazard <= 0))
    stop("The hazard-rate exponent must be finite and positive.", call. = FALSE)
  # Evaluate each distinct parameter value once. The half-normal and
  # exponential keys have closed forms evaluated for all values together;
  # only the hazard-rate key needs numerical quadrature.
  unique_par <- unique(parameter[!is.na(parameter)])
  probabilities <- matrix(0, length(unique_par), length(lower))
  ratio <- unique_par / radius
  # An infinite scale, or one too large for floating point, detects everything.
  flat <- is.infinite(ratio)
  if (key == "halfnorm") flat <- flat | (ratio > 0 & 0.5 / ratio^2 == 0)
  if (any(flat)) probabilities[flat, ] <- matrix(rep(upper^2 - lower^2,
    each = sum(flat)), nrow = sum(flat))
  # A zero scale detects nothing, so its row stays zero.
  active <- which(ratio > 0 & !flat)
  if (key == "halfnorm") {
    # Integral of 2 x exp(-a x^2) from l to u is
    # exp(-a l^2) (1 - exp(-a (u^2 - l^2))) / a, with a = 1 / (2 ratio^2).
    a <- 0.5 / ratio[active]^2
    # Scales below floating-point range keep their zero-probability limit.
    active <- active[is.finite(a)]
    a <- a[is.finite(a)]
    for (j in seq_along(lower)) probabilities[active, j] <-
      exp(-a * lower[j]^2) * (-expm1(-a * (upper[j]^2 - lower[j]^2))) / a
  } else if (key == "exp") {
    # Integral of 2 x exp(-x / s) is 2 s^2 times a shape-two gamma
    # probability. Lower-tail differences lose precision far into the tail,
    # so bins beyond the scale use upper-tail differences instead.
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
    # Hazard-rate key 1 - exp(-(r / scale)^(-b)), b the fitted exponent.
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

# Removal cell probabilities: the chance that an available individual is
# first detected in each interval. `p` holds per-time-unit probabilities and
# `lengths` the integer interval durations, as in unmarked's removalPiFun().
# Tracking the surviving fraction avoids dividing by a zero probability.
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

# Correct only integration and rounding excursions, never an invalid fit.
.msr_gdr_probability <- function(x, label) {
  if (any(!is.na(x) & (!is.finite(x) | x < -1e-10 | x > 1 + 1e-10)))
    stop(label, " must be a probability. Check the fitted detection model.", call. = FALSE)
  pmin(1, pmax(0, x))
}

# Check the model and prediction frame, and rebuild spatial site covariates
# from kernel inputs when they are supplied.
.msr_gdr_prepare <- function(object, newdata, kernel_inputs, join_by, include_fixed) {
  validate_scalar_logical(include_fixed, "include_fixed")
  spatial <- inherits(object, "multiScaleR")
  model <- if (spatial) object$opt_mod else object
  if (!methods::is(model, "unmarkedFitGDR"))
    stop("object must be an unmarkedFitGDR or a multiScaleR fit containing one.", call. = FALSE)
  if (methods::is(model, "unmarkedFitGDRField"))
    stop("GDR response fields carry a spatial abundance field; use predict() for them.", call. = FALSE)
  # Zero-inflated TMB fits and active offsets are not qualified (see
  # .msr_gdr_native_guard()), so their combined means are not reported.
  .msr_gdr_native_guard(model)
  frame <- if (is.null(newdata)) model@data else newdata
  if (!methods::is(frame, "unmarkedFrameGDR"))
    stop("newdata must be an explicit unmarkedFrameGDR with the prediction covariates and survey protocol.", call. = FALSE)
  if (!is.null(kernel_inputs)) {
    if (!spatial || is.null(newdata))
      stop("kernel_inputs require a multiScaleR fit and an explicit prediction frame.", call. = FALSE)
    covs <- intersect(rownames(object$scale_est), .unmarked_model_predictors(model))
    shadow <- intersect(covs, union(names(frame@yearlySiteCovs), names(frame@obsCovs)))
    if (length(shadow))
      stop("Optimized GDR site covariates are shadowed in yearlySiteCovs or obsCovs: ",
        paste(shadow, collapse = ", "), ". Use distinct names for distinct native inputs.",
        call. = FALSE)
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
  # Abundance from an abundance-output fit refers to its fitted circle. Do not
  # silently reinterpret it as a count in a different survey footprint.
  if (!isTRUE(all.equal(breaks, model@data@dist.breaks, tolerance = 0, check.attributes = FALSE)) ||
      !identical(frame@unitsIn, model@data@unitsIn))
    stop("Prediction distance breaks and input units must match the fitted GDR survey footprint.", call. = FALSE)
  list(model = model, frame = frame, m = m, tt = tt, spatial = spatial)
}

# Evaluate every component and combine them into whole-survey quantities.
.msr_gdr_evaluate <- function(object, newdata, kernel_inputs, join_by, include_fixed, re.form) {
  prep <- .msr_gdr_prepare(object, newdata, kernel_inputs, join_by, include_fixed)
  model <- prep$model
  frame <- prep$frame
  m <- prep$m; tt <- prep$tt
  components <- .msr_gdr_components(model)
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
  phi <- if ("phi" %in% components) native("phi", m * tt) else rep(1, m * tt)
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
  # Survey-circle area in the output units (square km, times 100 for ha),
  # as the native likelihood scales a density-output lambda.
  area <- pi * max(frame@dist.breaks)^2 * if (frame@unitsIn == "m") 1e-6 else 1
  if (model@unitsOut == "ha") area <- area * 100
  base_abundance <- if (model@output == "density") lambda * area else lambda
  alpha <- if (model@mixture == "P") NULL else as.numeric(unmarked::coef(model, type = "alpha"))
  if (!model@mixture %in% c("P", "NB", "ZIP") ||
      (model@mixture != "P" && (length(alpha) != 1L || !is.finite(alpha))))
    stop("A supported GDR mixture and a finite mixture parameter are required.", call. = FALSE)
  # The native likelihood uses size exp(alpha) for NB and a structural-zero
  # probability plogis(alpha) for ZIP.
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

#' Predict Whole-Survey Detection and Expected Counts for Distance-Removal Models
#'
#' Combine the fitted components of an [unmarked::gdistremoval()] model into
#' the quantities a survey produces: expected abundance and density in the
#' survey circle, the probability that an individual present is detected at
#' each visit, and the expected number of individuals counted. Fitted
#' coefficients and kernel scales are held at their estimates.
#'
#' This function is experimental. Its arguments and results may change in a
#' later release without a deprecation period.
#'
#' @param object A fitted `unmarkedFitGDR` model, or a `multiScaleR` fit from
#'   [multiScale_optim()] that contains one.
#' @param newdata An `unmarkedFrameGDR` describing the prediction sites, their
#'   visits, covariates, and removal interval durations. The default `NULL`
#'   uses the fitted frame. For new sites, supply response matrices of zeros
#'   with the same layout. Distance breaks and their units must match the
#'   fitted survey footprint. Removal durations must be positive integers and
#'   may differ from the fitted ones, so you can predict a longer or shorter
#'   count.
#' @param kernel_inputs Optional [kernel_prep()] output for the prediction
#'   sites. With a `multiScaleR` fit and explicit `newdata`, the spatial site
#'   covariates are rebuilt at the fitted scales with the training
#'   normalization, as in [prepare_unmarked_data()]. Without it, the
#'   covariates in `newdata` are used as they are.
#' @param join_by Data frame of prediction-site keys, required with
#'   `kernel_inputs`. See [prepare_unmarked_data()].
#' @param include_fixed Logical. Default `FALSE`. With `kernel_inputs`, set
#'   `TRUE` to rebuild and check fixed-radius metrics in the model.
#' @param re.form Random-effect policy passed to `unmarked::predict()`. The
#'   default `NULL` uses the fitted group effects. `NA` sets group deviations
#'   to zero, for example at sites in unsampled groups. Neither averages over
#'   the random-effect distribution.
#'
#' @return A data frame with one row per site and visit (primary period), in
#'   site order and then visit order, with columns:
#'   * `site`, `period`: row number in the frame and visit number.
#'   * `lambda`: the native abundance or density parameter. For a
#'     zero-inflated (`"ZIP"`) model it is the Poisson-component mean, before
#'     structural zeros.
#'   * `zero_probability`: the probability of a structural zero, 0 for
#'     `"P"` and `"NB"` models.
#'   * `abundance`: expected number of individuals in the survey circle,
#'     including structural zeros. When phi is modeled it counts the
#'     individuals that use the site across visits, not those present at one
#'     visit.
#'   * `density`: `abundance` divided by the survey-circle area, in
#'     individuals per hectare or square kilometer (the fit's `unitsOut`).
#'   * `phi`: probability that an individual is available at this visit, 1
#'     for a single-visit model.
#'   * `availability`: probability that an available individual is detected
#'     (gives a cue) in at least one removal interval.
#'   * `perceptibility`: probability that an individual giving a cue anywhere
#'     in the circle is detected, from the distance detection function.
#'   * `detection`: `phi * availability * perceptibility`, the probability
#'     that an individual present is counted at this visit.
#'   * `expected_count`: `abundance * detection`, the expected number of
#'     individuals counted at this visit.
#'
#'   Attributes `density_units`, `survey_area` (in those units), and
#'   `lambda_output` record the units and whether the native lambda is an
#'   abundance or a density. The `prediction_contract` attribute records what
#'   the predictions condition on (see [prediction_contract]). These are
#'   point estimates without standard errors or intervals.
#'
#' @details
#' The components combine as they do in the native likelihood. The removal
#' probability refers to one time unit, so an interval lasting `d` units has
#' a first-detection probability of `1 - (1 - p)^d`. Distance and removal
#' detection are independent given the covariates. Each visit uses its own
#' distance parameter. The half-normal and exponential keys use closed-form
#' integrals over the distance bins, and the hazard-rate key uses numerical
#' integration.
#'
#' `unmarked::getP()` and functions built on it, such as `fitted()` and
#' `simulate()`, are not used. In unmarked 1.5.2 they apply the first
#' visit's distance parameter to every visit, which matters when the
#' distance formula contains a visit-level (yearly site) covariate.
#'
#' For a density map, `predict(fit, raster, type = "lambda")` remains the
#' route (see [predict.multiScaleR()]). For a zero-inflated model that map
#' shows the Poisson-component mean. Use this function on a frame of grid
#' sites when you need the mean including structural zeros, or the expected
#' count under a planned survey protocol.
#'
#' Zero-inflated models fitted with `engine = "TMB"`, models with active
#' offsets, and GDR response fields from [fit_gdr_response_field()] are not
#' supported.
#'
#' @examples
#' y <- matrix(c(2, 1, 0, 1, 1, 0, 3, 1, 1, 2, 0, 0), ncol = 3, byrow = TRUE)
#' frame <- unmarked::unmarkedFrameGDR(y, y, dist.breaks = c(0, 25, 50, 100),
#'   unitsIn = "m", period.lengths = c(1, 2, 2))
#' model <- unmarked::gdistremoval(data = frame, K = 25, se = FALSE)
#' prediction <- predict_gdistremoval(model)
#' prediction[, c("density", "detection", "expected_count")]
#' # density corrects for detection; expected_count is what a survey records
#' @seealso [simulate_gdistremoval()], [prepare_unmarked_data()],
#'   [multiScale_optim()], [sim_dat_gdr()]
#' @export
predict_gdistremoval <- function(object, newdata = NULL, kernel_inputs = NULL,
                                 join_by = NULL, include_fixed = FALSE, re.form = NULL) {
  .msr_gdr_evaluate(object, newdata, kernel_inputs, join_by, include_fixed, re.form)$table
}

#' Simulate Distance-Removal Surveys from a Fitted Model
#'
#' Draw new distance and removal count histories from a fitted
#' [unmarked::gdistremoval()] model, holding its coefficients and kernel
#' scales at their estimates. Each visit uses its own detection
#' probabilities, and the distance and removal histories of a visit describe
#' the same simulated individuals. Use the replicates to study how the model
#' behaves or to build conditional predictive checks.
#'
#' This function is experimental. Its arguments and results may change in a
#' later release without a deprecation period.
#'
#' @inheritParams predict_gdistremoval
#' @param nsim Positive whole number of replicate frames. Default `1`.
#' @return A list of `nsim` `unmarkedFrameGDR` objects with the covariates and
#'   protocol of the prediction frame and simulated counts. A visit with a
#'   missing distance or removal history in the frame stays missing in both.
#'   The `prediction_contract` attribute records that coefficients and scales
#'   are fixed and that latent abundance and observation noise are simulated.
#'
#' @details
#' Latent abundance is drawn once per site and shared by its visits: Poisson
#' for `"P"`, negative binomial with the fitted size for `"NB"`, and a
#' structural zero with the fitted probability followed by a Poisson draw
#' for `"ZIP"`. The number counted at each visit is binomial with the
#' `detection` probability from [predict_gdistremoval()]. That total is then
#' split independently among the distance bins and among the removal
#' intervals with the fitted cell probabilities.
#'
#' Call [set.seed()] first for reproducible draws. Group effects follow
#' `re.form`; new group effects are not drawn. These conditional replicates
#' alone are not a calibrated bootstrap of random-effect variances, selected
#' scales, or model selection, which would need the whole fitting and
#' sampling process repeated.
#'
#' `unmarked::simulate()` is not used. In unmarked 1.5.2 it applies the
#' first visit's distance parameter to every visit and back-transforms the
#' zero-inflation parameter with `exp()`, whereas the likelihood uses
#' `plogis()`.
#'
#' @examples
#' y <- matrix(c(2, 1, 0, 1, 1, 0, 3, 1, 1, 2, 0, 0), ncol = 3, byrow = TRUE)
#' frame <- unmarked::unmarkedFrameGDR(y, y, dist.breaks = c(0, 25, 50, 100),
#'   unitsIn = "m")
#' model <- unmarked::gdistremoval(data = frame, K = 25, se = FALSE)
#' set.seed(14)
#' replicate <- simulate_gdistremoval(model)[[1]]
#' # Both histories count the same simulated individuals
#' rowSums(replicate@yDistance) == rowSums(replicate@yRemoval)
#' @seealso [predict_gdistremoval()], [sim_dat_gdr()]
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
  # A visit is simulated only where both of its histories are complete.
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
