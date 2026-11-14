#' Simulate Paired Distance-Removal Surveys with Retained Truth
#'
#' Simulate single-visit point counts in which each bird is recorded by
#' distance band and by the time interval of its first detection. Abundance
#' is Poisson, detection by distance follows a half-normal curve, and
#' detection over time follows a removal (first-detection) process. The two
#' count tables describe the same birds. The observed data and the true
#' values used to generate them are returned separately, so you can fit a
#' model to the data and compare its estimates with the truth.
#'
#' @param density Numeric vector of expected bird densities, in birds per
#'   square kilometer, one value per survey point. Zero is allowed. This is
#'   the response scale, not the log scale. Compute any ecological formula
#'   or kernel-weighted covariate effect first, for example
#'   `exp(log(400) + 0.5 * habitat)`.
#' @param distance_scale Positive number, or one per point. The standard
#'   deviation (sigma) of the half-normal detection function, in `units_in`.
#'   Detection falls to about 61% of its maximum at this distance. It is
#'   unrelated to the kernel scale of effect.
#' @param removal Probability, between 0 and 1, that a bird present is first
#'   detected during one unit of interval time. Supply one number, one value
#'   per point, or a matrix with points in rows and removal intervals in
#'   columns for interval-specific probabilities.
#' @param dist_breaks Increasing numeric distance boundaries starting at
#'   zero, in `units_in`. The last boundary is the survey radius.
#' @param interval_lengths Positive whole numbers giving the duration of each
#'   removal interval, in `interval_units`. Intervals may differ in length.
#' @param interval_units Character string naming the time unit of
#'   `interval_lengths`, for example `"minutes"`. It is recorded, not
#'   converted.
#' @param units_in Distance units, `"m"` (default) or `"km"`.
#' @param site_covs Optional data frame of site covariates, one row per point,
#'   such as kernel-weighted landscape covariates from [kernel_prep()].
#'   Columns are copied into the observed frame unchanged. As in `unmarked`,
#'   character columns are converted to factors with a warning.
#' @param event_keys Optional data frame with one row per point holding the
#'   original survey identifiers. It is returned unchanged, separate from the
#'   covariates, for use with [validate_unmarked_data()].
#' @param seed Nonnegative whole number used as the random seed. Your random
#'   number state is restored afterward, so the simulation does not change
#'   later random draws in your session.
#' @param return_joint Logical, default `FALSE`. If `TRUE`, also return the
#'   simulated counts by point, distance band and interval as `joint`.
#' @param provenance Optional named list of notes about how the inputs were
#'   made, such as the generating formula, coefficients or the kernel scale
#'   used for a covariate. It is stored in `design$provenance` without being
#'   evaluated.
#' @param x An object returned by `sim_dat_gdr()`.
#' @param ... Not used.
#' @return A list of class `msr_gdr_simulation`:
#'   * `observed`: the `unmarkedFrameGDR` to analyze, with the distance and
#'     removal tables, `site_covs`, distance breaks, units and interval
#'     lengths. It holds no true values.
#'   * `truth`: a data frame with one row per point. See "Read the truth".
#'   * `probabilities`: a list with `distance`, a points by distance-band
#'     matrix of the probability that a bird in the circle is detected in each
#'     band, and `removal`, a points by interval matrix of the probability that
#'     a bird is first detected in each interval.
#'   * `joint`: the points by distance band by interval array of simulated
#'     counts when `return_joint = TRUE`, otherwise `NULL`.
#'   * `event_keys`: as supplied.
#'   * `design`: the settings used, including `density_units`
#'     (`"birds/km2"`), `area_km2` (the survey circle), `distance_scale` and
#'     `removal` per point, the breaks, interval lengths and units, the `seed`
#'     and `RNGkind`, the model assumptions (`mixture = "P"`,
#'     `keyfun = "halfnorm"`, one primary period, closed population), the
#'     `unmarked` version and `provenance`.
#'
#'   Printing shows a short summary of the number of points, bands and
#'   intervals, the density range and the units, and returns the object
#'   invisibly. The full frame and truth stay in the object.
#' @details
#' # How the counts are generated
#'
#' The survey area is `pi * radius^2`, in square kilometers, and expected
#' abundance is density times that area. Detection by distance and by time
#' are independent for each bird. The probability of a first detection in an
#' interval accounts for not being detected in earlier intervals, and an
#' interval of duration `t` has detection probability `1 - (1 - p)^t` for
#' per-unit probability `p`. The function draws independent Poisson counts
#' for every combination of distance band and interval, plus one for birds
#' never detected. Their sum is the Poisson abundance in the circle, and both
#' count tables are totals of the same simulated counts, so their row totals
#' always match.
#'
#' # Read the truth
#'
#' | Quantity | Where | How to read it |
#' |:---|:---|:---|
#' | Expected density | `truth$expected_density` | Birds per square km, before imperfect detection |
#' | Expected abundance | `truth$expected_abundance` | Expected birds in the survey circle: density times `design$area_km2` |
#' | Realized abundance | `truth$realized_abundance` | Birds simulated in the circle, a Poisson draw around the expectation |
#' | Distance detection | `truth$distance_probability` | Chance that a bird in the circle is detected by distance, averaged over the circle |
#' | Removal detection | `truth$removal_probability` | Chance that a bird is first detected at some point during the count |
#' | Overall detection | `truth$observation_probability` | The product of the two, the share of birds an observer records |
#' | Expected observations | `truth$expected_observed` | Expected count after both detection processes |
#' | Realized observations | `truth$realized_observed` | Birds recorded, equal to the row totals of both count tables |
#' | Distance scale | `design$distance_scale` | Half-normal sigma, in `units_in` |
#' | Your notes | `design$provenance` | Whatever you passed in `provenance` |
#'
#' The truth never enters the observed frame, so it cannot leak into fitted
#' covariates. `truth$site_row` gives the frame row of each point.
#'
#' # What the simulator covers
#'
#' The simulator covers single-visit Poisson surveys with a half-normal
#' distance function and equal or unequal whole-number interval lengths. It
#' assumes a closed population inside the radius. It does not simulate
#' repeat visits with availability, negative binomial or zero-inflated
#' abundance, movement across the radius, or random variation in detection.
#' It fits no model and leaves [sim_dat_unmarked()] unchanged. See the
#' "Simulate a survey with a known scale of effect" section of
#' `vignette("v12-distance-removal", package = "multiScaleR")` for a full
#' example.
#' @seealso [validate_unmarked_data()], [prepare_gdr_histories()],
#'   [fit_gdr_response_field()]
#' @examples
#' # Three points: two with birds and one with density zero
#' sim <- sim_dat_gdr(density = c(100, 200, 0), distance_scale = 40,
#'   removal = 0.3, dist_breaks = c(0, 25, 75),
#'   interval_lengths = c(1, 2, 1), interval_units = "minutes", seed = 24)
#' sim                        # compact summary
#' # Truth behind the counts: about 1.8 and 3.5 birds expected in the circle
#' sim$truth[, c("expected_density", "expected_abundance",
#'               "realized_abundance", "realized_observed")]
#' # The same birds by distance band and by time interval
#' sim$observed@yDistance
#' sim$observed@yRemoval
#' @export
sim_dat_gdr <- function(density, distance_scale, removal, dist_breaks,
                        interval_lengths, interval_units, units_in = c("m", "km"),
                        site_covs = NULL, event_keys = NULL, seed,
                        return_joint = FALSE, provenance = list()) {
  if (!is.numeric(density) || !length(density) || !is.null(dim(density)) ||
      any(!is.finite(density) | density < 0))
    stop("density must be a finite nonnegative vector in birds per square kilometer.", call. = FALSE)
  n <- length(density)
  expand <- function(x, label) {
    if (!is.numeric(x) || !is.null(dim(x)) || !length(x) ||
        !length(x) %in% c(1L, n) || any(!is.finite(x)))
      stop(label, " must be a finite scalar or a vector with one value per site.", call. = FALSE)
    rep(x, length.out = n)
  }
  sigma <- expand(distance_scale, "distance_scale")
  if (any(sigma <= 0)) stop("distance_scale must be positive.", call. = FALSE)
  units_in <- match.arg(units_in)
  if (!is.numeric(dist_breaks) || !is.null(dim(dist_breaks)) || length(dist_breaks) < 2L ||
      any(!is.finite(dist_breaks)) || dist_breaks[1L] != 0 || any(diff(dist_breaks) <= 0))
    stop("dist_breaks must start at zero and increase strictly.", call. = FALSE)
  if (!is.numeric(interval_lengths) || !is.null(dim(interval_lengths)) || !length(interval_lengths) ||
      any(!is.finite(interval_lengths) | interval_lengths <= 0 | interval_lengths != round(interval_lengths)))
    stop("interval_lengths must contain positive integer durations.", call. = FALSE)
  if (!is.character(interval_units) || length(interval_units) != 1L ||
      is.na(interval_units) || !nzchar(trimws(interval_units)))
    stop("Supply an explicit, nonempty interval_units label.", call. = FALSE)
  if (!is.logical(return_joint) || length(return_joint) != 1L || is.na(return_joint))
    stop("return_joint must be one nonmissing logical value.", call. = FALSE)
  if (missing(seed) || !is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != round(seed))
    stop("Supply a nonnegative integer seed.", call. = FALSE)
  if (!is.list(provenance) || (length(provenance) && (is.null(names(provenance)) ||
      anyNA(names(provenance)) || any(!nzchar(names(provenance))) || anyDuplicated(names(provenance)))))
    stop("provenance must be a named list.", call. = FALSE)
  for (nm in c("site_covs", "event_keys")) {
    x <- if (nm == "site_covs") site_covs else event_keys
    if (!is.null(x) && (!is.data.frame(x) || nrow(x) != n || anyDuplicated(names(x))))
      stop(nm, " must have one row per site and unique column names.", call. = FALSE)
  }
  jd <- length(dist_breaks) - 1L; jr <- length(interval_lengths)
  if (is.matrix(removal)) {
    if (!is.numeric(removal) || !identical(dim(removal), c(as.integer(n), as.integer(jr))) ||
        any(!is.finite(removal)))
      stop("removal matrix dimensions must be sites by removal intervals, with finite probabilities.", call. = FALSE)
    p <- removal
  } else p <- matrix(expand(removal, "removal"), n, jr)
  if (any(p < 0 | p > 1)) stop("removal probabilities must lie between zero and one.", call. = FALSE)
  radius <- utils::tail(dist_breaks, 1L)
  area <- pi * (radius / if (units_in == "m") 1000 else 1)^2
  expected_n <- density * area
  if (!is.finite(area) || area <= 0 || any(!is.finite(expected_n)))
    stop("Survey area and expected abundance must be finite and representable.", call. = FALSE)
  # Stable half-normal annulus integration on the dimensionless radius scale.
  a <- 0.5 * (radius / sigma)^2
  lower <- utils::head(dist_breaks, -1L) / radius
  upper <- utils::tail(dist_breaks, -1L) / radius
  pd <- matrix(0, n, jd)
  active <- is.finite(a)
  for (j in seq_len(jd)) {
    width <- upper[j]^2 - lower[j]^2
    x <- a[active] * width
    relative <- rep(1, length(x))
    positive <- x > 0
    relative[positive] <- -expm1(-x[positive]) / x[positive]
    pd[active, j] <- exp(-a[active] * lower[j]^2) * width * relative
  }
  pr <- matrix(0, n, jr); survival <- rep(1, n)
  for (j in seq_len(jr)) {
    interval <- -expm1(interval_lengths[j] * log1p(-p[, j]))
    pr[, j] <- survival * interval
    survival <- survival * (1 - interval)
  }
  joint_p <- matrix(0, n, jd * jr)
  for (j in seq_len(jr)) joint_p[, (j - 1L) * jd + seq_len(jd)] <- pd * pr[, j]
  observation_p <- rowSums(joint_p)
  if (any(!is.finite(observation_p) | observation_p < 0 | observation_p > 1 + 1e-12))
    stop("Computed observation probabilities are invalid.", call. = FALSE)
  # Correct only floating-point excursions at the probability boundary.
  unseen_p <- pmax(0, 1 - observation_p)
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) old_rng <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_rng) assign(".Random.seed", old_rng, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  rng_kind <- RNGkind()
  set.seed(as.integer(seed))
  counts <- matrix(stats::rpois(length(joint_p), as.vector(joint_p * expected_n)), n, jd * jr)
  unseen <- stats::rpois(n, expected_n * unseen_p)
  if (any(!is.finite(counts)) || any(!is.finite(unseen)))
    stop("Poisson allocation produced nonfinite counts; review expected abundance.", call. = FALSE)
  yd <- matrix(0, n, jd); yr <- matrix(0, n, jr)
  for (j in seq_len(jd)) yd[, j] <- rowSums(counts[, seq.int(j, jd * jr, by = jd), drop = FALSE])
  for (j in seq_len(jr)) yr[, j] <- rowSums(counts[, (j - 1L) * jd + seq_len(jd), drop = FALSE])
  observed <- unmarked::unmarkedFrameGDR(yDistance = yd, yRemoval = yr,
    numPrimary = 1L, siteCovs = site_covs, dist.breaks = dist_breaks,
    unitsIn = units_in, period.lengths = interval_lengths)
  realized_observed <- rowSums(counts)
  structure(list(observed = observed, truth = data.frame(site_row = seq_len(n),
    expected_density = density, expected_abundance = expected_n,
    realized_abundance = realized_observed + unseen,
    distance_probability = rowSums(pd), removal_probability = rowSums(pr),
    observation_probability = observation_p, expected_observed = expected_n * observation_p,
    realized_observed = realized_observed),
    probabilities = list(distance = pd, removal = pr),
    joint = if (return_joint) array(counts, c(n, jd, jr)) else NULL,
    event_keys = event_keys,
    design = list(density_units = "birds/km2", abundance_units = "birds/survey_circle",
      distance_units = units_in, interval_units = interval_units, area_km2 = area,
      distance_breaks = dist_breaks, interval_lengths = interval_lengths,
      distance_scale = sigma, removal = p, seed = as.integer(seed), RNGkind = rng_kind,
      mixture = "P", keyfun = "halfnorm", survey = "point", primary_periods = 1L,
      closed_population = TRUE, unmarked_version = as.character(utils::packageVersion("unmarked")),
      provenance = provenance)), class = "msr_gdr_simulation")
}

#' @rdname sim_dat_gdr
#' @export
print.msr_gdr_simulation <- function(x, ...) {
  cat("Single-primary Poisson distance-removal simulation\n")
  cat(nrow(x$truth), "survey events;", ncol(x$observed@yDistance),
      "distance bins;", ncol(x$observed@yRemoval), "removal intervals\n")
  cat("Expected density:", paste(format(range(x$truth$expected_density), digits = 5), collapse = " to "),
      x$design$density_units, "\n")
  cat("Distance units:", x$design$distance_units,
      "; removal time units:", x$design$interval_units, "\n")
  cat("Full observations: $observed; generating truth: $truth\n")
  invisible(x)
}
