#' Simulate Paired Distance-Removal Surveys with Retained Truth
#'
#' Generate single-primary point surveys with Poisson abundance, a half-normal
#' distance key, and first-detection removal intervals. Paired margins describe
#' the same birds. Observed data and generating truth are returned separately.
#'
#' @param density Numeric vector of expected bird densities, in birds per square
#'   kilometer, one value per survey site. Zero density is allowed. Evaluate
#'   ecological formulas or prepared kernel covariates before supplying this
#'   vector; density is on the response scale, not the log scale.
#' @param distance_scale Positive numeric scalar or vector with one value per
#'   site. This is the half-normal detection standard deviation, in `units_in`.
#'   It is distinct from an ecological kernel's scale of effect.
#' @param removal Numeric per-time-unit singing/detection probabilities in
#'   `[0, 1]`. Supply a scalar, a vector with one value per site, or a matrix
#'   with sites in rows and removal intervals in columns. Matrix entries allow
#'   interval-specific probabilities, matching native GDR removal designs.
#' @param dist_breaks Increasing numeric distance boundaries starting at zero,
#'   in `units_in`. The last boundary is the point-survey radius.
#' @param interval_lengths Positive integer durations of the removal intervals,
#'   measured in `interval_units`. Unequal durations are supported.
#' @param interval_units One nonempty time-unit label, for example `"minutes"`.
#' @param units_in Distance units, `"m"` or `"km"`.
#' @param site_covs Optional data frame of actual or prepared site covariates,
#'   with one row per density. Factors, transformations supplied as columns,
#'   and year assignments are retained in the native frame.
#'   Native character-to-factor conversion and its warning are preserved.
#' @param event_keys Optional data frame with one row per site, retaining the
#'   original event identifiers separately from covariates. Supplying keys does
#'   not aggregate, reorder, or validate their scientific interpretation; use
#'   [validate_unmarked_data()] to audit them.
#' @param seed Nonnegative integer simulation seed. The caller's random-number
#'   state is restored, including its absence when no state existed initially.
#' @param return_joint Logical; return the site-by-distance-by-removal count
#'   array as `joint`. The default avoids retaining this intermediate array.
#' @param provenance Named list retaining caller-supplied generating formulas,
#'   coefficients, kernel definitions, fixed reference transformations, raster
#'   years, or other design information. Values are recorded without evaluation.
#' @return A list with class `msr_gdr_simulation` containing `observed`, an
#'   `unmarkedFrameGDR`; `truth`, a data frame
#'   with site row, expected density, expected circle abundance, realized latent
#'   abundance, distance and removal detection probabilities, total observation
#'   probability, expected observed count, and realized observed count;
#'   `probabilities`, distance-bin and first-removal-interval probability
#'   matrices; `joint`, the optional joint count array (otherwise `NULL`);
#'   `event_keys`; and `design`. The latter records units, survey area in square
#'   kilometers, distance/removal parameters, breaks/durations, seed/RNG kind,
#'   assumptions, native package version, and caller provenance. Truth and
#'   provenance are not inserted into the observed frame's covariates.
#'   Printing the result shows a compact survey and unit summary and returns
#'   the object invisibly. It does not call the native frame printer.
#' @details
#' Survey area is `pi * radius_km^2`. Expected abundance is density times this
#'   area. A bird's distance and removal categories are independent conditional
#'   on the supplied site parameters. First-removal probabilities account for
#'   survival through all preceding intervals; durations convert per-unit
#'   probability `p` to `1 - (1 - p)^duration`.
#'
#' Independent Poisson counts are drawn for every observed joint category and
#'   the unobserved category. Their sum is Poisson latent abundance, and their
#'   conditional allocation given that sum is multinomial. This avoids a loop
#'   over birds or repeated full multinomial draws. Both margins are derived
#'   from these same joint counts.
#'
#' This baseline assumes a closed population within the survey radius and one
#'   primary occasion. It does not simulate movement across the radius,
#'   negative-binomial or zero-inflated abundance, annual availability, or random
#'   detection heterogeneity. Those require separate models and validation.
#'   This helper does not change [sim_dat_unmarked()] or fit any model.
#' @examples
#' sim <- sim_dat_gdr(density = c(100, 200, 0), distance_scale = 40,
#'   removal = 0.3, dist_breaks = c(0, 25, 75),
#'   interval_lengths = c(1, 2, 1), interval_units = "minutes", seed = 24)
#' sim$truth
#' sim
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
#' @param x An object returned by `sim_dat_gdr()`.
#' @param ... Reserved for compatibility with the generic print function.
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
