#' Profile and Bootstrap Inference for an NB2 Scale
#'
#' Examine the full allowed scale range of a single-kernel NB2 fit, retain
#' separate likelihood regions, and optionally calibrate likelihood-ratio
#' tests by parametric bootstrap. This experimental diagnostic does not change
#' the fitted object or the default optimizer and uncertainty methods.
#'
#' @param object A \code{multiScaleR} fit from \code{backend = "profile_glm"}
#'   or \code{backend = "joint_nb2"},
#'   with NB2 errors and exactly one kernel covariate.
#' @param sigma_values Optional positive scale values in map units to examine.
#'   The search bounds and the refined optimum are always included.
#' @param level Confidence level between zero and one; default \code{0.95}.
#' @param grid_points Number of log-spaced scale values when
#'   \code{sigma_values} is omitted; default \code{41}. Increase this to check
#'   whether narrow or disconnected regions were missed.
#' @param bootstrap Number of parametric bootstrap samples per calibrated
#'   scale; zero (the default) skips calibration. Small values have substantial
#'   Monte Carlo uncertainty. At least 999 is advisable for final tail
#'   estimation, although this can be expensive.
#' @param bootstrap_values Optional subset of scales to calibrate, in map
#'   units. They are added to the profile grid. Omitting this argument
#'   calibrates every grid point. A subset supplies pointwise tests, not a
#'   complete bootstrap confidence set.
#' @param seed Nonnegative integer random seed, required when
#'   \code{bootstrap > 0}. The caller's random-number state is restored.
#' @param n_starts Number of additional dispersed starts for inference refits;
#'   default \code{12}. Original starts, the selected fit, and both scale
#'   bounds are also included. Failed starts cannot win.
#' @return A list of class \code{msr_scale_inference} with:
#' \describe{
#'   \item{estimate}{Refined scale estimate in map units.}
#'   \item{reference}{Refined log-parameter vector, negative log likelihood,
#'     start audit, original negative log likelihood, and improvement. The
#'     original fitted object remains unchanged.}
#'   \item{grid}{Scale values, profile negative log likelihood, likelihood-ratio
#'     statistic, chi-squared acceptance, bootstrap p-value and 95 percent
#'     Monte Carlo limits, bootstrap acceptance, and whether the Monte Carlo
#'     limits cross the rejection threshold. Uncalibrated entries are missing.}
#'   \item{profile_set}{Separate accepted scale regions. Columns \code{lower}
#'     and \code{upper} are root-refined endpoints in map units;
#'     \code{lower_outside} and \code{upper_outside} retain adjacent rejected
#'     grid points. The two \code{*_bound_limited} columns mark endpoints
#'     limited by the optimization bounds.}
#'   \item{bootstrap_set}{Accepted grid components with the same columns, or
#'     \code{NULL} unless every grid point was successfully calibrated.
#'     Endpoints are accepted grid values, not root-refined limits; outside
#'     brackets show their resolution.}
#'   \item{diagnostics}{Zero-effect likelihood deviance and compatibility with
#'     the profile cutoff; profile status (\code{localized}, \code{disconnected},
#'     \code{bound_limited}, or \code{unresolved}); original Wald availability;
#'     and bootstrap status, failure count, and Monte Carlo ambiguity count.}
#'   \item{bootstrap_records}{One row per attempted bootstrap refit pair:
#'     generating scale, replicate, success, likelihood-ratio statistic,
#'     fitted scale and theta, boundary indicator, both likelihoods,
#'     attempted/eligible start counts, errors, and warnings. No failed
#'     replicate is silently dropped.}
#'   \item{profile_audit}{Start records at each observed profile grid point.}
#'   \item{level, bootstrap, seed, bounds, call}{Requested confidence level,
#'     bootstrap count and seed, map-unit search bounds, and matched call.}
#' }
#' @details
#' Coefficients and NB2 theta are reestimated at every constrained scale.
#' Exponential-power shape is also reestimated. A stronger unrestricted refit
#' and a full-range grid can reveal an optimum missed by the original fit;
#' the returned reference records this without modifying \code{object}.
#' If conditional NB2 GLM iterations fail for a full-rank design, the backend
#' refines coefficients using the likelihood score and observed information,
#' checks convergence, and rebuilds the conditional GLM result.
#' For \code{joint_nb2} fits, these refits use the boundary-aware conditional
#' likelihood in \code{rho = log(1 + 1/theta)} coordinates, including exact
#' Poisson fits. Bootstrap samples use Poisson sampling when the constrained
#' dispersion estimate is on that boundary. This preserves the boundary in
#' inference without substituting a large finite theta. These refits currently
#' profile coefficients rather than using the joint fitting optimizer.
#'
#' The profile cutoff is \code{qchisq(level, 1)} for twice the log-likelihood
#' difference. Each detected accepted component is retained. Finite grids can
#' still miss narrow regions. A bound-limited endpoint is not a resolved
#' confidence limit.
#'
#' The zero-effect diagnostic fits the model without the kernel coefficient.
#' If that model is inside the profile likelihood cutoff, every admissible
#' scale has a feasible fit inside the cutoff. Scale is then reported as
#' unresolved. This is a likelihood diagnostic, not a calibrated test for
#' the presence of a kernel effect.
#'
#' Bootstrap samples are generated at each constrained scale using its fitted
#' coefficients and theta, with the landscape, fixed design, and offsets held
#' fixed. Every sample is refitted both freely and at the generating scale.
#' The tail estimate is \code{(1 + exceedances) / (bootstrap + 1)}. Monte Carlo
#' limits transform exact binomial limits by the same formula. A test accepts
#' when its tail estimate exceeds \code{1 - level}; the ambiguity flag is
#' separate. Any failed bootstrap pair makes that scale's calibration missing.
#' Only unit observation weights are currently supported for bootstrap.
#'
#' These are plug-in parametric bootstrap tests for independent observations,
#' not exact tests with known nuisance parameters. Neither bootstrap nor
#' positive curvature guarantees calibrated coverage under weak identification
#' or at parameter boundaries. Full-range, disconnected, bound-limited, or
#' Monte Carlo-ambiguous results should not be reduced to a precise scale.
#' Defaults of \code{multiScale_optim()}, \code{summary()}, and \code{predict()}
#' are unchanged.
#' @seealso \code{\link{multiScale_optim}}, \code{\link{profile_sigma}},
#'   \code{\link{diagnostics}}
#' @examples
#' \donttest{
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'   set.seed(42)
#'   habitat <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 600,
#'                          ymin = 0, ymax = 600, crs = "EPSG:32617")
#'   cells <- seq_len(terra::ncell(habitat))
#'   terra::values(habitat) <- sin(cells / 41) + cos(cells / 73)
#'   names(habitat) <- "hab"
#'   sites <- terra::vect(expand.grid(x = seq(150, 450, length.out = 6),
#'                                   y = seq(150, 450, length.out = 6)),
#'                        geom = c("x", "y"), crs = terra::crs(habitat))
#'   ki <- kernel_prep(sites, habitat, max_D = 120, bin = TRUE, verbose = FALSE)
#'   dat <- data.frame(ki$kernel_dat)
#'   dat$y <- rnbinom(nrow(dat), mu = exp(.5 + dat$hab), size = 2)
#'   initial <- glm(y ~ hab, data = dat, family = MASS::negative.binomial(2))
#'   fit <- multiScale_optim(initial, ki, backend = "profile_glm",
#'                            n_cores = 1, verbose = FALSE)
#'   assessment <- scale_inference(fit, grid_points = 9, n_starts = 2)
#'   assessment$profile_set   # Separate regions, with explicit search bounds.
#'   assessment$diagnostics  # A broad region may leave scale unresolved.
#' }
#' }
#' @export
scale_inference <- function(object, sigma_values = NULL, level = .95,
                            grid_points = 41L, bootstrap = 0L,
                            bootstrap_values = NULL, seed = NULL,
                            n_starts = 12L) {
  validate_scalar_numeric(level, "level", lower = 0, upper = 1,
                          inclusive_lower = FALSE, inclusive_upper = FALSE)
  validate_scalar_numeric(grid_points, "grid_points", lower = 3, integerish = TRUE)
  validate_scalar_numeric(bootstrap, "bootstrap", lower = 0, integerish = TRUE)
  setup <- .msr_inference_setup(object, n_starts)
  bounds <- exp(c(setup$lower[1L], setup$upper[1L])) * setup$unit
  check_values <- function(x, label) {
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
        any(x < bounds[1L] * (1 - 1e-12) | x > bounds[2L] * (1 + 1e-12)))
      stop(label, " must contain finite scales within the fitting bounds.", call. = FALSE)
    pmax(bounds[1L], pmin(bounds[2L], x))
  }
  if (is.null(sigma_values)) sigma_values <- exp(seq(log(bounds[1L]), log(bounds[2L]),
                                                    length.out = grid_points))
  sigma_values <- check_values(sigma_values, "sigma_values")
  if (!is.null(bootstrap_values))
    bootstrap_values <- check_values(bootstrap_values, "bootstrap_values")
  if (bootstrap > 0L) {
    if (any(setup$obj$weights != 1))
      stop("Parametric bootstrap currently requires unit observation weights.", call. = FALSE)
    if (1 / (bootstrap + 1) > 1 - level + 1e-12)
      stop("bootstrap is too small to resolve the requested rejection tail.", call. = FALSE)
    validate_scalar_numeric(seed, "seed", lower = 0, upper = .Machine$integer.max,
                            integerish = TRUE)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  } else if (!is.null(bootstrap_values)) {
    stop("bootstrap_values requires bootstrap > 0.", call. = FALSE)
  }
  reference <- .msr_inference_refit(setup)
  scales <- sort(unique(c(bounds, sigma_values, bootstrap_values,
                           exp(reference$par[1L]) * setup$unit)))
  cache <- new.env(parent = emptyenv())
  profile <- function(sigma) {
    key <- sprintf("%.17g", sigma)
    if (!exists(key, cache, inherits = FALSE)) {
      z <- .msr_inference_refit(setup, starts = rbind(reference$par, setup$starts),
                                 fixed_scale = log(sigma / setup$unit))
      assign(key, z, cache)
    }
    get(key, cache, inherits = FALSE)
  }
  profiles <- lapply(scales, profile)
  values <- vapply(profiles, function(z) z$value, numeric(1))
  if (min(values) < reference$value - 1e-5) {
    reference <- .msr_inference_refit(setup,
      starts = rbind(profiles[[which.min(values)]]$par, reference$par, setup$starts))
  }
  if (reference$value > min(values) + 1e-5)
    stop("The profile grid found a better optimum that could not be recovered.", call. = FALSE)
  # Include any newly refined optimum in the grid before constructing regions.
  scales <- sort(unique(c(scales, exp(reference$par[1L]) * setup$unit)))
  profiles <- lapply(scales, profile)
  values <- vapply(profiles, function(z) z$value, numeric(1))
  lr <- 2 * (values - reference$value)
  if (any(lr < -2e-5)) stop("Profile refinement found a better optimum; increase n_starts.", call. = FALSE)
  lr <- pmax(0, lr)
  cutoff <- stats::qchisq(level, 1)
  refine <- function(a, b) {
    fn <- function(z) {
      value <- profile(exp(z))$value
      if (value < reference$value - 1e-5)
        stop("Endpoint refinement found a better optimum; increase n_starts.", call. = FALSE)
      2 * (value - reference$value) - cutoff
    }
    exp(stats::uniroot(fn, log(c(a, b)), tol = 1e-6)$root)
  }
  profile_set <- .msr_scale_components(scales, lr <= cutoff, refine)
  zero_nll <- .msr_nb2_zero_effect(setup)
  zero_lr <- 2 * (zero_nll - reference$value)
  if (zero_lr < -2e-5)
    stop("The zero-effect model is better than the unrestricted fit; increase n_starts.", call. = FALSE)
  zero_lr <- max(0, zero_lr)
  unresolved <- zero_lr <= cutoff || all(lr <= cutoff)
  status <- if (unresolved) "unresolved" else
    if (nrow(profile_set) > 1L) "disconnected" else
      if (any(profile_set$lower_bound_limited | profile_set$upper_bound_limited))
        "bound_limited" else "localized"
  grid <- data.frame(sigma = scales, nll = values, lr = lr,
    profile_accept = lr <= cutoff, bootstrap_p = NA_real_,
    mc_lower = NA_real_, mc_upper = NA_real_, bootstrap_accept = NA,
    mc_ambiguous = NA)
  records <- list()
  if (bootstrap > 0L) {
    selected <- if (is.null(bootstrap_values)) seq_along(scales) else
      which(scales %in% bootstrap_values)
    for (i in selected) {
      z <- .msr_scale_bootstrap(setup, profiles[[i]], reference, bootstrap)
      z$generating_sigma <- scales[i]
      records[[length(records) + 1L]] <- z
      tail <- .msr_bootstrap_tail(z$lr, lr[i], level)
      grid[i, c("bootstrap_p", "mc_lower", "mc_upper")] <- tail
      grid$bootstrap_accept[i] <- tail["p"] > 1 - level
      grid$mc_ambiguous[i] <- tail["mc_lower"] <= 1 - level &&
        tail["mc_upper"] > 1 - level
    }
  }
  records <- if (length(records)) do.call(rbind, records) else data.frame()
  bootstrap_set <- if (bootstrap > 0L)
    .msr_scale_components(scales, grid$bootstrap_accept) else NULL
  failures <- if (nrow(records)) sum(!records$ok) else 0L
  ambiguous <- sum(grid$mc_ambiguous, na.rm = TRUE)
  bootstrap_status <- if (!bootstrap) "not_requested" else if (failures)
    "incomplete" else if (is.null(bootstrap_set)) "pointwise_only" else
      if (ambiguous) "monte_carlo_uncertain" else "grid_approximation"
  reference$original_value <- object$optim_results$value
  reference$improvement <- object$optim_results$value - reference$value
  structure(list(estimate = exp(reference$par[1L]) * setup$unit,
    reference = reference, grid = grid, profile_set = profile_set,
    bootstrap_set = bootstrap_set,
    diagnostics = list(zero_effect_deviance = zero_lr,
      zero_effect_compatible = zero_lr <= cutoff, profile_status = status,
      original_wald_available = object$uncertainty$available,
      bootstrap_status = bootstrap_status, bootstrap_failures = failures,
      monte_carlo_ambiguous = ambiguous),
    bootstrap_records = records,
    profile_audit = lapply(profiles, function(z) z$records),
    level = level, bootstrap = bootstrap, seed = seed, bounds = bounds,
    call = match.call()), class = "msr_scale_inference")
}
