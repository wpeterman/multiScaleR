# Joint estimation of kernel scales across independent data sets (experimental).
#
# Every data set keeps its own fitted model, prepared kernel inputs, and native
# parameters. The joint log-likelihood is the sum of the data-set
# log-likelihoods, each profiled over its native parameters by the same refit
# that multiScale_optim() uses (kernel_scale_fn()). The kernel scales enter
# through one outer parameter vector whose structure ties them across data
# sets: one common scale, one scale per data set, or one ecological scale seen
# through each data set's sampling footprint. Optionally the coefficient of a
# landscape covariate is also shared, entered as a fixed offset in each native
# model so that the native fitter estimates everything else.

#' Estimate scales of effect jointly across independent data sets
#'
#' @description
#' **Experimental.** `multiScale_joint()` combines two or more data sets, each
#' already fitted with [multiScale_optim()], and estimates the scale of effect
#' of each landscape covariate jointly. Each data set keeps its own response,
#' model class, and native parameters (for example a Poisson count model and an
#' `unmarked` occupancy model). What is shared is set by `scales` and, if you
#' choose, `share`:
#'
#' * `scales = "common"`: one scale for each landscape covariate, used by every
#'   data set that contains it.
#' * `scales = "separate"`: one scale per data set. This reproduces the separate
#'   fits and is the reference for testing whether a common scale is supported.
#' * `scales = "footprint"`: one ecological scale, seen through each data set's
#'   sampling footprint (see Details). The apparent scale in data set `d` is
#'   `sqrt(sigma^2 + footprint[d]^2)`.
#' * `scales = "footprint_scaled"`: as `"footprint"`, with the footprint SDs
#'   multiplied by an estimated factor `kappa`. `kappa = 0` is the common
#'   scale and `kappa = 1` is the stated footprint.
#'
#' `share` names landscape covariates whose regression coefficient is also
#' shared across data sets, so the data sets estimate one landscape effect but
#' keep their own intercepts and other coefficients.
#'
#' @param fits Named list of two or more fitted `multiScaleR` objects from
#'   [multiScale_optim()] with the default `backend = "refit"`, one per data
#'   set. The names identify the data sets in the output. Each fit supplies
#'   its starting model, its prepared kernel inputs, and starting values. All
#'   fits must use the same kernel, Gaussian or exponential. Landscape
#'   covariates are matched across data sets by name, so give the same raster
#'   layer the same name in every [kernel_prep()] call.
#' @param scales How the scale of each landscape covariate relates across data
#'   sets: one of `"common"` (default), `"separate"`, `"footprint"`, or
#'   `"footprint_scaled"` for every covariate, or a named character vector with
#'   one structure per covariate, such as `c(forest = "common", road =
#'   "separate")`. A covariate that occurs in only one data set always has its
#'   own scale.
#' @param footprint Named numeric vector with one sampling footprint SD per
#'   data set, in map units, required by the footprint structures. Use 0 for a
#'   data set whose response is measured at a point. [footprint_sd()] converts
#'   common sampling geometries to this SD.
#' @param share Optional character vector of landscape covariates whose
#'   coefficient is shared across the data sets that contain them. Each such
#'   covariate must enter every model as one plain linear term, and the models
#'   must use comparable links (all the same link, or a mix of log and
#'   complementary log-log, which share the log-intensity scale). Shared
#'   coefficients currently require models refitted by `update()`: `glm`,
#'   `MASS::glm.nb`, lme4, glmmTMB, and mgcv. Default `NULL` shares no
#'   coefficient.
#' @param cluster Optional site or block identifiers, used to check for and
#'   allow dependence between data sets that sample the same places. Either
#'   the name of a column present in every model's data, or a named list with
#'   one vector per data set aligned with that model's data rows (or `unmarked`
#'   sites). Default `NULL`.
#' @param uncertainty `"hessian"` (default) takes standard errors from the
#'   Hessian of the joint log-likelihood, which assumes the data sets are
#'   independent. `"jackknife"` deletes one group of clusters at a time and
#'   refits, which allows dependence within clusters, including between data
#'   sets that share a site. Its Wald limits use a t distribution with one
#'   fewer degree of freedom than the number of replicates. It requires
#'   `cluster` and is not yet available for `unmarked` models.
#' @param profile Logical. If `TRUE` (default), compute profile-likelihood
#'   intervals for every estimated parameter.
#' @param level Confidence level for the intervals. Default `0.95`.
#' @param groups Number of jackknife groups when `uncertainty = "jackknife"`.
#'   Clusters are assigned to groups at random. Default `20`, or the number of
#'   clusters if smaller.
#' @param seed Optional integer seed for the jackknife group assignment.
#' @param verbose Logical. If `TRUE` (default), report progress.
#'
#' @details
#' **The joint likelihood.** For a given set of scales, each data set's model
#' is refitted exactly as [multiScale_optim()] refits it, and its maximized
#' log-likelihood is returned. `multiScale_joint()` maximizes the sum of these
#' log-likelihoods over the scale parameters (and shared coefficients). With
#' `scales = "separate"` and no shared coefficient the sum separates, and the
#' result equals the separate fits. Because the joint fit treats the data sets
#' as independent pieces of evidence, it is valid when they are independent
#' given the covariates. Data sets that sample the same sites may share
#' unexplained site effects; supply `cluster` to measure the residual
#' correlation at shared sites, and use `uncertainty = "jackknife"` when it is
#' not negligible.
#'
#' **Comparing structures.** The returned object reports the likelihood-ratio
#' test of the fitted structure against separate scales. With two or three
#' data sets this test has little power, so a nonsignificant result is weak
#' evidence for a common scale. Fit several structures and compare them with
#' [stats::AIC()] (the fits have [stats::logLik()] methods) or [stats::anova()].
#' The test of `kappa = 0` in `"footprint_scaled"` lies on the boundary of the
#' parameter space, so its p-value is halved.
#'
#' **Sampling footprints.** A response is often measured over an area, not at a
#' point: counts within a fixed radius, detections by an acoustic recorder, a
#' camera viewshed, or animals detected by distance sampling. If the
#' ecological response follows the landscape smoothed with a Gaussian kernel of
#' scale `sigma`, then a response averaged over a footprint with per-axis SD
#' `tau` follows the landscape smoothed by the convolution of the two, whose
#' per-axis variance is `sigma^2 + tau^2`. Data sets with larger footprints
#' therefore show larger apparent scales even when the ecological scale is the
#' same. The relation is exact for a Gaussian footprint and for a covariate
#' that enters the linear predictor linearly; for a disk or square footprint it
#' matches the second moment, and on a log or logit link it is a first-order
#' approximation that holds while the landscape effect varies little within
#' one footprint. Footprint structures require the Gaussian kernel.
#'
#' **Shared coefficients.** For each shared covariate, every data set's model
#' is refitted with that term replaced by a fixed offset equal to the shared
#' coefficient times the covariate. The covariate is standardized with one
#' center and spread pooled across the data sets, evaluated at the ecological
#' scale, so the coefficient has the same meaning in every data set.
#'
#' **Component fits.** `$fits` holds one `multiScaleR` object per data set at
#' the jointly estimated scales, with standard errors from the joint fit. Use
#' them with [plot.multiScaleR()], [plot_marginal_effects()],
#' [kernel_scale.raster()], and `predict()` as for any other fit. Two cautions
#' apply: `summary(fit, profile = TRUE)` and [profile_sigma()] on a component
#' profile that data set alone, and a shared coefficient enters a component's
#' model as an offset named `.msr_joint_offset`, so effect curves and maps of
#' that covariate are not yet available from the component.
#'
#' @return An object of class `multiScaleR_joint`, a list with
#' * `parameters`: estimated parameters (scales in map units, `kappa`, shared
#'   coefficients) with standard errors and Wald and profile-likelihood limits.
#' * `scales`: the apparent scale of each covariate in each data set.
#' * `fits`: named list of component `multiScaleR` fits.
#' * `comparison`: log-likelihoods, parameter counts, and the likelihood-ratio
#'   test against separate scales.
#' * `residual_correlation`: correlation of residuals between data sets at
#'   shared clusters, when `cluster` is supplied.
#' * `logLik`, `df`, `nobs`, `vcov` (parameters on the optimizer scale: log
#'   scales, `kappa`, coefficients), `structure`, `footprint`, `share`,
#'   `standardization`, `uncertainty`, `optim`, and `call`.
#'
#' @seealso [multiScale_optim()], [footprint_sd()], [profile_sigma()]
#'
#' @examples
#' \donttest{
#' # Two independent surveys of one landscape: Poisson counts and binary
#' # detections. Both respond to `bin1` and `cont2` at the same true scales.
#' rs <- sim_rast(user_seed = 123, dim = 100, resolution = 30)
#' rs <- terra::subset(rs, c("bin1", "cont2"))
#'
#' counts <- sim_dat(alpha = 0.5, beta = c(0.75, -0.75), kernel = "gaussian",
#'                   sigma = c(75, 150), type = "count", n_points = 60,
#'                   raster_stack = rs, max_D = 450, user_seed = 1)
#' occ <- sim_dat(alpha = 0, beta = c(1, -1), kernel = "gaussian",
#'                sigma = c(75, 150), type = "occ", n_points = 80,
#'                raster_stack = rs, max_D = 450, user_seed = 2)
#'
#' in_counts <- kernel_prep(counts$pts, rs, max_D = 450, verbose = FALSE)
#' in_occ <- kernel_prep(occ$pts, rs, max_D = 450, verbose = FALSE)
#'
#' fit_counts <- multiScale_optim(
#'   glm(y ~ bin1 + cont2, family = poisson,
#'       data = data.frame(y = counts$obs, in_counts$kernel_dat)),
#'   in_counts, verbose = FALSE)
#' fit_occ <- multiScale_optim(
#'   glm(y ~ bin1 + cont2, family = binomial,
#'       data = data.frame(y = occ$obs, in_occ$kernel_dat)),
#'   in_occ, verbose = FALSE)
#'
#' joint <- multiScale_joint(list(counts = fit_counts, occ = fit_occ),
#'                           scales = "common", verbose = FALSE)
#' joint
#'
#' # Compare with one scale per data set
#' separate <- multiScale_joint(list(counts = fit_counts, occ = fit_occ),
#'                              scales = "separate", profile = FALSE,
#'                              verbose = FALSE)
#' AIC(joint, separate)
#' }
#' @export
multiScale_joint <- function(fits,
                             scales = "common",
                             footprint = NULL,
                             share = NULL,
                             cluster = NULL,
                             uncertainty = c("hessian", "jackknife"),
                             profile = TRUE,
                             level = 0.95,
                             groups = 20L,
                             seed = NULL,
                             verbose = TRUE) {
  uncertainty <- match.arg(uncertainty)
  validate_scalar_logical(profile, "profile")
  validate_scalar_logical(verbose, "verbose")
  validate_scalar_numeric(level, "level", lower = 0, upper = 1,
                          inclusive_lower = FALSE, inclusive_upper = FALSE)

  sets <- .msr_joint_sets(fits)
  structure <- .msr_joint_structure(scales, sets)
  tau <- .msr_joint_footprint(footprint, structure, sets)
  shared <- .msr_joint_shared(share, structure, sets)
  clusters <- .msr_joint_clusters(cluster, sets)
  if (identical(uncertainty, "jackknife")) {
    if (is.null(clusters))
      stop("uncertainty = 'jackknife' requires `cluster`.", call. = FALSE)
    validate_scalar_numeric(groups, "groups", integerish = TRUE, lower = 2)
  }
  spec <- .msr_joint_spec(sets, structure, tau, shared)
  fn <- .msr_joint_objective(spec)

  if (verbose) message("Fitting ", nrow(spec$parameters), " joint parameter(s) across ",
                       length(sets), " data sets...")
  opt <- .msr_joint_optimize(spec, fn)
  theta <- stats::setNames(opt$par, spec$parameters$parameter)

  hess <- .msr_joint_hessian(fn, theta, spec$parameters)
  vcov_theta <- hess$vcov
  jack <- NULL
  if (identical(uncertainty, "jackknife")) {
    if (verbose) message("Jackknife over cluster groups...")
    jack <- .msr_joint_jackknife(spec, theta, clusters, groups, seed)
    vcov_theta <- jack$vcov
  }
  dimnames(vcov_theta) <- list(names(theta), names(theta))

  prof <- NULL
  if (profile) {
    if (verbose) message("Profile-likelihood intervals...")
    prof <- .msr_joint_profile(fn, theta, opt$value, spec$parameters, level)
  }

  # A jackknife variance from m replicates has m - 1 degrees of freedom.
  quantile <- if (is.null(jack)) stats::qnorm(1 - (1 - level) / 2) else
    stats::qt(1 - (1 - level) / 2, jack$replicates - 1L)
  parameters <- .msr_joint_parameter_table(spec, theta, vcov_theta, prof, quantile, hess$at_bound)
  sigma <- .msr_joint_sigma(theta, spec)
  scale_table <- .msr_joint_scale_table(spec, theta, vcov_theta, prof, quantile)
  beta <- stats::setNames(theta[spec$beta_index], spec$shared)
  std <- if (length(shared)) .msr_joint_standardize(theta, sigma, spec) else NULL

  info <- list(structure = structure, footprint = tau, share = shared,
               data_sets = names(sets), uncertainty = uncertainty)
  components <- lapply(spec$sets, function(s) {
    rows <- scale_table$data_set == s$name
    se <- stats::setNames(scale_table$SE[rows], scale_table$covariate[rows])
    .msr_joint_component(s, sigma[[s$name]], se, beta, std,
                         utils::modifyList(info, list(data_set = s$name)))
  })
  names(components) <- names(sets)

  native_k <- vapply(components, function(x) .msr_parameter_count(x$opt_mod), numeric(1))
  ll <- -opt$value
  k <- sum(native_k) + nrow(spec$parameters)
  n <- sum(vapply(sets, `[[`, numeric(1), "nobs"))
  reference <- .msr_joint_separate_reference(spec, sigma)
  comparison <- .msr_joint_comparison(spec, ll, k, reference)

  residual_correlation <- if (is.null(clusters)) NULL else
    .msr_joint_residual_correlation(components, clusters)

  out <- list(parameters = parameters,
              scales = scale_table,
              fits = components,
              comparison = comparison,
              residual_correlation = residual_correlation,
              logLik = ll, df = k, nobs = n,
              vcov = vcov_theta,
              structure = .msr_joint_structure_table(spec),
              footprint = tau,
              share = shared,
              standardization = std,
              uncertainty = list(method = uncertainty,
                                 at_bound = names(theta)[hess$at_bound],
                                 hessian_vcov = hess$vcov,
                                 jackknife = jack),
              level = level,
              kernel = spec$kernel,
              optim = opt,
              call = match.call())
  dimnames(out$uncertainty$hessian_vcov) <- list(names(theta), names(theta))
  class(out) <- "multiScaleR_joint"
  improved <- comparison$improved
  if (!is.null(improved) && nrow(improved))
    warning("The separate fit of ", paste(sprintf("`%s` (logLik %.3f to %.3f at sigma %s)",
              improved$data_set, improved$logLik_input, improved$logLik_refit,
              improved$sigma), collapse = "; "),
            " was not at its maximum. The comparison with separate scales uses the improved ",
            "values; refit those data sets before reporting them alone.", call. = FALSE)
  out
}


#' Convert a sampling geometry to a footprint SD
#'
#' @description
#' **Experimental.** Returns the per-axis standard deviation, in map units, of
#' the area over which a response is measured. [multiScale_joint()] uses it to
#' relate the apparent scale of effect in a data set to the ecological scale
#' (`scales = "footprint"`).
#'
#' @param type Sampling geometry: `"point"` (no footprint), `"disk"` (counts or
#'   detections within a fixed radius), `"square"` (a square plot or grid
#'   cell), `"rectangle"` (a plot or transect strip), or `"halfnormal"` (animals
#'   detected with a half-normal detection function, as in distance sampling).
#' @param size Radius for `"disk"`, side length for `"square"`, length for
#'   `"rectangle"`, and the detection-function scale for `"halfnormal"`, in map
#'   units. Ignored for `"point"`.
#' @param width Width for `"rectangle"`, in map units.
#'
#' @details
#' The SD is the square root of the per-axis variance of a point drawn at
#' random from the footprint, weighted by detection where that applies: `r / 2`
#' for a disk of radius `r`, `L / sqrt(12)` for a square of side `L`,
#' `sqrt((L^2 + W^2) / 24)` for an `L` by `W` rectangle (its two axes averaged,
#' since the Gaussian kernel is isotropic), and the detection scale itself for a
#' half-normal detection function, whose detection weights over the plane form
#' a bivariate normal distribution.
#'
#' @return A numeric vector of footprint SDs in map units.
#' @seealso [multiScale_joint()]
#' @examples
#' footprint_sd("disk", 100)          # 100 m fixed-radius point count: 50
#' footprint_sd("halfnormal", 60)     # distance sampling, sigma = 60 m: 60
#' footprint_sd("rectangle", 200, 20) # 200 m by 20 m transect strip
#' @export
footprint_sd <- function(type = c("point", "disk", "square", "rectangle", "halfnormal"),
                         size = NULL, width = NULL) {
  type <- match.arg(type)
  if (type == "point") return(0)
  if (!is.numeric(size) || !length(size) || any(!is.finite(size) | size <= 0))
    stop("`size` must contain positive distances in map units.", call. = FALSE)
  if (type == "rectangle" &&
      (!is.numeric(width) || !length(width) || any(!is.finite(width) | width <= 0)))
    stop("`width` must contain positive distances in map units for a rectangle.",
         call. = FALSE)
  switch(type,
         disk = size / 2,
         square = size / sqrt(12),
         rectangle = sqrt((size^2 + width^2) / 24),
         halfnormal = size)
}


# ---- Inputs ------------------------------------------------------------------

# One record per data set: the separate fit and everything its objective needs.
.msr_joint_sets <- function(fits) {
  if (!is.list(fits) || inherits(fits, "multiScaleR") || length(fits) < 2L)
    stop("`fits` must be a list of two or more multiScaleR fits, one per data set.",
         call. = FALSE)
  nm <- names(fits)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm))
    stop("`fits` must be a named list with one unique name per data set.", call. = FALSE)
  sets <- lapply(nm, function(d) {
    x <- fits[[d]]
    if (!inherits(x, "multiScaleR"))
      stop("Every element of `fits` must be a fitted multiScaleR object from ",
           "multiScale_optim(); `", d, "` is not.", call. = FALSE)
    if (!is.null(x$joint))
      stop("`", d, "` is a component of a joint fit. Supply the separate fits from ",
           "multiScale_optim().", call. = FALSE)
    if (!is.null(x$backend) && !identical(x$backend, "refit"))
      stop("multiScale_joint() currently combines fits from the default refit backend; `",
           d, "` used backend = '", x$backend, "'.", call. = FALSE)
    if (!is.null(x$diagnostics$radius_profile))
      stop("multiScale_joint() cannot yet combine radius_grid fits (`", d, "`).",
           call. = FALSE)
    ki <- x$kernel_inputs
    ctx <- x$opt_context
    if (is.null(ctx) || is.null(x$fitted_mod_original) || is.null(ki$kernel) ||
        (is.null(ki$d_list) && is.null(ki$binned)))
      stop("`", d, "` lacks the stored optimization inputs. Refit it with multiScale_optim().",
           call. = FALSE)
    if (!ki$kernel %in% c("gaussian", "exp"))
      stop("multiScale_joint() currently supports Gaussian and exponential kernels; `", d,
           "` uses kernel = '", ki$kernel, "'.", call. = FALSE)
    covs <- ctx$covs
    hard <- .msr_hard_radius_covariates(ki$scale_vars, ki$kernel, covs)
    if (length(hard))
      stop("Hard-radius metric covariates cannot yet be combined across data sets (`", d,
           "`: ", paste(hard, collapse = ", "), ").", call. = FALSE)
    list(name = d, fit = x, model = x$fitted_mod_original, ki = ki, ctx = ctx,
         covs = covs, join_by = x$join_by, kernel = ki$kernel, unit = ki$unit_conv,
         min_D = x$min_D, max_D = x$max_D,
         sigma_sep = stats::setNames(x$scale_est[covs, "Mean"], covs),
         loglik_sep = .msr_model_loglik(x$opt_mod),
         k_sep = .msr_parameter_count(x$opt_mod) + length(covs),
         nobs = .msr_model_nobs(x$opt_mod),
         model_class = class(.analysis_model(x$opt_mod))[1L])
  })
  names(sets) <- nm
  kernels <- unique(vapply(sets, `[[`, "", "kernel"))
  if (length(kernels) != 1L)
    stop("All data sets must use the same kernel; found ", paste(kernels, collapse = " and "),
         ".", call. = FALSE)
  sets
}


.msr_joint_members <- function(sets) {
  covariates <- unique(unlist(lapply(sets, `[[`, "covs"), use.names = FALSE))
  stats::setNames(lapply(covariates, function(j)
    names(sets)[vapply(sets, function(s) j %in% s$covs, logical(1))]), covariates)
}


.msr_joint_structure <- function(scales, sets) {
  allowed <- c("common", "separate", "footprint", "footprint_scaled")
  members <- .msr_joint_members(sets)
  covariates <- names(members)
  if (!is.character(scales) || !length(scales) || anyNA(scales))
    stop("`scales` must be one of ", paste(sprintf("'%s'", allowed), collapse = ", "),
         ", or a named character vector of these.", call. = FALSE)
  if (is.null(names(scales))) {
    if (length(scales) != 1L)
      stop("`scales` must be one structure, or a named vector with one structure per ",
           "landscape covariate.", call. = FALSE)
    out <- stats::setNames(rep(scales, length(covariates)), covariates)
  } else {
    unknown <- setdiff(names(scales), covariates)
    if (length(unknown))
      stop("`scales` names covariates that no fit uses: ", paste(unknown, collapse = ", "),
           ".", call. = FALSE)
    shared_covs <- covariates[lengths(members) > 1L]
    missing <- setdiff(shared_covs, names(scales))
    if (length(missing))
      stop("`scales` must give a structure for every covariate used by two or more data ",
           "sets; missing: ", paste(missing, collapse = ", "), ".", call. = FALSE)
    out <- stats::setNames(rep("separate", length(covariates)), covariates)
    out[names(scales)] <- scales
  }
  bad <- setdiff(unique(out), allowed)
  if (length(bad))
    stop("Unknown scale structure: ", paste(bad, collapse = ", "), ". Use ",
         paste(sprintf("'%s'", allowed), collapse = ", "), ".", call. = FALSE)
  # A covariate measured in one data set has nothing to share.
  out[lengths(members) == 1L] <- "separate"
  out
}


.msr_joint_footprint <- function(footprint, structure, sets) {
  if (!any(structure %in% c("footprint", "footprint_scaled"))) {
    if (!is.null(footprint))
      warning("`footprint` is ignored because no covariate uses a footprint structure.",
              call. = FALSE)
    return(NULL)
  }
  if (!identical(sets[[1L]]$kernel, "gaussian"))
    stop("Footprint structures require the Gaussian kernel, for which a sampling footprint ",
         "adds to the kernel variance.", call. = FALSE)
  if (!is.numeric(footprint) || is.null(names(footprint)) ||
      !setequal(names(footprint), names(sets)) || anyDuplicated(names(footprint)) ||
      any(!is.finite(footprint) | footprint < 0))
    stop("`footprint` must be a named numeric vector of nonnegative footprint SDs in map ",
         "units, one per data set (", paste(names(sets), collapse = ", "),
         "). See footprint_sd().", call. = FALSE)
  footprint <- footprint[names(sets)]
  members <- .msr_joint_members(sets)
  for (j in names(structure)[structure == "footprint_scaled"]) {
    if (length(unique(footprint[members[[j]]])) < 2L)
      stop("scales = 'footprint_scaled' needs at least two different footprint sizes among ",
           "the data sets that contain `", j, "`.", call. = FALSE)
  }
  footprint
}


.msr_joint_shared <- function(share, structure, sets) {
  if (is.null(share)) return(character())
  members <- .msr_joint_members(sets)
  if (!is.character(share) || !length(share) || anyNA(share) || anyDuplicated(share))
    stop("`share` must be a character vector of landscape covariate names.", call. = FALSE)
  unknown <- setdiff(share, names(members))
  if (length(unknown))
    stop("`share` names covariates that no fit uses: ", paste(unknown, collapse = ", "),
         ".", call. = FALSE)
  single <- share[lengths(members[share]) < 2L]
  if (length(single))
    stop("A shared coefficient needs two or more data sets; only one uses: ",
         paste(single, collapse = ", "), ".", call. = FALSE)
  for (j in share) {
    links <- character()
    for (d in members[[j]]) {
      s <- sets[[d]]
      model <- .analysis_model(s$model)
      if (methods::is(model, "unmarkedFit") || inherits(model, "coxph") ||
          !is.null(s$ctx$refit_fn) || identical(s$ctx$mod_class, "gls"))
        stop("Shared coefficients currently need models refitted by update(), such as ",
             "glm, MASS::glm.nb, lme4, glmmTMB, and mgcv; `", d, "` is a ",
             class(model)[1L], " fit.", call. = FALSE)
      labels <- attr(stats::terms(stats::formula(model)), "term.labels")
      if (!j %in% labels)
        stop("A shared covariate must enter each model as one plain linear term; `", j,
             "` does not in `", d, "`.", call. = FALSE)
      inside <- vapply(setdiff(labels, j), function(t)
        j %in% all.vars(str2lang(t)), logical(1))
      if (any(inside))
        stop("A shared covariate cannot also appear in another term (interaction, ",
             "polynomial, or smooth); `", j, "` does in `", d, "`.", call. = FALSE)
      link <- tryCatch(stats::family(model)$link, error = function(e) NA_character_)
      if (!is.character(link) || length(link) != 1L || is.na(link))
        stop("Could not determine the link function of `", d, "`.", call. = FALSE)
      links[d] <- link
    }
    if (length(unique(links)) > 1L && !all(links %in% c("log", "cloglog")))
      stop("A shared coefficient needs comparable links. `", j, "` enters models with links ",
           paste(sprintf("%s (%s)", names(links), links), collapse = ", "),
           ". Use the same link, or log for counts with cloglog for binary data.",
           call. = FALSE)
  }
  share
}


.msr_joint_clusters <- function(cluster, sets) {
  if (is.null(cluster)) return(NULL)
  ids <- lapply(sets, function(s) {
    values <- if (is.character(cluster) && length(cluster) == 1L) {
      dat <- .msr_joint_model_rows(s)
      if (cluster %in% names(dat)) dat[[cluster]] else {
        # Native model data keep only formula variables; look up the column in
        # the data frame the model was fitted to, matched by row name.
        model <- .analysis_model(s$model)
        full <- tryCatch(eval(stats::getCall(model)$data,
                              envir = environment(stats::formula(model))),
                         error = function(e) NULL)
        if (!is.data.frame(full) || !cluster %in% names(full))
          stop("`cluster` column '", cluster, "' is not in the data of `", s$name, "`.",
               call. = FALSE)
        idx <- match(row.names(dat), row.names(full))
        if (anyNA(idx) && nrow(full) == nrow(dat)) idx <- seq_len(nrow(dat))
        full[[cluster]][idx]
      }
    } else if (is.list(cluster) && !is.null(names(cluster)) &&
               setequal(names(cluster), names(sets))) {
      cluster[[s$name]]
    } else {
      stop("`cluster` must be a column name present in every model's data, or a named list ",
           "with one vector per data set.", call. = FALSE)
    }
    n <- nrow(.msr_joint_model_rows(s))
    if (length(values) != n || anyNA(values))
      stop("Cluster identifiers for `", s$name, "` must have one nonmissing value per ",
           "model data row or site (", n, ").", call. = FALSE)
    as.character(values)
  })
  ids
}


# Rows of a data set's model: native data rows, or unmarked sites.
.msr_joint_model_rows <- function(s) {
  if (identical(s$ctx$mod_class, "unmarked")) s$ctx$umf_template@siteCovs else
    s$ctx$data_template
}


# ---- Parameters --------------------------------------------------------------

.msr_joint_spec <- function(sets, structure, tau, shared) {
  members <- .msr_joint_members(sets)
  rows <- list()
  add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
  for (j in names(structure)) {
    in_sets <- members[[j]]
    lo <- vapply(sets[in_sets], `[[`, numeric(1), "min_D")
    hi <- vapply(sets[in_sets], `[[`, numeric(1), "max_D")
    sep <- vapply(sets[in_sets], function(s) s$sigma_sep[[j]], numeric(1))
    st <- structure[[j]]
    if (st == "separate") {
      for (d in in_sets)
        add(parameter = paste0("sigma:", j, ":", d), type = "log_sigma", covariate = j,
            data_set = d, lower = log(lo[[d]]), upper = log(hi[[d]]), start = log(sep[[d]]))
    } else if (st == "common") {
      if (max(lo) >= min(hi))
        stop("The search ranges of the data sets do not overlap for `", j, "`. Prepare ",
             "their kernel inputs with overlapping min_D and max_D.", call. = FALSE)
      add(parameter = paste0("sigma:", j), type = "log_sigma", covariate = j,
          data_set = NA_character_, lower = log(max(lo)), upper = log(min(hi)),
          start = mean(log(sep)))
    } else {
      t <- tau[in_sets]
      if (any(hi <= t))
        stop("A footprint SD must be smaller than max_D of its data set (`", j, "`).",
             call. = FALSE)
      floor_sigma <- min(lo) / 10
      lower <- max(floor_sigma, sqrt(pmax(lo^2 - t^2, 0)))
      upper <- min(sqrt(hi^2 - t^2))
      if (lower >= upper)
        stop("No ecological scale for `", j, "` keeps every data set's apparent scale ",
             "within its min_D and max_D.", call. = FALSE)
      eco <- sqrt(pmax(sep^2 - t^2, floor_sigma^2))
      add(parameter = paste0("sigma:", j), type = "log_sigma", covariate = j,
          data_set = NA_character_, lower = log(lower), upper = log(upper),
          start = mean(log(eco)))
      if (st == "footprint_scaled")
        add(parameter = paste0("kappa:", j), type = "kappa", covariate = j,
            data_set = NA_character_, lower = 0, upper = 10, start = 1)
    }
  }
  for (j in shared) {
    start <- mean(vapply(sets[members[[j]]], function(s) {
      p <- tryCatch(insight::get_parameters(.analysis_model(s$fit$opt_mod)),
                    error = function(e) NULL)
      v <- if (is.data.frame(p)) p$Estimate[match(j, p$Parameter)] else NA_real_
      if (length(v) == 1L && is.finite(v)) v else 0
    }, numeric(1)))
    add(parameter = paste0("beta:", j), type = "beta", covariate = j,
        data_set = NA_character_, lower = -Inf, upper = Inf, start = start)
  }
  par <- do.call(rbind, rows)
  par$start <- pmin(pmax(par$start, par$lower), par$upper)
  spec <- list(sets = sets, structure = structure, tau = tau, shared = shared,
               members = members, parameters = par, kernel = sets[[1L]]$kernel,
               beta_index = which(par$type == "beta"))
  # Data sets with a shared coefficient refit an offset version of their model.
  for (d in names(sets)) {
    shared_d <- intersect(shared, sets[[d]]$covs)
    if (length(shared_d)) {
      spec$sets[[d]]$shared <- shared_d
      spec$sets[[d]]$offset_model <- .msr_joint_offset_model(sets[[d]], shared_d)
    }
  }
  spec
}


# Apparent scale (map units) of every covariate in every data set.
.msr_joint_sigma <- function(theta, spec, eco = FALSE) {
  par <- spec$parameters
  out <- lapply(spec$sets, function(s) stats::setNames(rep(NA_real_, length(s$covs)), s$covs))
  for (j in names(spec$structure)) {
    st <- spec$structure[[j]]
    rows <- which(par$covariate == j & par$type != "beta")
    if (st == "separate") {
      for (r in rows) out[[par$data_set[r]]][[j]] <- exp(theta[[r]])
      next
    }
    sigma <- exp(theta[[rows[par$type[rows] == "log_sigma"]]])
    k <- if (eco || st == "common") 0 else if (st == "footprint") 1 else
      theta[[rows[par$type[rows] == "kappa"]]]
    for (d in spec$members[[j]])
      out[[d]][[j]] <- if (k == 0) sigma else sqrt(sigma^2 + (k * spec$tau[[d]])^2)
  }
  out
}


# One center and spread per shared covariate, pooled across data sets at the
# ecological scale (or at each data set's own scale under "separate").
.msr_joint_standardize <- function(theta, sigma, spec) {
  ref <- .msr_joint_sigma(theta, spec, eco = TRUE)
  center <- spread <- stats::setNames(numeric(length(spec$shared)), spec$shared)
  for (j in spec$shared) {
    values <- unlist(lapply(spec$members[[j]], function(d) {
      s <- spec$sets[[d]]
      par <- unname(ref[[d]][s$covs]) / s$unit
      .msr_kernel_cov_w(par, s$ki$d_list, s$ki$raw_cov, s$kernel, s$ctx, covariates = j)[, 1]
    }), use.names = FALSE)
    if (anyNA(values) || length(values) < 3L || !is.finite(stats::sd(values)) ||
        stats::sd(values) <= 0) return(NULL)
    center[[j]] <- mean(values)
    spread[[j]] <- stats::sd(values)
  }
  list(center = center, spread = spread)
}


# ---- Objective -----------------------------------------------------------------

.msr_joint_penalty <- 1e60

.msr_joint_objective <- function(spec) {
  function(theta) {
    sigma <- .msr_joint_sigma(theta, spec)
    for (s in spec$sets) {
      v <- sigma[[s$name]]
      if (any(!is.finite(v)) || any(v < s$min_D * (1 - 1e-8)) || any(v > s$max_D * (1 + 1e-8)))
        return(.msr_joint_penalty)
    }
    beta <- stats::setNames(theta[spec$beta_index], spec$shared)
    std <- NULL
    if (length(spec$shared)) {
      std <- tryCatch(.msr_joint_standardize(theta, sigma, spec), error = function(e) NULL)
      if (is.null(std)) return(.msr_joint_penalty)
    }
    total <- 0
    for (s in spec$sets) {
      # Native warnings at trial values (for example glm.fit convergence at an
      # extreme offset) say nothing about the final fit, whose refit keeps them.
      v <- tryCatch(suppressWarnings(.msr_joint_set_value(s, sigma[[s$name]], beta, std)),
                    error = function(e) .msr_joint_penalty)
      if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v >= 1e50)
        return(.msr_joint_penalty)
      total <- total + v
    }
    total
  }
}


# Negative log-likelihood of one data set at its apparent scales, or its final
# refitted model when return_model = TRUE.
.msr_joint_set_value <- function(s, sigma, beta, std, return_model = FALSE) {
  par <- unname(sigma[s$covs]) / s$unit
  if (is.null(s$shared)) {
    return(kernel_scale_fn(par = par, d_list = s$ki$d_list, cov_df = s$ki$raw_cov,
                           kernel = s$kernel, fitted_mod = s$model, join_by = s$join_by,
                           mod_return = if (return_model) TRUE else NULL,
                           opt_context = s$ctx))
  }
  fail <- function(msg) if (return_model) stop(msg, call. = FALSE) else .msr_joint_penalty
  raw <- .msr_kernel_cov_w(par, s$ki$d_list, s$ki$raw_cov, s$kernel, s$ctx)
  if (anyNA(raw) || any(!is.finite(raw)))
    return(fail(paste0("An optimized covariate is missing for an observation in `", s$name, "`.")))
  free <- setdiff(s$covs, s$shared)
  scl <- raw
  center <- spread <- stats::setNames(rep(NA_real_, length(s$covs)), s$covs)
  if (length(free)) {
    z <- scale(raw[, free, drop = FALSE])
    scl[, free] <- z
    center[free] <- attr(z, "scaled:center")
    spread[free] <- attr(z, "scaled:scale")
  }
  for (j in s$shared) {
    scl[, j] <- (raw[, j] - std$center[[j]]) / std$spread[[j]]
    center[[j]] <- std$center[[j]]
    spread[[j]] <- std$spread[[j]]
  }
  dat <- s$ctx$data_template
  scl_dat <- align_scaled_covariates(scl, dat, s$ctx)
  dat[s$covs] <- as.data.frame(scl_dat[, s$covs, drop = FALSE])
  dat$.msr_joint_offset <- drop(scl_dat[, s$shared, drop = FALSE] %*% beta[s$shared])
  context <- list(mod_class = s$ctx$mod_class, return_model = return_model)
  refit_error <- NULL
  mod <- tryCatch(.refit_model(s$offset_model, data = dat, opt_context = context),
                  error = function(e) { refit_error <<- conditionMessage(e); NULL })
  if (is.null(mod))
    return(fail(paste0("Failed to refit `", s$name, "` with the shared coefficient: ",
                       refit_error)))
  if (!return_model) {
    value <- tryCatch(-.msr_model_loglik(mod), error = function(e) NA_real_)
    return(if (is.finite(value)) value else .msr_joint_penalty)
  }
  list(mod = mod, scl_params = list(mean = center, sd = spread))
}


# The native model with each shared covariate replaced by a fixed offset.
.msr_joint_offset_model <- function(s, shared_d) {
  model <- .analysis_model(s$model)
  drop <- paste(shared_d, collapse = " - ")
  new_formula <- stats::as.formula(paste(". ~ . -", drop, "+ offset(.msr_joint_offset)"))
  dat <- s$ctx$data_template
  dat$.msr_joint_offset <- 0
  fitted <- tryCatch(stats::update(model, formula. = new_formula, data = dat),
                     error = function(e) e)
  if (inherits(fitted, "error"))
    stop("Could not refit `", s$name, "` with its shared coefficient as an offset: ",
         conditionMessage(fitted), call. = FALSE)
  fitted
}


# ---- Optimization and uncertainty ---------------------------------------------

.msr_joint_optimize <- function(spec, fn) {
  par <- spec$parameters
  starts <- list(par$start)
  alt <- par$start
  scale_rows <- par$type == "log_sigma"
  alt[scale_rows] <- (par$lower[scale_rows] + par$upper[scale_rows]) / 2
  if (!isTRUE(all.equal(alt, par$start))) starts <- c(starts, list(alt))
  best <- NULL
  for (start in starts) {
    res <- tryCatch(stats::optim(start, fn, method = "L-BFGS-B", lower = par$lower,
                                 upper = par$upper, control = list(maxit = 500)),
                    error = function(e) NULL)
    if (!is.null(res) && is.finite(res$value) && res$value < .msr_joint_penalty &&
        (is.null(best) || res$value < best$value)) best <- res
  }
  if (is.null(best))
    stop("Joint optimization failed: no starting value gave a finite joint log-likelihood. ",
         "Check that each separate fit refits at its own estimates.", call. = FALSE)
  # Native refits converge to a tolerance, so the joint objective carries small
  # numerical noise. L-BFGS-B can then stop with a line-search code at the
  # optimum. Restart once from the solution; an unchanged value is a converged
  # search.
  if (best$convergence %in% c(51L, 52L)) {
    again <- tryCatch(stats::optim(best$par, fn, method = "L-BFGS-B", lower = par$lower,
                                   upper = par$upper, control = list(maxit = 500)),
                      error = function(e) NULL)
    if (!is.null(again) && is.finite(again$value) && again$value <= best$value + 1e-6) {
      again$message <- paste0("restarted after line-search stop: ", best$message)
      if (again$convergence %in% c(51L, 52L) && abs(again$value - best$value) < 1e-6)
        again$convergence <- 0L
      best <- again
    }
  }
  if (best$convergence != 0L)
    warning("The joint optimizer reported convergence code ", best$convergence, " (",
            best$message, "). Check the estimates against the separate fits.", call. = FALSE)
  best
}


.msr_joint_at_bound <- function(theta, par) {
  tol <- 1e-5
  (is.finite(par$lower) & theta - par$lower < tol) |
    (is.finite(par$upper) & par$upper - theta < tol)
}


# Outer Hessian on the optimizer scale. Parameters on a bound are held fixed;
# the remaining block gives SEs conditional on them.
.msr_joint_hessian <- function(fn, theta, par) {
  p <- length(theta)
  at_bound <- .msr_joint_at_bound(theta, par)
  V <- matrix(NA_real_, p, p)
  free <- which(!at_bound)
  if (length(free)) {
    f_free <- function(x) { th <- theta; th[free] <- x; fn(th) }
    H <- tryCatch(stats::optimHess(theta[free], f_free), error = function(e) NULL)
    if (!is.null(H) && all(is.finite(H)) && max(abs(H)) < 1e30) {
      Vf <- tryCatch(solve(H), error = function(e) NULL)
      if (!is.null(Vf) && all(is.finite(Vf)) && all(diag(Vf) > 0)) V[free, free] <- Vf
    }
  }
  list(vcov = V, at_bound = at_bound)
}


# Profile-likelihood limits for every parameter, on the optimizer scale.
.msr_joint_profile <- function(fn, theta, value, par, level) {
  crit <- stats::qchisq(level, 1)
  p <- length(theta)
  limits <- matrix(NA_real_, p, 2L, dimnames = list(names(theta), c("lower", "upper")))
  limited <- matrix(FALSE, p, 2L, dimnames = dimnames(limits))
  for (i in seq_len(p)) {
    others <- setdiff(seq_len(p), i)
    warm <- theta[others]
    profile_value <- function(psi) {
      th <- theta
      th[i] <- psi
      if (!length(others)) return(fn(th))
      g <- function(x) { th[others] <- x; fn(th) }
      r <- tryCatch(stats::optim(warm, g, method = "L-BFGS-B", lower = par$lower[others],
                                 upper = par$upper[others], control = list(maxit = 200)),
                    error = function(e) NULL)
      if (is.null(r) || !is.finite(r$value)) return(.msr_joint_penalty)
      warm <<- r$par
      r$value
    }
    h <- function(psi) 2 * (profile_value(psi) - value) - crit
    for (side in 1:2) {
      warm <- theta[others]
      end <- if (side == 1L) par$lower[i] else par$upper[i]
      if (!is.finite(end)) {
        step <- 0.5
        end <- theta[[i]] + if (side == 1L) -step else step
        tries <- 0L
        while (h(end) < 0 && tries < 12L) {
          step <- step * 2
          end <- theta[[i]] + if (side == 1L) -step else step
          tries <- tries + 1L
        }
      }
      h_end <- h(end)
      if (!is.finite(h_end) || h_end < 0) {
        limits[i, side] <- end
        limited[i, side] <- TRUE
        next
      }
      warm <- theta[others]
      root <- tryCatch(stats::uniroot(h, sort(c(theta[[i]], end)), tol = 1e-4)$root,
                       error = function(e) NA_real_)
      limits[i, side] <- root
    }
  }
  list(limits = limits, bound_limited = limited)
}


.msr_joint_natural <- function(x, type) {
  ifelse(type == "log_sigma", exp(x), x)
}


.msr_joint_parameter_table <- function(spec, theta, V, prof, z, at_bound) {
  par <- spec$parameters
  se_theta <- sqrt(pmax(diag(V), 0))
  est <- .msr_joint_natural(theta, par$type)
  se <- ifelse(par$type == "log_sigma", est * se_theta, se_theta)
  out <- data.frame(parameter = par$parameter, type = par$type, covariate = par$covariate,
                    data_set = par$data_set, estimate = unname(est), SE = unname(se),
                    wald_lower = unname(.msr_joint_natural(theta - z * se_theta, par$type)),
                    wald_upper = unname(.msr_joint_natural(theta + z * se_theta, par$type)),
                    stringsAsFactors = FALSE)
  out$wald_lower[par$type == "kappa"] <- pmax(out$wald_lower[par$type == "kappa"], 0)
  if (!is.null(prof)) {
    out$profile_lower <- unname(.msr_joint_natural(prof$limits[, 1], par$type))
    out$profile_upper <- unname(.msr_joint_natural(prof$limits[, 2], par$type))
    out$bound_limited <- unname(apply(prof$bound_limited, 1L, any))
  }
  out$at_bound <- unname(at_bound)
  out
}


# Apparent scale in each data set with delta-method SEs. Profile limits carry
# over when the apparent scale is a monotone function of one parameter.
.msr_joint_scale_table <- function(spec, theta, V, prof, z) {
  flat <- function(th) unlist(.msr_joint_sigma(th, spec), use.names = FALSE)
  sigma <- .msr_joint_sigma(theta, spec)
  keys <- do.call(rbind, lapply(names(sigma), function(d)
    data.frame(data_set = d, covariate = names(sigma[[d]]), stringsAsFactors = FALSE)))
  value <- flat(theta)
  J <- matrix(0, length(value), length(theta))
  for (i in seq_along(theta)) {
    h <- 1e-6 * max(1, abs(theta[[i]]))
    up <- dn <- theta
    up[i] <- up[i] + h
    dn[i] <- dn[i] - h
    J[, i] <- (flat(up) - flat(dn)) / (2 * h)
  }
  Vz <- V
  Vz[is.na(Vz)] <- 0
  used <- colSums(abs(J)) > 0
  se <- sqrt(pmax(rowSums((J %*% Vz) * J), 0))
  # A scale that depends on a parameter without a variance has no SE.
  undefined <- as.vector(abs(J[, used, drop = FALSE]) %*%
                           as.numeric(is.na(diag(V))[used]) > 0)
  se[undefined] <- NA_real_
  out <- keys
  out$sigma <- value
  out$SE <- se
  out$wald_lower <- value * exp(-z * se / value)
  out$wald_upper <- value * exp(z * se / value)
  out$footprint <- if (is.null(spec$tau)) 0 else unname(spec$tau[out$data_set])
  out$structure <- unname(spec$structure[out$covariate])
  if (!is.null(prof)) {
    out$profile_lower <- NA_real_
    out$profile_upper <- NA_real_
    par <- spec$parameters
    for (r in seq_len(nrow(out))) {
      j <- out$covariate[r]
      d <- out$data_set[r]
      st <- out$structure[r]
      if (st == "footprint_scaled") next
      i <- if (st == "separate") which(par$type == "log_sigma" & par$covariate == j &
                                          par$data_set %in% d) else
        which(par$type == "log_sigma" & par$covariate == j)
      if (length(i) != 1L) next
      lim <- exp(prof$limits[i, ])
      add <- if (st == "footprint") spec$tau[[d]]^2 else 0
      out$profile_lower[r] <- sqrt(lim[[1L]]^2 + add)
      out$profile_upper[r] <- sqrt(lim[[2L]]^2 + add)
    }
  }
  rownames(out) <- NULL
  out
}


.msr_joint_structure_table <- function(spec) {
  data.frame(covariate = names(spec$structure), structure = unname(spec$structure),
             data_sets = vapply(spec$members[names(spec$structure)], paste, character(1),
                                collapse = ", "),
             shared_coefficient = names(spec$structure) %in% spec$shared,
             stringsAsFactors = FALSE)
}


# Maximize the separate-scales model itself, from the separate fits and from
# the joint estimates, so a separate fit stopped at a poor optimum cannot make
# the joint structure look better than the model that nests it.
.msr_joint_separate_reference <- function(spec, sigma_joint) {
  sets <- spec$sets
  for (d in names(sets)) sets[[d]]$shared <- sets[[d]]$offset_model <- NULL
  structure <- stats::setNames(rep("separate", length(spec$structure)), names(spec$structure))
  sep <- tryCatch(.msr_joint_spec(sets, structure, NULL, character()), error = function(e) NULL)
  if (is.null(sep)) return(NULL)
  fn <- .msr_joint_objective(sep)
  par <- sep$parameters
  from_joint <- vapply(seq_len(nrow(par)), function(i)
    log(sigma_joint[[par$data_set[i]]][[par$covariate[i]]]), numeric(1))
  from_joint <- pmin(pmax(from_joint, par$lower), par$upper)
  best <- NULL
  for (start in list(par$start, from_joint)) {
    res <- tryCatch(stats::optim(start, fn, method = "L-BFGS-B", lower = par$lower,
                                 upper = par$upper, control = list(maxit = 500)),
                    error = function(e) NULL)
    if (!is.null(res) && is.finite(res$value) && res$value < .msr_joint_penalty &&
        (is.null(best) || res$value < best$value)) best <- res
  }
  if (is.null(best)) return(NULL)
  sigma <- .msr_joint_sigma(best$par, sep)
  # Per-data-set log-likelihoods at the refitted separate scales.
  ll_set <- vapply(sets, function(s) -suppressWarnings(
    .msr_joint_set_value(s, sigma[[s$name]], numeric(), NULL)), numeric(1))
  list(logLik = -best$value, logLik_set = ll_set, sigma = sigma)
}


.msr_joint_comparison <- function(spec, ll, k, reference = NULL) {
  ll_input <- vapply(spec$sets, `[[`, numeric(1), "loglik_sep")
  ll_sep <- sum(ll_input)
  improved <- NULL
  if (!is.null(reference) && reference$logLik > ll_sep) {
    ll_sep <- reference$logLik
    gain <- reference$logLik_set - ll_input
    up <- names(gain)[gain > 1e-3 * pmax(1, abs(ll_input))]
    if (length(up))
      improved <- data.frame(
        data_set = up, logLik_input = unname(ll_input[up]),
        logLik_refit = unname(reference$logLik_set[up]),
        sigma = vapply(up, function(d) paste(sprintf("%s = %.4g", names(reference$sigma[[d]]),
                                                     reference$sigma[[d]]), collapse = ", "),
                       character(1)),
        stringsAsFactors = FALSE)
  }
  k_sep <- sum(vapply(spec$sets, `[[`, numeric(1), "k_sep"))
  df <- k_sep - k
  stat <- 2 * (ll_sep - ll)
  boundary <- any(spec$structure == "footprint_scaled")
  p <- if (df > 0) stats::pchisq(max(stat, 0), df, lower.tail = FALSE) else NA_real_
  list(logLik = c(joint = ll, separate = ll_sep), df = c(joint = k, separate = k_sep),
       statistic = max(stat, 0), df_difference = df, p_value = p,
       improved = improved, boundary = boundary)
}


# ---- Components ----------------------------------------------------------------

# A separate fit moved to the joint scales, so the existing plot, map and
# prediction tools apply to each data set.
.msr_joint_component <- function(s, sigma_d, se_d, beta, std, info) {
  final <- .msr_joint_set_value(s, sigma_d, beta, std, return_model = TRUE)
  x <- s$fit
  covs <- s$covs
  n <- length(covs)
  x$scale_est <- data.frame(Mean = unname(sigma_d[covs]), SE = unname(se_d[covs]),
                            row.names = covs)
  x$optim_results$par <- unname(sigma_d[covs]) / s$unit
  x$optim_results$par_unscale <- unname(sigma_d[covs])
  # The separate fit's Hessian no longer describes these scales.
  x$optim_results$hessian <- x$optim_results$hessian_unscale <- matrix(NA_real_, n, n)
  x$opt_mod <- .msr_native_fixed_metadata(final$mod)
  x$scl_params <- .msr_merge_scl_params(primary = final$scl_params,
                                        fallback = s$ki$scl_params)
  x$profile_scale_est <- NULL
  x$diagnostics$refit_hessian <- list(code = "refit_hessian", triggered = FALSE,
                                      status = "joint", uncertainty = "joint")
  x$diagnostics$sigma_precision <- NULL
  x$diagnostics$start_stall <- NULL
  x$warn_message <- 0
  x$joint <- info
  scale_D <- tryCatch(kernel_dist(x, dimension = "1d"), error = function(e) NULL)
  if (!is.null(scale_D))
    x$diagnostics$max_distance <- .max_distance_diagnostic(scale_D = scale_D,
                                                           max_D = s$max_D)
  x
}


# ---- Shared sites --------------------------------------------------------------

# Correlation between data sets of mean standardized residuals at clusters that
# both data sets sampled. Large positive values mean shared unexplained site
# effects, which the independence assumption ignores.
.msr_joint_residual_correlation <- function(components, clusters) {
  res <- lapply(names(components), function(d) {
    r <- .msr_joint_site_residuals(components[[d]]$opt_mod)
    if (is.null(r) || length(r) != length(clusters[[d]])) return(NULL)
    tapply(r, clusters[[d]], mean, na.rm = TRUE)
  })
  names(res) <- names(components)
  pairs <- utils::combn(names(components), 2L, simplify = FALSE)
  rows <- lapply(pairs, function(pr) {
    a <- res[[pr[1L]]]
    b <- res[[pr[2L]]]
    shared <- if (is.null(a) || is.null(b)) character() else intersect(names(a), names(b))
    r <- if (length(shared) >= 5L) suppressWarnings(stats::cor(a[shared], b[shared],
                                                               method = "spearman")) else NA_real_
    data.frame(data_set_1 = pr[1L], data_set_2 = pr[2L], shared_clusters = length(shared),
               spearman = r, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}


.msr_joint_site_residuals <- function(model) {
  model <- .analysis_model(model)
  if (methods::is(model, "unmarkedFit")) {
    r <- tryCatch(unmarked::residuals(model), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    r <- as.matrix(r)
    return(rowMeans(r, na.rm = TRUE) / stats::sd(r, na.rm = TRUE))
  }
  tryCatch(stats::residuals(model, type = "pearson"), error = function(e) NULL)
}


# Delete-a-group jackknife over cluster groups. Each replicate refits every
# data set's starting model without the deleted rows, then refits the joint
# model from the full-data estimates.
.msr_joint_jackknife <- function(spec, theta, clusters, groups, seed) {
  unmarked <- vapply(spec$sets, function(s) identical(s$ctx$mod_class, "unmarked"),
                     logical(1))
  if (any(unmarked))
    stop("uncertainty = 'jackknife' is not yet available for unmarked models (",
         paste(names(spec$sets)[unmarked], collapse = ", "), ").", call. = FALSE)
  ids <- sort(unique(unlist(clusters, use.names = FALSE)))
  g <- min(as.integer(groups), length(ids))
  if (g < 2L) stop("The jackknife needs at least two clusters.", call. = FALSE)
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
      assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  assignment <- stats::setNames(sample(rep_len(seq_len(g), length(ids))), ids)
  par <- spec$parameters
  estimates <- matrix(NA_real_, g, length(theta), dimnames = list(NULL, names(theta)))
  for (k in seq_len(g)) {
    drop_ids <- names(assignment)[assignment == k]
    reduced <- spec
    ok <- TRUE
    for (d in names(spec$sets)) {
      s <- spec$sets[[d]]
      keep <- !clusters[[d]] %in% drop_ids
      if (all(keep)) next
      if (sum(keep) < 3L) { ok <- FALSE; break }
      dat <- s$ctx$data_template[keep, , drop = FALSE]
      model <- tryCatch(.refit_model(s$model, data = dat, opt_context = s$ctx),
                        error = function(e) NULL)
      if (is.null(model)) { ok <- FALSE; break }
      ctx <- tryCatch(build_opt_context(fitted_mod = model, cov_df = s$ki$raw_cov,
                                        join_by = s$join_by, refit_fn = s$ctx$refit_fn,
                                        scale_vars = s$ki$scale_vars, unit_conv = s$ki$unit_conv,
                                        resolution = s$ki$resolution, n_cols = s$ki$n_cols,
                                        binned = s$ki$binned, d_list = s$ki$d_list,
                                        surface_binned = s$ki$surface_binned),
                      error = function(e) NULL)
      if (is.null(ctx)) { ok <- FALSE; break }
      reduced$sets[[d]]$model <- model
      reduced$sets[[d]]$ctx <- ctx
      if (!is.null(s$shared))
        reduced$sets[[d]]$offset_model <- tryCatch(
          .msr_joint_offset_model(reduced$sets[[d]], s$shared), error = function(e) NULL)
    }
    if (!ok) next
    reduced$parameters$start <- theta
    fit <- tryCatch(suppressWarnings(.msr_joint_optimize(reduced, .msr_joint_objective(reduced))),
                    error = function(e) NULL)
    if (!is.null(fit) && fit$convergence != 0L) fit <- NULL
    if (!is.null(fit)) estimates[k, ] <- fit$par
  }
  complete <- stats::complete.cases(estimates)
  m <- sum(complete)
  if (m < 2L) stop("Fewer than two jackknife replicates converged.", call. = FALSE)
  est <- estimates[complete, , drop = FALSE]
  centered <- sweep(est, 2L, colMeans(est))
  V <- (m - 1) / m * crossprod(centered)
  list(vcov = V, estimates = estimates, groups = assignment, replicates = m)
}
