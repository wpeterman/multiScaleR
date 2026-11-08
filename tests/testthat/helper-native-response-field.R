# One common declared conformance dataset, not a recovery or calibration panel.
native_response_field_fixture <- local({
  common <- NULL
  cache <- list()
  function(family = c("poisson", "gaussian")) {
    family <- match.arg(family)
    if (!is.null(cache[[family]])) return(cache[[family]])
    if (is.null(common)) {
      base <- refit_spatial_fixture()
      locations <- terra::crds(base$points)
      index <- rep(seq_len(nrow(locations)), each = 3L)
      points <- terra::vect(locations[index, , drop = FALSE], crs = terra::crs(base$points))
      inputs <- kernel_prep(points, base$rasters, max_D = 40, kernel = "gaussian",
        scale_vars = msr_vars(mean = kernel_var("hab"),
          rough = surface_var("elev", "sq", weighted = TRUE),
          fixed = surface_var("elev", "sa", radius = 25), unused = kernel_var("hab")),
        bin = FALSE, verbose = FALSE)
      par <- c(19, 29) / inputs$unit_conv
      context <- list(covs = c("mean", "rough"), n_covs = 2L, scale_vars = inputs$scale_vars,
        unit_conv = inputs$unit_conv, resolution = inputs$resolution, n_cols = inputs$n_cols,
        binned = inputs$binned)
      raw <- .msr_kernel_cov_w(par, inputs$d_list, inputs$raw_cov, inputs$kernel, context)
      scaled <- scale(raw)
      data <- inputs$kernel_dat
      data[c("mean", "rough")] <- as.data.frame(scaled)
      data$x <- locations[index, 1L]; data$y_coordinate <- locations[index, 2L]
      data$field_pos <- glmmTMB::numFactor(data$x, data$y_coordinate)
      data$field_group <- factor(rep("shared", nrow(data)))
      data$z <- sin(seq_len(nrow(data)) / 3)
      data$exposure <- 1 + seq_len(nrow(data)) / (2 * nrow(data))
      present <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
      previous <- if (present) get(".Random.seed", .GlobalEnv) else NULL
      on.exit(if (present) assign(".Random.seed", previous, .GlobalEnv) else
        if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
      set.seed(861104)
      covariance <- .6^2 * exp(-as.matrix(dist(locations)) / 55)
      field <- drop(t(chol(covariance)) %*% rnorm(nrow(locations)))
      eta <- .3 * data$mean - .2 * data$rough + .2 * data$z +
        log(data$exposure) + field[index]
      data$gaussian_response <- 1 + eta + rnorm(nrow(data), sd = .2)
      data$poisson_response <- rpois(nrow(data), exp(2 + eta))
      common <<- list(data = data, inputs = inputs, par = par, locations = locations,
        points = points, rasters = base$rasters, scaled = scaled)
    }
    data <- common$data
    data$y <- data[[paste0(family, "_response")]]
    warnings <- character()
    path <- gdr_test_audit_path(paste0("response-field-native-", family, "-"))
    fit <- tryCatch(withCallingHandlers({
      if (family == "poisson") glmmTMB::glmmTMB(
        y ~ mean + rough + z + offset(log(exposure)) + exp(field_pos + 0 | field_group),
        data = data, family = stats::poisson(), start = list(theta = c(log(.6), log(55))),
        control = glmmTMB::glmmTMBControl(parallel = 1L, rank_check = "stop",
          optCtrl = list(iter.max = 1000, eval.max = 1500))) else glmmTMB::glmmTMB(
        y ~ mean + rough + z + offset(log(exposure)) + exp(field_pos + 0 | field_group),
        data = data, family = stats::gaussian(), start = list(theta = c(log(.6), log(55))),
        control = glmmTMB::glmmTMBControl(parallel = 1L, rank_check = "stop",
          optCtrl = list(iter.max = 1000, eval.max = 1500)))
    }, warning = function(w) warnings <<- c(warnings, conditionMessage(w))),
      error = function(e) {
        .msr_gdr_atomic_audit(list(status = "native_constructor_failed", family = family,
          warnings = warnings, error = conditionMessage(e), seed = 861104L), paste0(path, ".rds"))
        stop(e)
      })
    .msr_gdr_atomic_audit(list(status = "native_constructor_returned", family = family,
      warnings = warnings, convergence = fit$fit$convergence, pdHess = fit$sdr$pdHess,
      optimizer = fit$fit[c("par", "objective", "convergence", "message", "evaluations")],
      logLik = as.numeric(logLik(fit)), K = attr(logLik(fit), "df"), seed = 861104L,
      data_sha256 = .msr_gdr_hash(data), fields = tryCatch(.msr_glmmTMB_response_fields(fit),
        error = function(e) list(metadata_error = conditionMessage(e)))), paste0(path, ".rds"))
    object <- structure(list(opt_mod = fit, backend = "refit", kernel_inputs = common$inputs,
      scale_est = data.frame(Mean = common$par * common$inputs$unit_conv,
        SE = NA_real_, row.names = c("mean", "rough")), optim_results = list(par = common$par),
      min_D = common$inputs$min_D, max_D = common$inputs$max_D,
      scl_params = .msr_merge_scl_params(list(mean = attr(common$scaled, "scaled:center"),
        sd = attr(common$scaled, "scaled:scale")), common$inputs$scl_params)), class = "multiScaleR")
    cache[[family]] <<- list(model = fit, object = object, data = data,
      inputs = common$inputs, par = common$par, points = common$points, rasters = common$rasters)
    cache[[family]]
  }
})
