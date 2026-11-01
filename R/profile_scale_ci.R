#' @title Profile-Likelihood Helpers for Scale Parameters
#' @description Internal helpers for estimating profile-likelihood intervals
#'   for `sigma` while re-optimizing nuisance parameters.
#' @keywords internal
#' @noRd

.profile_scale_cache <- new.env(parent = emptyenv())

.msr_profile_cache_bytes <- function() {
  keys <- ls(.profile_scale_cache, all.names = TRUE)
  sum(vapply(keys, function(key) as.numeric(utils::object.size(
    get(key, envir = .profile_scale_cache, inherits = FALSE))), numeric(1)))
}

.msr_clear_profile_cache <- function() {
  rm(list = ls(.profile_scale_cache, all.names = TRUE), envir = .profile_scale_cache)
  invisible(NULL)
}

scale_ci_table <- function(object,
                           df,
                           min_D = object$min_D,
                           names = row.names(object$scale_est),
                           profile = FALSE) {

  if (identical(object$backend, "joint_mixed")) {
    if (isTRUE(profile)) stop("Scale profiling is not validated for joint_mixed.", call. = FALSE)
    tab <- ci_func(object$scale_est, df = df, min_D = min_D, names = names)
    attr(tab, "interval_method") <- stats::setNames(rep("unavailable", nrow(tab)), rownames(tab))
    return(tab)
  }

  if (.msr_analytic_backend(object)) {
    tab <- .msr_profile_wald_table(object)
    if (isTRUE(profile)) {
      pr <- .msr_profile_likelihood_intervals(object)
      tab[, c("2.5%", "97.5%")] <- pr$intervals
      attr(tab, "interval_method") <- pr$method
      attr(tab, "bound_limited") <- pr$bound_limited
    }
    return(tab)
  }
  tab <- ci_func(object$scale_est,
                 df = df,
                 min_D = min_D,
                 names = names)

  methods <- ifelse(is.finite(tab$SE), "wald", "unavailable")
  names(methods) <- row.names(tab)

  if (!isTRUE(profile)) {
    attr(tab, "interval_method") <- methods
    return(tab)
  }

  budget <- getOption("multiScaleR.profile_cache_max_bytes", 64 * 1024^2)
  validate_scalar_numeric(budget, "multiScaleR.profile_cache_max_bytes", lower = 0)
  if (.msr_profile_cache_bytes() > budget) .msr_clear_profile_cache()
  cache_enabled <- budget > 0 && is.null(object$opt_context$refit_fn)
  cache_info <- list(status = if (cache_enabled) "enabled" else "disabled",
    reason = if (budget == 0) "zero_budget" else if (!cache_enabled) "custom_refit" else NULL)
  cache_key <- NULL
  entry <- NULL
  if (cache_enabled) {
    cache_key <- tryCatch(profile_scale_cache_key(object, min_D, names),
      error = function(e) {
        cache_info <<- list(status = "unavailable", reason = conditionMessage(e))
        warning("Profile cache unavailable; computing intervals without caching: ",
          conditionMessage(e), call. = FALSE)
        NULL
      })
    if (is.null(cache_key)) cache_enabled <- FALSE else {
      # Release oversized key payloads before the numerical profile starts.
      key_bytes <- as.numeric(utils::object.size(list(
        payload = attr(cache_key, "payload"), result = NULL)))
      if (key_bytes > budget) {
        cache_enabled <- FALSE
        cache_info <- list(status = "disabled", reason = "payload_exceeds_budget")
        cache_key <- NULL
      }
    }
  }
  if (cache_enabled && exists(cache_key, envir = .profile_scale_cache, inherits = FALSE))
    entry <- get(cache_key, envir = .profile_scale_cache, inherits = FALSE)
  # Verify the full payload as well as the digest; cache corruption is a miss.
  if (is.list(entry) && identical(entry$payload, attr(cache_key, "payload")) &&
      is.list(entry$result) && !is.null(entry$result$intervals)) {
    profile_tab <- entry$result
  } else {
    profile_tab <- try(profile_sigma_intervals(object), silent = TRUE)
    if (inherits(profile_tab, "try-error") || is.null(profile_tab)) {
      detail <- if (inherits(profile_tab, "try-error")) {
        conditionMessage(attr(profile_tab, "condition"))
      } else {
        "The fitted object lacks information needed for profiling."
      }
      stop("Profile-likelihood intervals could not be computed: ", detail,
           call. = FALSE)
    }
    if (cache_enabled) {
      candidate <- list(payload = attr(cache_key, "payload"), result = profile_tab)
      bytes <- as.numeric(utils::object.size(candidate))
      if (bytes <= budget) {
        if (.msr_profile_cache_bytes() + bytes > budget) .msr_clear_profile_cache()
        assign(cache_key, candidate, envir = .profile_scale_cache)
      } else cache_info <- list(status = "disabled", reason = "result_exceeds_budget")
    }
  }

  if (!inherits(profile_tab, "try-error") && !is.null(profile_tab$intervals)) {
    common_rows <- intersect(row.names(tab), row.names(profile_tab$intervals))
    if (length(common_rows) > 0) {
      tab[common_rows, c("2.5%", "97.5%")] <-
        profile_tab$intervals[common_rows, c("2.5%", "97.5%")]
      methods[common_rows] <- profile_tab$method[common_rows]
    }
  }

  attr(tab, "interval_method") <- methods
  attr(tab, "profile_cache") <- cache_info
  tab
}


profile_scale_cache_key <- function(object,
                                    min_D = object$min_D,
                                    names = row.names(object$scale_est)) {

  key_payload <- list(
    call = deparse(object$call, width.cutoff = 500L),
    kernel = object$kernel_inputs$kernel,
    min_D = min_D,
    max_D = object$max_D,
    unit_conv = object$kernel_inputs$unit_conv,
    scale_est = unclass(object$scale_est),
    optim_par = object$optim_results$par,
    optim_value = object$optim_results$value,
    names = names
  )

  # Data and refit controls can change while the call and optimum stay equal.
  model <- object$fitted_mod_original
  frame <- .msr_selection_frame(.analysis_model(model))
  if (!is.null(frame)) attr(frame, "terms") <- NULL
  key_payload$training_inputs <- object$kernel_inputs
  key_payload$training_frame <- frame
  key_payload$response <- if (is.list(model)) model$y else NULL
  key_payload$weights <- if (is.list(model)) model$prior.weights else NULL
  key_payload$join_by <- object$join_by
  key_payload$model_class <- class(model)
  key_payload$family <- tryCatch(stats::family(model)[c("family", "link")],
    error = function(e) NULL)
  key_payload$native_data <- if (methods::is(model, "unmarkedFit")) model@data else NULL
  context <- object$opt_context
  context$fitted_mod <- context$refit_fn <- NULL
  key_payload$context <- context
  # Formula environments may contain the caller's mutable workspace. Do not
  # serialize them: the retained numeric frame and prepared cells bind the data.
  key_payload$formula <- tryCatch(paste(deparse(stats::formula(model)), collapse = " "),
    error = function(e) NULL)
  payload <- serialize(key_payload, NULL, version = 2)
  path <- tempfile("msr-profile-key-")
  on.exit(unlink(path), add = TRUE)
  writeBin(payload, path)
  structure(unname(tools::md5sum(path)), payload = payload)
}


profile_sigma_intervals <- function(object,
                                    level = 0.95,
                                    max_expand = 5,
                                    refine_steps = 4) {

  if (.msr_analytic_backend(object))
    return(.msr_profile_likelihood_intervals(object, level))
  if (!inherits(object, "multiScaleR")) {
    stop("`object` must inherit from class 'multiScaleR'.", call. = FALSE)
  }

  if (is.null(object$optim_results$par) ||
      is.null(object$max_D) ||
      is.null(object$opt_context) ||
      is.null(object$fitted_mod_original)) {
    return(NULL)
  }

  n_covs <- nrow(object$scale_est)
  if (is.null(n_covs) || n_covs < 1) {
    return(NULL)
  }

  bounds <- parameter_bounds(object)
  full_par <- object$optim_results$par
  # Verify that the stored call can still refit the model before a profile
  # interval is constructed. A lost formula or data symbol must not turn a
  # constant optimizer penalty into an apparently precise confidence interval.
  kernel_scale_fn(par = full_par,
                  d_list = object$kernel_inputs$d_list,
                  cov_df = object$kernel_inputs$raw_cov,
                  kernel = object$kernel_inputs$kernel,
                  fitted_mod = object$fitted_mod_original,
                  join_by = object$join_by,
                  mod_return = TRUE,
                  opt_context = object$opt_context)
  mle_fit <- list(value = object$optim_results$value,
                  par_full = full_par)
  cutoff <- object$optim_results$value + stats::qchisq(level, df = 1) / 2

  intervals <- matrix(NA_real_,
                      nrow = n_covs,
                      ncol = 2,
                      dimnames = list(row.names(object$scale_est),
                                      c("2.5%", "97.5%")))
  method <- rep("wald", n_covs)
  names(method) <- row.names(object$scale_est)

  for (idx in seq_len(n_covs)) {
    prof <- try(
      profile_single_sigma(object = object,
                           fixed_index = idx,
                           mle_fit = mle_fit,
                           lower = bounds$lower,
                           upper = bounds$upper,
                           cutoff = cutoff,
                           max_expand = max_expand,
                           refine_steps = refine_steps),
      silent = TRUE
    )

    if (!inherits(prof, "try-error") && !is.null(prof$interval)) {
      intervals[idx, ] <- prof$interval * object$kernel_inputs$unit_conv
      method[idx] <- prof$method
    }
  }

  list(intervals = as.data.frame(intervals),
       method = method)
}


profile_single_sigma <- function(object,
                                 fixed_index,
                                 mle_fit,
                                 lower,
                                 upper,
                                 cutoff,
                                 max_expand = 5,
                                 refine_steps = 4) {

  mle_sigma <- mle_fit$par_full[fixed_index]
  cache <- new.env(parent = emptyenv())

  eval_point <- function(candidate, start_full) {
    key <- format(signif(candidate, 12), scientific = FALSE, trim = TRUE)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }

    fit <- profile_eval(object = object,
                        fixed_index = fixed_index,
                        fixed_value = candidate,
                        start_full = start_full,
                        lower = lower,
                        upper = upper)
    assign(key, fit, envir = cache)
    fit
  }

  lower_side <- profile_sigma_side(direction = "lower",
                                   mle_sigma = mle_sigma,
                                   boundary = lower[fixed_index],
                                   cutoff = cutoff,
                                   mle_fit = mle_fit,
                                   eval_point = eval_point,
                                   max_expand = max_expand,
                                   refine_steps = refine_steps)

  upper_side <- profile_sigma_side(direction = "upper",
                                   mle_sigma = mle_sigma,
                                   boundary = upper[fixed_index],
                                   cutoff = cutoff,
                                   mle_fit = mle_fit,
                                   eval_point = eval_point,
                                   max_expand = max_expand,
                                   refine_steps = refine_steps)

  list(
    interval = c(lower_side$endpoint, upper_side$endpoint),
    method = if (identical(lower_side$method, "profile") &&
                 identical(upper_side$method, "profile")) {
      "profile"
    } else {
      "profile_bound"
    }
  )
}


profile_sigma_side <- function(direction,
                               mle_sigma,
                               boundary,
                               cutoff,
                               mle_fit,
                               eval_point,
                               max_expand = 5,
                               refine_steps = 4,
                               factor = 1.5) {

  tol <- sqrt(.Machine$double.eps)
  if (abs(mle_sigma - boundary) <= tol) {
    return(list(endpoint = boundary, method = "profile_bound"))
  }

  inside_value <- mle_sigma
  inside_fit <- mle_fit
  candidate <- next_profile_point(mle_sigma, boundary, direction, factor)

  if (abs(candidate - inside_value) <= tol) {
    return(list(endpoint = boundary, method = "profile_bound"))
  }

  outside_value <- NA_real_
  outside_fit <- NULL

  for (step in seq_len(max_expand)) {
    fit <- eval_point(candidate, inside_fit$par_full)

    if (!is.finite(fit$value) || fit$value >= cutoff) {
      outside_value <- candidate
      outside_fit <- fit
      break
    }

    inside_value <- candidate
    inside_fit <- fit

    if (abs(candidate - boundary) <= tol) {
      return(list(endpoint = boundary, method = "profile_bound"))
    }

    next_candidate <- next_profile_point(candidate, boundary, direction, factor)
    if (abs(next_candidate - candidate) <= tol) {
      return(list(endpoint = boundary, method = "profile_bound"))
    }
    candidate <- next_candidate
  }

  if (is.null(outside_fit)) {
    boundary_fit <- eval_point(boundary, inside_fit$par_full)
    if (is.finite(boundary_fit$value) && boundary_fit$value >= cutoff) {
      outside_value <- boundary
    } else {
      return(list(endpoint = boundary, method = "profile_bound"))
    }
  }

  for (iter in seq_len(refine_steps)) {
    mid <- midpoint_profile(inside_value, outside_value)
    if (abs(mid - inside_value) <= tol || abs(mid - outside_value) <= tol) {
      break
    }

    mid_fit <- eval_point(mid, inside_fit$par_full)
    if (is.finite(mid_fit$value) && mid_fit$value < cutoff) {
      inside_value <- mid
      inside_fit <- mid_fit
    } else {
      outside_value <- mid
    }
  }

  list(endpoint = midpoint_profile(inside_value, outside_value),
       method = "profile")
}


profile_eval <- function(object,
                         fixed_index,
                         fixed_value,
                         start_full,
                         lower,
                         upper) {

  total_params <- seq_along(start_full)
  free_index <- total_params[-fixed_index]
  profile_warnings <- character()

  objective <- function(par_free) {
    full_par <- inject_fixed_par(start_full = start_full,
                                 free_index = free_index,
                                 fixed_index = fixed_index,
                                 par_free = par_free,
                                 fixed_value = fixed_value)

    value <- try(withCallingHandlers(
      kernel_scale_fn(par = full_par,
                      d_list = object$kernel_inputs$d_list,
                      cov_df = object$kernel_inputs$raw_cov,
                      kernel = object$kernel_inputs$kernel,
                      fitted_mod = object$fitted_mod_original,
                      join_by = object$join_by,
                      opt_context = object$opt_context),
      warning = function(w) {
        profile_warnings <<- c(profile_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ), silent = TRUE)

    if (inherits(value, "try-error") || !is.finite(value)) {
      return(1e20)
    }

    as.numeric(value)[1]
  }

  if (length(free_index) == 0) {
    full_par <- start_full
    full_par[fixed_index] <- fixed_value
    value <- objective(numeric(0))
    return(list(value = value,
                par_full = full_par,
                convergence = 0,
                warnings = unique(profile_warnings)))
  }

  start_free <- pmin(pmax(start_full[free_index], lower[free_index]),
                     upper[free_index])

  fit <- try(
    optim(par = start_free,
          fn = objective,
          lower = lower[free_index],
          upper = upper[free_index],
          method = "L-BFGS-B",
          control = list(maxit = 250)),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    full_par <- inject_fixed_par(start_full = start_full,
                                 free_index = free_index,
                                 fixed_index = fixed_index,
                                 par_free = start_free,
                                 fixed_value = fixed_value)
    return(list(value = 1e20,
                par_full = full_par,
                convergence = 1,
                warnings = unique(profile_warnings)))
  }

  full_par <- inject_fixed_par(start_full = start_full,
                               free_index = free_index,
                               fixed_index = fixed_index,
                               par_free = fit$par,
                               fixed_value = fixed_value)

  list(value = fit$value,
       par_full = full_par,
       convergence = fit$convergence,
       warnings = unique(profile_warnings))
}


inject_fixed_par <- function(start_full,
                             free_index,
                             fixed_index,
                             par_free,
                             fixed_value) {
  out <- start_full
  out[fixed_index] <- fixed_value
  if (length(free_index) > 0) {
    out[free_index] <- par_free
  }
  out
}


parameter_bounds <- function(object) {
  n_covs <- nrow(object$scale_est)
  lower <- rep(object$min_D / object$kernel_inputs$unit_conv, n_covs)
  upper <- rep(object$max_D / object$kernel_inputs$unit_conv, n_covs)

  if (identical(object$kernel_inputs$kernel, "expow")) {
    lower <- c(lower, rep(0.75, n_covs))
    upper <- c(upper, rep(50, n_covs))
  }

  list(lower = lower, upper = upper)
}


next_profile_point <- function(value,
                               boundary,
                               direction,
                               factor = 1.5) {
  if (identical(direction, "lower")) {
    max(boundary, value / factor)
  } else {
    min(boundary, value * factor)
  }
}


midpoint_profile <- function(x, y) {
  if (all(c(x, y) > 0)) {
    exp(mean(log(c(x, y))))
  } else {
    mean(c(x, y))
  }
}
