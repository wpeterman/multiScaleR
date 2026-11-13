# Build the joint grid in model order so exact ties resolve reproducibly.
.msr_radius_grid_points <- function(radius_grid, hard, min_D, max_D) {
  if (is.numeric(radius_grid) && length(hard) == 1L)
    radius_grid <- stats::setNames(list(radius_grid), hard)
  if (!is.list(radius_grid) || is.data.frame(radius_grid) ||
      is.null(names(radius_grid)) || anyDuplicated(names(radius_grid)) ||
      !setequal(names(radius_grid), hard))
    stop("Name one radius_grid vector for every optimized hard-radius predictor.", call. = FALSE)
  candidates <- lapply(radius_grid[hard], function(x) {
    if (!is.numeric(x) || !is.null(dim(x)) || length(unique(x)) < 2L ||
        any(!is.finite(x)) || any(x < min_D | x > max_D))
      stop("Each radius_grid vector needs at least two distinct finite radii within min_D and max_D (in map units).", call. = FALSE)
    sort(unique(as.numeric(x)))
  })
  # expand.grid() varies the first column fastest. Reverse that order so
  # the first hard predictor has priority when likelihoods tie exactly.
  # Pass columns without names, which can otherwise collide with formals.
  grid <- expand.grid(unname(rev(candidates)), KEEP.OUT.ATTRS = FALSE)
  grid <- grid[rev(seq_along(candidates))]
  names(grid) <- hard
  grid
}

.msr_refit_radius_grid <- function(radius_grid, r_vars, inputs, model,
                                   join_by, refit_fn, n_cores, PSOCK,
                                   start_strategy) {
  hard <- intersect(r_vars, .msr_hard_radius_covariates(inputs$scale_vars,
    inputs$kernel, r_vars))
  if (!length(hard))
    stop("radius_grid requires at least one optimized hard-radius predictor.", call. = FALSE)
  grid <- .msr_radius_grid_points(radius_grid, hard, inputs$min_D, inputs$max_D)
  if (!is.null(join_by) || !is.null(refit_fn) ||
      (!is.null(n_cores) && n_cores != 1L) || PSOCK || start_strategy != "single")
    stop("radius_grid uses serial native refits with no join_by, refit_fn or start screening.", call. = FALSE)

  native <- .analysis_model(model)
  if (inherits(native, "coxph")) {
    .msr_choice_native_frame(native)
  } else if (!class(native)[1L] %in% c("lm", "glm", "negbin") ||
             !is.data.frame(native$model) || !is.null(native$call$subset)) {
    stop("radius_grid requires a retained native lm, glm, glm.nb or complete one-case clogit model frame, with no subset argument.", call. = FALSE)
  }
  ids <- rownames(inputs$kernel_dat)
  if (!.msr_valid_observation_ids(ids, length(inputs$d_list)) ||
      !identical(names(inputs$d_list), ids) || !identical(names(inputs$raw_cov), ids) ||
      !all(rownames(native$model) %in% ids))
    stop("radius_grid requires matching model row IDs and named prepared raw-cell buffers.", call. = FALSE)
  grid
}

# Bound-aware finite differences for the native likelihood on log parameters.
.msr_radius_gradient <- function(theta, objective, lower, upper, step = 1e-4) {
  step <- rep_len(step, length(theta))
  vapply(seq_along(theta), function(i) {
    left <- right <- theta
    left[i] <- max(theta[i] - step[i], lower[i])
    right[i] <- min(theta[i] + step[i], upper[i])
    (objective(right) - objective(left)) / (right[i] - left[i])
  }, numeric(1))
}

.msr_radius_profile_hessian <- function(par, smooth, lower, upper, evaluate) {
  hessian <- matrix(NA_real_, length(par), length(par))
  theta <- log(par[smooth])
  # Match the refit Wald boundary threshold and stay inside the fitted bounds.
  room <- pmin(theta - lower, upper - theta)
  if (!length(theta) || any(room <= 1e-5))
    return(hessian)
  step <- pmin(2e-4, room / 4)
  objective <- function(log_par) {
    par[smooth] <- exp(log_par)
    -evaluate(par)$log_likelihood
  }
  block <- tryCatch({
    log_hessian <- stats::optimHess(theta, objective,
      control = list(ndeps = step))
    log_gradient <- .msr_radius_gradient(theta, objective, lower, upper, step)
    # Remove the first-derivative contribution before changing coordinates.
    log_block <- log_hessian - diag(log_gradient, nrow = length(theta))
    log_block <- (log_block + t(log_block)) / 2
    # Nested differences lose precision through cancellation. Bound the
    # roundoff scale in log coordinates using the smallest actual step.
    resolution <- 32 * length(theta) * .Machine$double.eps *
      max(1, abs(objective(theta))) / min(step)^2
    curvature <- eigen(log_block, symmetric = TRUE, only.values = TRUE)$values
    if (min(curvature) <= resolution) NULL else log_block / tcrossprod(par[smooth])
  }, error = function(e) NULL)
  if (!is.null(block) && all(is.finite(block)))
    hessian[smooth, smooth] <- (block + t(block)) / 2
  hessian
}

.msr_refit_radius_profile <- function(grid, inputs, context, model, par,
                                      lower, upper, user_supplied_par) {
  covariates <- context$covs
  hard <- names(grid)
  radius_index <- match(hard, covariates)
  n <- length(covariates)
  has_shape <- inputs$kernel == "expow"
  is_smooth <- !covariates %in% hard
  smooth <- which(rep(is_smooth, length.out = length(par)))
  parameter_names <- c(paste0("sigma:", covariates),
    if (has_shape) paste0("shape:", covariates))
  report_scale <- c(rep(inputs$unit_conv, n), if (has_shape) rep(1, n))

  # A hard radius has no shape. Keep the native refit parameter layout.
  if (has_shape) par[n + radius_index] <- 1
  if (length(smooth)) {
    if (!user_supplied_par) {
      prepared <- match(covariates, .msr_optimized_covariates(inputs$scale_vars,
        cov_df = inputs$raw_cov))
      initial <- c(inputs$sigma[prepared], if (has_shape) inputs$shape[prepared])
      if (length(initial) == length(par) && all(is.finite(initial[smooth])))
        par[smooth] <- initial[smooth]
    } else if (any(!is.finite(par[smooth]) | par[smooth] < lower[smooth] |
                   par[smooth] > upper[smooth])) {
      stop("Smooth par starting values must be finite and within the optimization bounds.", call. = FALSE)
    }
    lo <- log(lower[smooth])
    hi <- log(upper[smooth])
    first <- log(pmin(pmax(par[smooth], lower[smooth]), upper[smooth]))
    starts <- unique(rbind(first, lo + .25 * (hi - lo), lo + .75 * (hi - lo)))
    colnames(starts) <- parameter_names[smooth]
  }

  native <- .analysis_model(model)
  choice <- inherits(native, "coxph")
  target <- if (choice) .msr_choice_native_frame(native)$targets else native$model
  values <- function(p, vars) .msr_kernel_cov_w(p, inputs$d_list, inputs$raw_cov,
    inputs$kernel, context, covariates = vars)
  likelihood_values <- function(x) {
    # Row IDs are audited separately. Native recovery can change names and
    # discard scale() provenance without changing a model's actual inputs.
    if (is.null(x)) return(NULL)
    x <- unname(x)
    attr(x, "scaled:center") <- attr(x, "scaled:scale") <- NULL
    x
  }
  audited <- function(p, hard_values) {
    data <- if (any(is_smooth)) cbind(hard_values, values(p, covariates[is_smooth])) else hard_values
    result <- kernel_scale_fn(p, inputs$d_list, inputs$raw_cov, inputs$kernel,
      model, mod_return = TRUE, opt_context = context,
      cov_w = data[, covariates, drop = FALSE])
    fit <- .analysis_model(result$mod)
    if (choice) {
      .msr_choice_identity(fit, list(analysis = target$row_id), target, native$method)
    } else {
      same <- function(get) identical(likelihood_values(get(fit$model)),
        likelihood_values(get(target)))
      if (!identical(rownames(fit$model), rownames(target)) ||
          !same(stats::model.response) || !same(stats::model.weights))
        stop("The refit changed the modeled rows, responses or weights.")
    }
    if (!identical(likelihood_values(stats::model.offset(fit$model)),
                   likelihood_values(stats::model.offset(native$model))))
      stop("The refit changed the model offset.")
    ll <- .msr_model_loglik(fit)
    if (length(ll) != 1L || !is.finite(ll) ||
        any(!is.finite(stats::coef(fit))) ||
        (inherits(fit, "glm") && !isTRUE(fit$converged)))
      stop("The native refit did not produce a finite, converged likelihood and coefficient vector.")
    list(result = result, log_likelihood = ll, par = p)
  }
  strict <- function(expr) withCallingHandlers(expr,
    warning = function(w) stop(conditionMessage(w), call. = FALSE))

  profile <- grid
  names(profile) <- if (ncol(grid) == 1L) "radius" else paste0("radius:", hard)
  for (nm in parameter_names[smooth]) profile[[nm]] <- NA_real_
  profile$log_likelihood <- NA_real_
  profile$status <- "failed"
  profile$message <- NA_character_
  profile$selected <- FALSE
  if (length(smooth)) profile$smooth_starts <- rep(list(NULL), nrow(grid))
  selected <- NA_integer_
  best_loglik <- -Inf
  evaluations <- 0L
  for (i in seq_len(nrow(grid))) {
    candidate <- tryCatch(strict({
      p <- par
      p[radius_index] <- unlist(grid[i, , drop = FALSE], use.names = FALSE) /
        inputs$unit_conv
      # Cache only this combination's fixed metric columns. Each native refit
      # recomputes smooth predictors and estimates all native coefficients.
      fixed <- values(p, hard)
      validate_covariates_before_scale(fixed, context = "optimized covariates",
        scale_vars = context$scale_vars)
      if (anyNA(fixed))
        stop("A hard-radius predictor is missing for an observation in the fitted sample.")
      if (length(smooth)) {
        evaluate <- function(theta) {
          p[smooth] <- exp(theta)
          evaluations <<- evaluations + 1L
          strict(audited(p, fixed))
        }
        fn <- function(theta) tryCatch(-evaluate(theta)$log_likelihood,
          error = function(e) 1e60)
        gr <- function(theta) .msr_radius_gradient(theta, fn, lo, hi)
        cold <- function(theta) list(value = -evaluate(theta)$log_likelihood,
          gradient = gr(theta))
        optim <- .msr_multistart(list(fn = fn, gr = gr, cold_evaluate = cold),
          starts, lo, hi, gradient_tolerance = max(1e-3, nrow(native$model) * 1e-5))
        records <- optim$records[setdiff(names(optim$records), "seconds")]
        records[paste0("initial:", parameter_names[smooth])] <-
          as.data.frame(sweep(exp(starts), 2L, report_scale[smooth], "*"))
        profile$smooth_starts[i] <- list(records)
        if (is.null(optim$best))
          stop("No smooth optimization start passed convergence and native likelihood checks.")
        result <- evaluate(optim$best$par)
        if (abs(result$log_likelihood + optim$best$objective) > 1e-5)
          stop("The final native likelihood disagrees with the smooth optimization objective.")
        result
      } else {
        evaluations <- evaluations + 1L
        audited(p, fixed)
      }
    }), error = function(e) e)
    if (inherits(candidate, "error")) {
      profile$message[i] <- conditionMessage(candidate)
      next
    }
    if (length(smooth))
      profile[i, parameter_names[smooth]] <- as.list(candidate$par[smooth] * report_scale[smooth])
    profile$log_likelihood[i] <- candidate$log_likelihood
    profile$status[i] <- "scored"
    if (candidate$log_likelihood > best_loglik) {
      selected <- i
      best <- candidate
      best_fixed <- fixed
      best_loglik <- candidate$log_likelihood
    }
  }
  if (is.na(selected))
    stop(structure(list(message = "No radius_grid candidate produced an eligible native refit.",
      call = NULL, radius_profile = profile),
      class = c("msr_radius_profile_error", "error", "condition")))
  profile$selected[selected] <- TRUE
  hessian <- .msr_radius_profile_hessian(best$par, smooth,
    log(lower[smooth]), log(upper[smooth]), function(p) strict(audited(p, best_fixed)))
  list(par = best$par, value = -best_loglik, convergence = 0L,
    counts = c("function" = evaluations, gradient = NA_integer_),
    message = "Evaluated the declared hard-radius grid.",
    hessian = hessian, radius_profile = profile, final_mod = best$result)
}
