# Marginalize one scale-group SD instead of approximating its posterior with
# a Gaussian. Conditional Laplace volumes account for the other parameters.
.msr_bayes_marginal <- function(out, draws, seed, control) {
  out$inference_metadata <- NULL
  if (!requireNamespace("RTMB", quietly = TRUE))
    stop("Marginal Laplace inference requires optional package 'RTMB'.", call. = FALSE)
  .msr_scale_check_capability(out$specification, "posterior", "marginal_laplace")
  validate_scalar_numeric(draws, "laplace_draws", integerish = TRUE, lower = 2)
  if (!is.null(seed)) validate_scalar_numeric(seed, "seed", integerish = TRUE,
    lower = 1, upper = .Machine$integer.max)
  allowed <- c("starts", "iter.max", "max_grid", "cdf_tol", "quantile_tol", "logz_tol",
               "correction", "importance_min", "importance_max", "importance_df", "importance_refinements", "nested_response", "conditional_beta")
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyDuplicated(names(control)) || any(!names(control) %in% allowed))))
    stop("Unsupported marginal control; see multiScale_bayes help for available settings.", call. = FALSE)
  ctl <- utils::modifyList(list(starts = 3L, iter.max = 1000L, max_grid = 129L,
    cdf_tol = .005, quantile_tol = .02, logz_tol = .01,
    correction = "none", importance_min = 16384L, importance_max = 65536L, importance_df = 7,
    importance_refinements = 1L,
    nested_response = FALSE, conditional_beta = FALSE), control)
  if (!"nested_response" %in% names(control))
    ctl$nested_response <- identical(ctl$correction, "importance") && !isTRUE(ctl$conditional_beta)
  if (!is.character(ctl$correction) || length(ctl$correction) != 1L ||
      !ctl$correction %in% c("none", "importance"))
    stop("correction must be 'none' or 'importance'.", call. = FALSE)
  if (!is.logical(ctl$nested_response) || length(ctl$nested_response) != 1L || is.na(ctl$nested_response))
    stop("nested_response must be TRUE or FALSE.", call. = FALSE)
  if (ctl$nested_response && ctl$correction != "importance")
    stop("nested_response currently requires correction = 'importance'.", call. = FALSE)
  if (!is.logical(ctl$conditional_beta) || length(ctl$conditional_beta) != 1L || is.na(ctl$conditional_beta))
    stop("conditional_beta must be TRUE or FALSE.", call. = FALSE)
  if (ctl$conditional_beta && ctl$correction != "importance")
    stop("conditional_beta requires correction = 'importance'.", call. = FALSE)
  if (ctl$conditional_beta && ctl$nested_response)
    stop("conditional_beta requires nested_response = FALSE; omit nested_response to select this automatically.", call. = FALSE)
  if (length(ctl$importance_max) != 1L || !ctl$importance_max %in% c(16384L, 32768L, 65536L, 131072L))
    stop("importance_max must be 16384, 32768, 65536, or 131072.", call. = FALSE)
  validate_scalar_numeric(ctl$importance_max, "importance_max", integerish = TRUE, lower = 16384)
  if (length(ctl$importance_min) != 1L || !ctl$importance_min %in% c(16384L, 32768L, 65536L, 131072L))
    stop("importance_min must be 16384, 32768, 65536, or 131072.", call. = FALSE)
  validate_scalar_numeric(ctl$importance_min, "importance_min", integerish = TRUE, lower = 16384)
  if (ctl$importance_min > ctl$importance_max)
    stop("importance_min must not exceed importance_max.", call. = FALSE)
  validate_scalar_numeric(ctl$importance_df, "importance_df", lower = 5, upper = 7)
  if (!ctl$importance_df %in% c(5, 7))
    stop("importance_df must be 5 or 7.", call. = FALSE)
  validate_scalar_numeric(ctl$importance_refinements, "importance_refinements",
    integerish = TRUE, lower = 0, upper = 2)
  if (ctl$correction == "importance" && !requireNamespace("posterior", quietly = TRUE))
    stop("Importance correction requires optional package 'posterior'.", call. = FALSE)
  if (ctl$correction == "importance" && utils::packageVersion("posterior") < "1.6.0")
    stop("Importance correction requires posterior >= 1.6.0 for Pareto diagnostics.", call. = FALSE)
  for (nm in c("starts", "iter.max"))
    validate_scalar_numeric(ctl[[nm]], nm, integerish = TRUE, lower = 1)
  if (length(ctl$max_grid) != 1L || !ctl$max_grid %in% c(33L, 65L, 129L, 257L))
    stop("max_grid must be 33, 65, 129, or 257.", call. = FALSE)
  for (nm in c("cdf_tol", "quantile_tol", "logz_tol"))
    validate_scalar_numeric(ctl[[nm]], nm, lower = 0, inclusive_lower = FALSE)
  spec <- out$specification
  obj <- RTMB::MakeADFun(.msr_bayes_objective(spec, ad = TRUE),
                         spec$parameters, silent = TRUE)
  conditional <- obj
  if (ctl$nested_response && length(spec$parameters$response_z)) {
    conditional <- RTMB::MakeADFun(.msr_bayes_objective(spec, ad = TRUE),
      spec$parameters, random = "response_z", silent = TRUE)
    # RTMB does not provide an analytic Hessian after random-effect integration.
    # Differentiate its analytic outer gradient; do not repair the curvature.
    conditional$he <- function(p) stats::optimHess(p, conditional$fn, conditional$gr)
  }
  hyper <- which(names(obj$par) == "log_scale_sd")
  response_factorized <- ctl$nested_response && ctl$correction == "importance" &&
    length(spec$response$random) == 1L
  fixed_indices <- which(names(obj$par) != "response_z")
  grid_hyper <- if (response_factorized) which(names(conditional$par) == "log_scale_sd") else hyper
  prior_sd <- spec$priors$scale_group_sd
  evaluate <- function(x, warm) {
    tau <- if (x == 0) prior_sd * 1e-8 else if (x == 1) prior_sd * 8 else
      prior_sd * stats::qnorm((1 + x) / 2)
    if (identical(conditional, obj))
      return(.msr_bayes_marginal_node(obj, hyper, log(tau), prior_sd, warm, ctl))
    if (response_factorized)
      return(.msr_bayes_marginal_node(conditional, grid_hyper, log(tau), prior_sd, warm, ctl))
    fixed <- which(names(obj$par) != "response_z")
    small_warm <- if (is.null(warm)) NULL else list(mode = warm$mode[fixed])
    node <- .msr_bayes_marginal_node(conditional, hyper, log(tau), prior_sd, small_warm, ctl)
    .msr_bayes_nested_expand(node, conditional, obj, hyper)
  }
  grid <- .msr_bayes_marginal_grid(evaluate, prior_sd, ctl)
  if (ctl$correction == "none") {
    samples <- .msr_bayes_marginal_draws(grid, hyper, draws, seed)
    ends <- cumsum(lengths(spec$parameters))
    begins <- ends - lengths(spec$parameters)
    out$draws <- lapply(seq_along(ends), function(j)
      samples[, begins[j] + seq_len(ends[j] - begins[j]), drop = FALSE])
    names(out$draws) <- names(spec$parameters)
  }
  competing <- sum(grid$weights * vapply(grid$nodes, function(z) z$competing, logical(1)))
  upper_mass <- 1 - .msr_bayes_marginal_cdf(grid, .9999)
  out$inference <- "marginal Laplace with scale-SD quadrature"
  out$marginal <- c(grid, list(hyper_index = grid_hyper,
    parameter_indices = if (response_factorized) fixed_indices else seq_along(obj$par),
    draws = draws, seed = seed,
    control = ctl, coordinates = if (response_factorized)
      "Prior-CDF scale SD; response effects integrated; Gaussian outer-parameter nodes" else
      "Prior-CDF quadrature for scale SD; conditional noncentered Gaussian for other parameters"))
  passed <- grid$converged && competing <= .01 && upper_mass <= .001
  out$diagnostics <- list(passed = passed, type = "Marginal Laplace numerical checks",
    grid_converged = grid$converged, grid_nodes = length(grid$x),
    refinement = grid$refinement, competing_mass = competing, upper_tail_mass = upper_mass,
    gradient_max = max(vapply(grid$nodes, function(z) z$gradient_max, 0)),
    condition_number = max(vapply(grid$nodes, function(z) z$condition_number, 0)),
    criteria = paste("Conditional gradients <= 1e-4; positive-definite Hessians; condition <= 1e12;",
      "grid CDF/SD-quantile/log-normalizer convergence; competing-mode mass <= .01;",
      "upper prior-CDF tail mass <= .001. Numerical checks do not establish calibration."))
  if (!passed && ctl$correction == "none")
    warning("Marginal Laplace numerical checks need attention; inspect grid convergence, conditional modes, and tail mass.", call. = FALSE)
  out$inference_metadata <- .msr_scale_inference_metadata("posterior", "marginal_laplace")
  if (ctl$correction == "importance")
    out <- .msr_bayes_importance(out, obj, draws, seed, ctl)
  out
}

.msr_bayes_marginal_node <- function(obj, hyper, log_tau, prior_sd, warm, control) {
  initial <- obj$par
  initial[hyper] <- log_tau
  nuisance <- setdiff(seq_along(initial), hyper)
  unpack <- function(u) { q <- initial; q[nuisance] <- u; q }
  fn <- function(u) obj$fn(unpack(u))
  gr <- function(u) as.numeric(obj$gr(unpack(u)))[nuisance]
  he <- function(u) obj$he(unpack(u))[nuisance, nuisance, drop = FALSE]
  attempts <- lapply(seq_len(control$starts), function(i) {
    p <- initial[nuisance]
    if (i == 1L && !is.null(warm)) p <- warm$mode[nuisance]
    if (i > 1L) {
      p <- p + .1 * sin(seq_along(p) + i)
      p[names(p) == "scale_beta"] <- p[names(p) == "scale_beta"] +
        rep(c(-1, 1, -.5, .5), length.out = control$starts - 1L)[i - 1L]
    }
    tryCatch({
      a <- stats::nlminb(p, fn, gr, control = list(iter.max = control$iter.max,
        eval.max = 2 * control$iter.max, rel.tol = 1e-10))
      for (step in seq_len(8L)) {
        gradient <- gr(a$par)
        if (max(abs(gradient)) <= 1e-6) { a$convergence <- 0L; break }
        R <- tryCatch(chol(he(a$par)), error = function(e) NULL)
        if (is.null(R)) break
        direction <- as.numeric(backsolve(R, forwardsolve(t(R), gradient)))
        accepted <- FALSE
        for (half in 0:12) {
          candidate <- a$par - direction * 2^-half
          value <- fn(candidate)
          if (is.finite(value) && value <= a$objective + 1e-12 * max(1, abs(a$objective)) &&
              max(abs(gr(candidate))) < max(abs(gradient))) {
            a$par <- candidate; a$objective <- value; accepted <- TRUE; break
          }
        }
        if (!accepted) break
      }
      a$gradient_max <- max(abs(gr(a$par)))
      if (a$gradient_max <= 1e-6) a$convergence <- 0L
      a$usable <- a$convergence == 0L && is.finite(a$gradient_max) && a$gradient_max <= 1e-4
      a
    }, error = function(e) list(usable = FALSE, initial = p, objective = Inf,
                                message = conditionMessage(e)))
  })
  good <- which(vapply(attempts, function(z) isTRUE(z$usable), logical(1)))
  if (!length(good)) .msr_bayes_laplace_failure("No stationary conditional Laplace mode.", attempts)
  best <- good[which.min(vapply(attempts[good], function(z) z$objective, 0))]
  mode <- attempts[[best]]$par
  H <- he(mode); H <- (H + t(H)) / 2
  geometry <- tryCatch(.msr_bayes_laplace_geometry(H),
    error = function(e) .msr_bayes_laplace_failure(conditionMessage(e), attempts))
  distance <- vapply(attempts[good], function(z) {
    delta <- z$par - mode
    sqrt(max(0, drop(crossprod(delta, H %*% delta))))
  }, 0)
  competing <- any(vapply(attempts[good], function(z)
    z$objective - attempts[[best]]$objective <= 2, logical(1)) & distance > 1)
  # f includes the half-normal prior density in log(tau) coordinates. Remove
  # it because integration against the prior CDF already supplies that prior.
  log_prior_h <- stats::dnorm(exp(log_tau), 0, prior_sd, log = TRUE) + log(2) + log_tau
  list(mode = unpack(mode), root = geometry$root, covariance = geometry$covariance,
    hessian = H, log_density = -attempts[[best]]$objective +
      length(nuisance) / 2 * log(2 * pi) - sum(log(diag(chol(H)))) - log_prior_h,
    gradient_max = attempts[[best]]$gradient_max, condition_number = geometry$condition_number,
    competing = competing, attempts = attempts, selected_start = best)
}

.msr_bayes_marginal_grid <- function(evaluate, prior_sd, control) {
  cache <- new.env(parent = emptyenv())
  get_node <- function(x) {
    key <- sprintf("%.15f", x)
    if (exists(key, cache, inherits = FALSE)) return(get(key, cache))
    keys <- ls(cache)
    warm <- if (length(keys)) get(keys[which.min(abs(as.numeric(keys) - x))], cache) else NULL
    z <- tryCatch(evaluate(x, warm), error = function(e) {
      stop(structure(list(message = paste("Conditional quadrature failure at prior CDF", x,
        ":", conditionMessage(e)), call = NULL, x = x, starts = e$starts,
        completed_nodes = as.list(cache)), class = c("msr_laplace_error", "error", "condition")))
    })
    assign(key, z, cache); z
  }
  # Start near the prior median, then approach both boundaries with warm starts.
  get_node(.5)
  previous <- NULL; history <- list(); converged <- FALSE
  for (n in c(17L, 33L, 65L, 129L, 257L)) {
    if (n > control$max_grid) break
    x <- sort(unique(c(seq(0, 1, length.out = n), .99, .999, .9999, .99999)))
    nodes <- lapply(x, get_node)
    grid <- .msr_bayes_marginal_weights(x, vapply(nodes, function(z) z$log_density, 0))
    grid$nodes <- nodes; grid$prior_sd <- prior_sd
    quantiles <- .msr_bayes_marginal_quantile(grid, c(.025, .5, .975))
    if (is.null(previous)) {
      cdf_error <- quantile_error <- logz_error <- Inf
    } else {
      cdf_error <- max(abs(.msr_bayes_marginal_cdf(grid, x) - .msr_bayes_marginal_cdf(previous, x)))
      quantile_error <- max(abs(quantiles - .msr_bayes_marginal_quantile(previous, c(.025, .5, .975)))) / prior_sd
      logz_error <- abs(grid$log_normalizer - previous$log_normalizer)
      converged <- cdf_error <= control$cdf_tol && quantile_error <= control$quantile_tol && logz_error <= control$logz_tol
    }
    history[[length(history) + 1L]] <- data.frame(nodes = length(x), cdf_error = cdf_error,
      quantile_error = quantile_error, logz_error = logz_error)
    previous <- grid
    if (converged) break
  }
  grid$converged <- converged
  grid$refinement <- do.call(rbind, history)
  grid
}

.msr_bayes_marginal_weights <- function(x, log_density) {
  if (any(!is.finite(log_density))) stop("Nonfinite conditional Laplace marginal density.", call. = FALSE)
  density <- exp(log_density - max(log_density))
  dx <- diff(x)
  area <- dx * (utils::head(density, -1L) + utils::tail(density, -1L)) / 2
  z <- sum(area)
  basis <- (c(0, dx) + c(dx, 0)) / 2
  list(x = x, density = density / z, cumulative = c(0, cumsum(area) / z),
       weights = basis * density / z, log_normalizer = max(log_density) + log(z))
}

.msr_bayes_marginal_cdf <- function(grid, x) {
  i <- findInterval(x, grid$x, all.inside = TRUE)
  dx <- grid$x[i + 1L] - grid$x[i]
  t <- pmin(dx, pmax(0, x - grid$x[i]))
  grid$cumulative[i] + grid$density[i] * t +
    (grid$density[i + 1L] - grid$density[i]) * t^2 / (2 * dx)
}

.msr_bayes_marginal_quantile <- function(grid, probs) {
  x <- vapply(probs, function(p) stats::uniroot(function(t)
    .msr_bayes_marginal_cdf(grid, t) - p, c(0, 1), tol = 1e-12)$root, 0)
  grid$prior_sd * stats::qnorm((1 + x) / 2)
}

.msr_bayes_marginal_draws <- function(grid, hyper, draws, seed) {
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  ids <- sample.int(length(grid$x), draws, replace = TRUE, prob = grid$weights)
  ans <- matrix(NA_real_, draws, length(grid$nodes[[1L]]$mode))
  nuisance <- setdiff(seq_len(ncol(ans)), hyper)
  for (i in unique(ids)) {
    rows <- which(ids == i)
    left <- if (i > 1L) grid$x[i] - grid$x[i - 1L] else 0
    right <- if (i < length(grid$x)) grid$x[i + 1L] - grid$x[i] else 0
    go_left <- stats::runif(length(rows)) < left / (left + right)
    fraction <- sqrt(stats::runif(length(rows)))
    # Continuous triangular basis around the selected quadrature node.
    x <- ifelse(go_left, grid$x[i] - left + left * fraction,
                         grid$x[i] + right - right * fraction)
    x <- pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps, x))
    ans[rows, hyper] <- log(grid$prior_sd * stats::qnorm((1 + x) / 2))
    node <- grid$nodes[[i]]
    ans[rows, nuisance] <- .msr_bayes_laplace_draws(node$mode[nuisance], node$root, length(rows), NULL)
  }
  ans
}
