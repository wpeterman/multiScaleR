# One-case conditional likelihood. Spatial preparation and choice-set validation
# are performed once, outside the joint coefficient/scale objective.
.msr_clogit_groups <- function(group, y = NULL) {
  if (!length(group) || anyNA(group)) stop("Choice-set IDs cannot be missing.", call. = FALSE)
  ids <- match(as.character(group), unique(as.character(group)))
  groups <- split(seq_along(ids), ids)
  if (any(lengths(groups) < 2L))
    stop("Every choice set requires at least two alternatives.", call. = FALSE)
  if (!is.null(y) && (length(y) != length(group) || anyNA(y) ||
      any(!y %in% c(0, 1)) || any(vapply(groups, function(i) sum(y[i]), 0) != 1)))
    stop("Every choice set must contain exactly one used step (0/1 response).", call. = FALSE)
  groups
}

.msr_clogit_contract <- function(model, inputs, implementation = "reference") {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("conditional_logit requires the optional survival package.", call. = FALSE)
  m <- .analysis_model(model)
  if (!inherits(m, "coxph") || !is.data.frame(m$model))
    stop("conditional_logit requires clogit/coxph or amt fit_clogit with model = TRUE.", call. = FALSE)
  if (!inputs$kernel %in% c("gaussian", "exp"))
    stop("conditional_logit requires Gaussian or exponential kernels.", call. = FALSE)
  if (!is.null(m$na.action) || !is.null(m$call$subset) || !is.null(m$call$weights) ||
      !is.null(m$weights) || !is.null(m$naive.var) || !is.null(m$call$cluster) ||
      !is.null(m$call$offset))
    stop("Missing-row omission, subsets, weights, robust clustering and separate offsets are unsupported; prepare complete choice sets explicitly.", call. = FALSE)
  tt <- stats::terms(m)
  labels <- attr(tt, "term.labels")
  sg <- grep("^strata\\(", labels)
  specials <- attr(tt, "specials")
  if (length(sg) != 1L || length(str2lang(labels[sg])) != 2L ||
      !is.symbol(str2lang(labels[sg])[[2L]]) ||
      inherits(m, "coxph.penal") || any(lengths(specials[setdiff(names(specials), "strata")]) > 0L) ||
      any(grepl("cluster\\(|frailty\\(|tt\\(|pspline\\(|strata\\(.*:", labels)))
    stop("Use one bare strata(choice_id) term; random effects and special Cox terms are unsupported.", call. = FALSE)
  strata_name <- as.character(str2lang(labels[sg])[[2L]])
  y <- stats::model.response(m$model)
  if (!inherits(y, "Surv") || ncol(y) != 2L || !identical(attr(y, "type"), "right") ||
      any(!is.finite(y)) || any(y[, 1L] != 1))
    stop("Only the constant-time Surv(... = 1, case) conditional-logistic response is supported.", call. = FALSE)
  groups <- .msr_clogit_groups(m$model[[labels[sg]]], y[, 2L])
  if (!m$method %in% c("exact", "efron", "breslow"))
    stop("Unsupported conditional-logistic likelihood method.", call. = FALSE)
  specs <- inputs$scale_vars
  surface_names <- if (is.null(specs)) character() else specs$covariate[
    specs$type == "surface" & specs$metric == "sq" & specs$weighted & specs$optimize]
  covs <- intersect(all.vars(stats::delete.response(tt)),
                    c(names(inputs$binned$vsum), surface_names))
  if (!length(covs) || strata_name %in% covs)
    stop("Require additive kernel predictors distinct from choice-set IDs.", call. = FALSE)
  for (label in labels) {
    ex <- str2lang(label); involved <- intersect(all.vars(ex), covs)
    if (length(involved) && !(length(involved) == 1L && identical(ex, as.name(involved))))
      stop("Kernel predictors must be bare additive terms; kernel interactions/transforms are unsupported.", call. = FALSE)
  }
  for (idx in attr(tt, "offset"))
    if (length(intersect(all.vars(attr(tt, "variables")[[idx + 1L]]), covs)))
      stop("Kernel predictors cannot occur in offsets.", call. = FALSE)
  if (!is.null(specs)) {
    used <- specs[specs$covariate %in% all.vars(tt), , drop = FALSE]
    optimized <- used[used$optimize, , drop = FALSE]
    if (!all(optimized$covariate %in% covs))
      stop("conditional_logit optimizes smooth kernel means and weighted SQ only; hard-radius surface optimization requires the refit backend.", call. = FALSE)
  }
  dat <- .model_data(m)
  required <- all.vars(stats::formula(m))
  if (!is.data.frame(dat) ||
      anyDuplicated(rownames(dat)) || anyNA(match(rownames(m$model), rownames(dat))))
    stop("Cannot recover explicit formula variables and retained row IDs.", call. = FALSE)
  dat <- dat[match(rownames(m$model), rownames(dat)), , drop = FALSE]
  response <- stats::formula(m)[[2L]]
  response_vars <- all.vars(response)
  if (!is.call(response) || length(response_vars) != 1L || !is.symbol(response[[length(response)]]))
    stop("Use a direct case response rather than a transformed response.", call. = FALSE)
  # insight may expose this validated constant-time response under `x`.
  if (!response_vars %in% names(dat)) dat[[response_vars]] <- as.numeric(y[, 2L])
  if (!all(required %in% names(dat)) || !identical(as.numeric(dat[[response_vars]]), as.numeric(y[, 2L])))
    stop("Cannot recover explicit formula variables and the original case response.", call. = FALSE)
  if (!identical(.msr_clogit_groups(dat[[strata_name]], y[, 2L]), groups))
    stop("Recovered data do not preserve the fitted choice sets.", call. = FALSE)
  X <- stats::model.matrix(m)
  if (anyNA(stats::coef(m)) || any(!is.finite(X)) || anyDuplicated(colnames(X)) ||
      !all(covs %in% colnames(X)))
    stop("Require finite, nonaliased coefficients and additive numeric kernel columns.", call. = FALSE)
  fixed <- setdiff(colnames(X), covs)
  offset <- stats::model.offset(m$model)
  if (is.null(offset)) offset <- rep(0, nrow(X))
  if (any(!is.finite(offset))) stop("Offsets must be finite.", call. = FALSE)
  rows <- .msr_profile_input_rows(m$model, inputs)
  engine <- if (any(covs %in% surface_names))
    .msr_prepare_surface_engine(inputs, covs, rows, implementation = implementation) else
    .msr_prepare_kernel_engine(inputs$binned, covs, rows, inputs$kernel,
                               implementation = implementation)
  # Rebuild prediction terms on the retained data: drop.terms() can leave Cox
  # predvars pointing at a removed strata term instead of a formula offset.
  offsets <- vapply(attr(tt, "offset"), function(i)
    paste(deparse(attr(tt, "variables")[[i + 1L]]), collapse = ""), "")
  pred_formula <- stats::reformulate(c(labels[-sg], offsets),
    intercept = attr(tt, "intercept") == 1L, env = environment(tt))
  pred_terms <- stats::terms(stats::model.frame(pred_formula, dat, na.action = stats::na.fail))
  list(model = m, data = dat, frame = m$model, covariates = covs,
    strata_name = strata_name, groups = groups, y = as.numeric(y[, 2L]),
    X_fixed = X[, fixed, drop = FALSE], offset = offset, engine = engine,
    terms = pred_terms, fixed_columns = fixed,
    contrasts = m$contrasts, xlevels = m$xlevels[setdiff(names(m$xlevels), labels[sg])],
    original_columns = colnames(X))
}

# Stable normalization also defines prediction probabilities on a finite set.
.msr_clogit_softmax <- function(eta, groups) {
  probability <- log_probability <- numeric(length(eta))
  for (s in seq_along(groups)) {
    i <- groups[[s]]; shift <- max(eta[i]); w <- exp(eta[i] - shift)
    log_probability[i] <- eta[i] - shift - log(sum(w)); probability[i] <- w / sum(w)
  }
  list(probability = probability, log_probability = log_probability)
}

.msr_clogit_objective <- function(contract) {
  e <- contract$engine; nk <- length(e$parameter_names)
  nb <- ncol(contract$X_fixed) + e$p
  cached <- NULL; previous <- NULL
  evaluate <- function(par) {
    if (identical(par, previous)) return(cached)
    if (length(par) != nk + nb || any(!is.finite(par))) stop("Invalid joint parameters.")
    kd <- .msr_eval_kernel_engine(e, unname(par[seq_len(nk)]))
    X <- cbind(contract$X_fixed, kd$value); beta <- par[nk + seq_len(nb)]
    eta <- drop(X %*% beta) + contract$offset
    if (any(!is.finite(eta))) stop("Nonfinite linear predictor.")
    z <- .msr_clogit_softmax(eta, contract$groups)
    residual <- z$probability - contract$y
    gradient <- c(colSums(kd$d_scale * residual) * utils::tail(beta, e$p),
                  drop(crossprod(X, residual)))
    value <- -sum(z$log_probability[contract$y == 1])
    previous <<- par
    cached <<- list(value = value, gradient = gradient, X = X,
      probability = z$probability, beta = stats::setNames(beta, colnames(X)))
    cached
  }
  list(evaluate = evaluate, fn = function(p) evaluate(p)$value,
       gr = function(p) evaluate(p)$gradient, nk = nk, nb = nb)
}

# Always initialize cold and independently validate terminal points with the
# public survival fitter. No internal survival routines or Hessian bypasses.
.msr_clogit_native <- function(contract, kp) {
  nd <- contract$data; kd <- .msr_eval_kernel_engine(contract$engine, unname(kp))
  for (nm in contract$covariates) nd[[nm]] <- kd$value[, nm]
  for (nm in names(contract$contrasts))
    if (nm %in% names(nd) && is.factor(nd[[nm]])) stats::contrasts(nd[[nm]]) <- contract$contrasts[[nm]]
  fit <- survival::coxph(stats::formula(contract$model), data = nd,
    method = contract$model$method, model = TRUE, x = TRUE, y = TRUE,
    na.action = stats::na.fail,
    control = survival::coxph.control(iter.max = 100L, eps = 1e-10))
  fit$call$formula <- stats::formula(contract$model)
  fit$call$method <- contract$model$method
  fit
}

#' @exportS3Method stats::update
#' @noRd
update.msr_conditional_clogit <- function(object, ..., data = NULL) {
  if (is.null(data)) stop("Supply explicit data with already standardized kernel covariates, or refit the original model and rerun multiScale_optim().", call. = FALSE)
  stats::update.default(object, ..., data = data)
}

.msr_fit_clogit_backend <- function(model, inputs, join_by, par, start_strategy,
                                    n_cores, refit_fn, control, verbose, call) {
  if (!is.null(join_by) || !is.null(refit_fn) || start_strategy != "single" ||
      (!is.null(n_cores) && n_cores != 1L))
    stop("conditional_logit uses retained rows and serial independent starts; join_by/refit_fn are unsupported.", call. = FALSE)
  allowed <- c("starts", "n_starts", "iter.max", "gradient_tolerance", "kernel_engine")
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      any(!names(control) %in% allowed) || anyDuplicated(names(control)))))
    stop("conditional_logit controls: starts, n_starts, iter.max, gradient_tolerance, kernel_engine.", call. = FALSE)
  ctl <- .msr_profile_control(control)
  contract <- .msr_clogit_contract(model, inputs, ctl$kernel_engine)
  e <- contract$engine; obj <- .msr_clogit_objective(contract)
  b <- .msr_profile_starts(list(objective = list(engine = e, family = "conditional_logit",
    parameter_names = e$parameter_names)), inputs, par, ctl, contract$model)
  if (any(sweep(b$starts, 2L, b$lower, "<") | sweep(b$starts, 2L, b$upper, ">")))
    stop("All starts must be inside the scale bounds.", call. = FALSE)
  lower <- c(b$lower, rep(-Inf, obj$nb)); upper <- c(b$upper, rep(Inf, obj$nb))
  records <- vector("list", nrow(b$starts)); fits <- vector("list", length(records))
  for (i in seq_along(records)) {
    t0 <- proc.time()[["elapsed"]]; warnings <- character()
    ans <- tryCatch(withCallingHandlers({
      kp <- b$starts[i, ]
      X0 <- obj$evaluate(c(kp, rep(0, obj$nb)))$X
      for (g in contract$groups) X0[g, ] <- sweep(X0[g, , drop = FALSE], 2L, colMeans(X0[g, , drop = FALSE]))
      coefficient_units <- sqrt(colMeans(X0^2))
      if (any(!is.finite(coefficient_units) | coefficient_units <= 0))
        stop("Predictors must vary within choice sets.")
      # Optimize dimensionless coefficients, avoiding meter-scale step lengths
      # dominating the search. Eligibility still uses the original-unit score.
      initial <- stats::optim(rep(0, obj$nb), function(beta) obj$fn(c(kp, beta / coefficient_units)),
        function(beta) utils::tail(obj$gr(c(kp, beta / coefficient_units)), obj$nb) / coefficient_units, method = "BFGS",
        control = list(maxit = ctl$iter.max, reltol = 1e-12))
      if (initial$convergence != 0L) stop("Cold coefficient initialization did not converge.")
      units <- c(rep(1, obj$nk), coefficient_units)
      safe <- function(p) tryCatch({
        z <- obj$evaluate(p / units); z$gradient <- z$gradient / units; z
      }, error = function(e)
        list(value = 1e100, gradient = rep(0, length(p))))
      fit <- stats::optim(c(kp, initial$par), function(p) safe(p)$value,
        function(p) safe(p)$gradient, method = "L-BFGS-B", lower = lower, upper = upper,
        control = list(maxit = ctl$iter.max, factr = 10, pgtol = 1e-9))
      fit$par <- fit$par / units
      exact <- obj$evaluate(fit$par)
      pg <- max(abs(.msr_projected_gradient(fit$par, exact$gradient, lower, upper)))
      native <- .msr_clogit_native(contract, fit$par[seq_len(obj$nk)])
      beta <- stats::coef(native)[names(exact$beta)]
      gap <- abs(exact$value + as.numeric(stats::logLik(native)))
      coefficient_gap <- max(abs(beta - exact$beta))
      V <- stats::vcov(native)
      centered <- exact$X
      for (g in contract$groups) centered[g, ] <- sweep(exact$X[g, , drop = FALSE], 2L,
        colSums(exact$X[g, , drop = FALSE] * exact$probability[g]))
      H <- crossprod(centered, centered * exact$probability)
      inverse <- .msr_information_inverse(H)
      covariance_gap <- if (is.null(inverse$vcov)) Inf else
        max(abs(inverse$vcov - V[colnames(exact$X), colnames(exact$X)]))
      # The returned model uses cold native coefficients. Refine the terminal
      # coefficient point only after auditing the original joint solution.
      # Preserve that score and all original agreement diagnostics for review.
      joint_pg <- pg
      cold_pg <- Inf
      refined <- FALSE
      if (fit$convergence == 0L && is.finite(gap) && gap <= 1e-5 &&
          coefficient_gap <= 1e-4 && is.finite(covariance_gap) &&
          covariance_gap <= 1e-4 && !length(warnings)) {
        cold_par <- c(fit$par[seq_len(obj$nk)], unname(beta))
        cold <- obj$evaluate(cold_par)
        cold_pg <- max(abs(.msr_projected_gradient(cold_par, cold$gradient, lower, upper)))
        if (cold$value <= exact$value && cold_pg < pg) {
          fit$par <- cold_par
          exact <- cold
          pg <- cold_pg
          refined <- TRUE
        }
      }
      eligible <- isTRUE(fit$convergence == 0L && pg <= ctl$gradient_tolerance &&
        is.finite(gap) && gap <= 1e-5 && coefficient_gap <= 1e-4 &&
        is.finite(covariance_gap) && covariance_gap <= 1e-4 && !length(warnings))
      list(par = fit$par, value = exact$value, convergence = fit$convergence,
        joint_gradient = joint_pg, cold_gradient = cold_pg, refined = refined,
        gradient = pg, gap = gap, coefficient_gap = coefficient_gap,
        covariance_gap = covariance_gap, eligible = eligible, native = native,
        message = if (is.null(fit$message)) "" else fit$message)
    }, warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
    error = function(e) list(par = NULL, value = NA_real_, convergence = NA_integer_,
      joint_gradient = Inf, cold_gradient = Inf, refined = FALSE,
      gradient = Inf, gap = NA_real_, coefficient_gap = NA_real_, covariance_gap = NA_real_,
      eligible = FALSE, message = conditionMessage(e)))
    fits[[i]] <- ans
    records[[i]] <- data.frame(start = i, value = ans$value, convergence = ans$convergence,
      projected_gradient = ans$gradient, native_gap = ans$gap, coefficient_gap = ans$coefficient_gap,
      joint_projected_gradient = ans$joint_gradient, cold_projected_gradient = ans$cold_gradient,
      terminal_coefficients_refined = ans$refined,
      covariance_gap = ans$covariance_gap, eligible = ans$eligible,
      seconds = proc.time()[["elapsed"]] - t0, message = ans$message,
      warnings = paste(unique(warnings), collapse = "; "))
  }
  records <- do.call(rbind, records); good <- which(records$eligible)
  if (!length(good)) stop(structure(list(message = "No conditional_logit start passed score and cold native audits.",
    call = NULL, start_results = records), class = c("msr_clogit_fit_error", "error", "condition")))
  winner <- good[which.min(records$value[good])]; best <- fits[[winner]]
  kp <- unname(best$par[seq_len(obj$nk)]); kd <- .msr_eval_kernel_engine(e, kp)
  class(best$native) <- c("msr_conditional_clogit", class(best$native))
  opt_mod <- best$native
  if (inherits(model, "fit_clogit")) {opt_mod <- model; opt_mod$model <- best$native}
  attr(opt_mod, "scale_uncertainty") <- "conditional on optimized kernel scales"
  scale <- exp(kp) * inputs$unit_conv
  # Store only the reusable data contract, not closures retaining all start fits.
  contract$objective <- list(engine = e, X_fixed = contract$X_fixed, offset = contract$offset)
  if (inherits(e, "msr_surface_engine")) {
    # Raw cells already live in kernel_inputs. Rebuild the filtered prediction
    # buffers on demand rather than retaining a second copy in every fit.
    contract$objective$engine$buffers <- NULL
    contract$objective$engine$means <- NULL
    contract$objective$engine$rebuild_from_inputs <- TRUE
  }
  # The original and selected native models remain in the public result.
  # Retain only prediction state here, avoiding repeated serialization of the
  # fitting engine, original model/data and design. Prediction needs row IDs,
  # not another complete model frame. Older full contracts remain readable.
  contract <- contract[setdiff(names(contract), c("model", "data", "engine", "X_fixed", "y"))]
  contract$frame <- contract$frame[, FALSE, drop = FALSE]
  H <- matrix(NA_real_, obj$nk, obj$nk)
  out <- list(backend = "conditional_logit", scale_est = data.frame(Mean = scale, SE = NA_real_, row.names = e$covariates),
    shape_est = NULL, opt_mod = opt_mod, fitted_mod_original = model,
    optim_results = list(par = exp(kp), par_unscale = scale, value = best$value,
      convergence = best$convergence, message = best$message, hessian = H, hessian_unscale = H),
    min_D = inputs$min_D, max_D = inputs$max_D,
    kernel_inputs = inputs[setdiff(names(inputs), c("min_D", "max_D"))],
    scl_params = list(mean = stats::setNames(kd$center, e$covariates), sd = stats::setNames(kd$spread, e$covariates)),
    profile_fit = list(contract = contract, par = kp, starts = b$starts, start_results = records,
      winner = winner, lower = b$lower, upper = b$upper, control = ctl),
    uncertainty = list(available = FALSE, reasons = "Step-selection scale uncertainty is not calibrated."),
    join_by = NULL, opt_context = NULL, profile_scale_est = NULL,
    diagnostics = .empty_diagnostics(), next_run = NULL, warn_message = 0, call = call)
  out$diagnostics$conditional_logit <- list(attempted_starts = nrow(records), eligible_starts = length(good),
    boundary = any(kp - b$lower < 1e-5 | b$upper - kp < 1e-5),
    uncertainty_available = FALSE, reasons = out$uncertainty$reasons)
  class(out) <- "multiScaleR"
  if (verbose) message("Conditional logit: ", length(good), "/", nrow(records), " starts passed native audits; scale SEs withheld.")
  out
}

.msr_predict_clogit <- function(object, newdata, kernel_inputs, type, se.fit, uncertainty, interval) {
  if (interval != "none" || (se.fit && uncertainty != "conditional"))
    stop("conditional_logit supports conditional SEs only; joint uncertainty and intervals are unavailable.", call. = FALSE)
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) stop("se.fit must be logical.")
  contract <- object$profile_fit$contract
  groups <- contract$groups
  if (!is.null(newdata) || !is.null(kernel_inputs)) {
    if (!is.data.frame(newdata) || !contract$strata_name %in% names(newdata))
      stop("New predictions require explicit choice-set IDs in newdata.", call. = FALSE)
    groups <- .msr_clogit_groups(newdata[[contract$strata_name]])
  }
  d <- .msr_profile_prediction_design(object, newdata, kernel_inputs)
  if (inherits(object, "multiScaleR_prediction")) {
    beta <- object$coefficients
    V <- object$conditional_vcov
  } else {
    m <- .analysis_model(object$opt_mod)
    beta <- stats::coef(m)[colnames(d$X)]
    V <- if (se.fit) stats::vcov(m)[colnames(d$X), colnames(d$X), drop = FALSE] else NULL
  }
  eta <- drop(d$X %*% beta) + d$offset
  probability <- .msr_clogit_softmax(eta, groups)$probability
  fit <- if (type == "link") eta else probability
  names(fit) <- d$rows
  if (!se.fit) return(fit)
  J <- d$X
  if (type == "response") for (g in groups)
    J[g, ] <- probability[g] * sweep(d$X[g, , drop = FALSE], 2L,
      colSums(d$X[g, , drop = FALSE] * probability[g]))
  list(fit = fit, se.fit = sqrt(pmax(0, rowSums((J %*% V) * J))),
       type = type, uncertainty = "conditional")
}
