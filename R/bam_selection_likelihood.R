# Recover the normalized count likelihood and its actual nuisance count for
# native mgcv negative-binomial BAMs. Some provider versions retain the outer
# scale-estimated flag after fitting the family at unit scale.
.msr_bam_nb_selection <- function(model) {
  if (!inherits(model, "bam") || !inherits(model$family, "extended.family") ||
      !is.character(model$family$family) || length(model$family$family) != 1L ||
      !grepl("^negative binomial", model$family$family, ignore.case = TRUE)) return(NULL)
  if (!requireNamespace("mgcv", quietly = TRUE))
    stop("Native BAM likelihood verification requires mgcv.", call. = FALSE)
  method <- NULL
  for (kind in class(model)) {
    method <- utils::getS3method("logLik", kind, optional = TRUE)
    if (!is.null(method)) break
  }
  if (!identical(method, get("logLik.gam", asNamespace("mgcv")))) return(NULL)
  fail <- function() stop("Native negative-binomial BAM likelihood metadata do not match its normalized response likelihood and parameter count.", call. = FALSE)
  family <- model$family
  n_theta <- family$n.theta
  if (!is.numeric(n_theta) || length(n_theta) != 1L || is.na(n_theta) ||
      !n_theta %in% c(0, 1) || !is.function(family$getTheta)) fail()
  theta <- family$getTheta(TRUE)
  if (!is.numeric(theta) || length(theta) != 1L || !is.finite(theta) || theta <= 0) fail()
  native_family <- mgcv::nb(theta = if (n_theta == 0) theta else NULL)
  for (name in c("aic", "variance", "getTheta", "Dd")) {
    if (!is.function(family[[name]]) ||
        !identical(body(family[[name]]), body(native_family[[name]])) ||
        !identical(formals(family[[name]]), formals(native_family[[name]])))
      stop("Native negative-binomial BAM family implementation changed; its likelihood cannot be verified.", call. = FALSE)
  }
  phi <- model[["sig2", exact = TRUE]]
  if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi) || phi != 1) fail()
  y <- model[["y", exact = TRUE]]; mu <- model[["fitted.values", exact = TRUE]]
  weights <- model[["prior.weights", exact = TRUE]]
  scalar_column <- function(x) is.numeric(x) && is.null(dim(x)) && all(is.finite(x))
  if (!scalar_column(y) || !scalar_column(mu) || !scalar_column(weights) ||
      !length(y) || length(mu) != length(y) || length(weights) != length(y) ||
      any(y < 0) || any(y != floor(y)) || any(mu <= 0) || any(weights < 0) || !any(weights > 0)) fail()
  used <- weights > 0
  density <- stats::dnbinom(y[used], size = theta, mu = mu[used], log = TRUE)
  direct <- sum(weights[used] * density)
  native <- stats::logLik(model); native_k <- attr(native, "df")
  edf <- model[["edf2", exact = TRUE]]
  if (is.null(edf)) edf <- model[["edf", exact = TRUE]]
  if (!scalar_column(edf) || length(edf) != length(stats::coef(model)) ||
      !is.numeric(native_k) || length(native_k) != 1L || !is.finite(native_k) ||
      length(native) != 1L || !is.finite(native) || !is.finite(direct)) fail()
  count <- min(sum(edf), length(stats::coef(model))) + n_theta
  if (!is.finite(count) || count < 0) fail()
  close <- function(x, target) isTRUE(all.equal(as.numeric(x), target, tolerance = 1e-8, scale = 1))
  delta_likelihood <- as.numeric(native) - direct
  delta_count <- as.numeric(native_k) - count
  if (close(delta_likelihood, 0) && close(delta_count, 0))
    return(list(loglik = as.numeric(native), K = as.numeric(native_k), adjustment = NULL))
  scale_flag <- isTRUE(model[["scale.estimated", exact = TRUE]])
  roundoff_verified <- FALSE
  if (!scale_flag || !close(delta_likelihood, 1) || !close(delta_count, 1)) {
    # Subtracting large weighted likelihoods can erase a one-unit offset.
    # Keep the original comparison tolerance, and verify the provider's AIC
    # arithmetic separately from the stable normalized density calculation.
    if (!((scale_flag && close(delta_count, 1)) ||
          (!scale_flag && close(delta_count, 0))) || max(weights) <= 1) fail()
    raw_edf <- model[["edf", exact = TRUE]]
    native_aic <- model[["aic", exact = TRUE]]
    if (!scalar_column(raw_edf) || length(raw_edf) != length(stats::coef(model)) ||
        !is.numeric(native_aic) || length(native_aic) != 1L || !is.finite(native_aic)) fail()
    # Evaluate verified native code with explicit theta and base arithmetic.
    # Separate family closures can otherwise retain inconsistent theta state.
    theta_raw <- family$getTheta(FALSE)
    if (!is.numeric(theta_raw) || length(theta_raw) != 1L || !is.finite(theta_raw) ||
        !identical(as.numeric(exp(theta_raw)), as.numeric(theta))) fail()
    canonical_aic <- native_family$aic
    environment(canonical_aic) <- baseenv()
    response_aic <- canonical_aic(y = y, mu = mu, theta = theta_raw, wt = weights, dev = 0)
    expected_aic <- response_aic + 2 * sum(raw_edf)
    expected_loglik <- sum(raw_edf) + as.numeric(scale_flag) - expected_aic / 2
    rescaled <- weights[used] / max(weights)
    response_scaled <- -canonical_aic(y = y[used], mu = mu[used], theta = theta_raw, wt = rescaled, dev = 0) / 2
    if (!is.finite(expected_aic) || !is.finite(response_scaled) ||
        !close(native_aic, expected_aic) || !close(native, expected_loglik) ||
        !close(response_scaled, sum(rescaled * density))) fail()
    if (!scale_flag)
      return(list(loglik = as.numeric(native), K = as.numeric(native_k), adjustment = NULL))
    roundoff_verified <- TRUE
  }
  adjustment <- list(method = "verified_native_nb_unit_scale", native_logLik = as.numeric(native),
    logLik = direct, native_K = as.numeric(native_k), K = count,
    estimated_theta_parameters = n_theta, estimated_scale_parameters = 0L)
  if (roundoff_verified) adjustment$verification <- "native_aic_and_rescaled_density"
  list(loglik = direct, K = count, adjustment = adjustment)
}

.msr_selection_adjustments <- function(context) {
  out <- lapply(seq_along(context), function(i) {
    value <- if (is.null(context[[i]])) NULL else context[[i]]$adjustment
    if (is.null(value)) return(NULL)
    c(list(input_position = as.integer(i)), value)
  })
  names(out) <- names(context)
  out[!vapply(out, is.null, logical(1))]
}
