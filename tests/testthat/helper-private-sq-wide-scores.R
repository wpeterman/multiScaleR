# Same centered signed integration as the original SD0.3/2 software audit.
# The independent wide-scale normalizer supplies the conditional density.
wide_sq_score_moments <- function(g, p) {
  tape <- .msr_factorized_group_score_tape(g, p)
  normalizer <- wide_sq_integral(g, p)
  bound <- if (g$family == "gaussian")
    -length(g$y) * (p$log_dispersion + .5 * log(2 * pi)) else 0
  cache <- new.env(parent = emptyenv())
  node <- function(z, order = 0L) {
    key <- sprintf("%a", z)
    record <- if (exists(key, cache, inherits = FALSE)) get(key, cache, inherits = FALSE) else list()
    v <- tape$obj$par; v[tape$z] <- z
    if (is.null(record$value)) record$value <- as.numeric(tape$obj$fn(v))
    if (order >= 1L && is.null(record$score)) record$score <- as.numeric(tape$obj$gr(v))[tape$outer]
    if (order >= 2L && is.null(record$hessian)) record$hessian <- tape$obj$he(v)[tape$outer, tape$outer, drop = FALSE]
    assign(key, record, cache); record
  }
  weight <- function(z) {
    if (exp(bound + dnorm(z, log = TRUE) + normalizer) == 0) return(0)
    exp(-node(z)$value + normalizer)
  }
  integrate_moment <- function(fn, absolute) integrate(function(z) vapply(z, fn, 0),
    -Inf, Inf, rel.tol = 1e-9, abs.tol = absolute,
    subdivisions = 500L, stop.on.error = TRUE)$value
  mass <- integrate_moment(weight, 0)
  center <- node(0, 1L)$score
  difference <- vapply(seq_along(center), function(j) integrate_moment(function(z) {
    w <- weight(z); if (w == 0) return(0)
    w * (node(z, 1L)$score[j] - center[j])
  }, 1e-10), 0)
  score <- center + difference
  H <- matrix(0, length(center), length(center))
  for (j in seq_along(center)) for (k in seq_len(j)) {
    H[j, k] <- H[k, j] <- integrate_moment(function(z) {
      w <- weight(z); if (w == 0) return(0)
      a <- node(z, 2L); delta <- a$score - center
      w * (a$hessian[j, k] - delta[j] * delta[k])
    }, 1e-8) + difference[j] * difference[k]
  }
  # Return compact numerical summaries, never the runtime tape or node cache.
  list(normalizer = normalizer, mass = mass, score = score, hessian = H,
    parameter_names = names(tape$obj$par)[tape$outer])
}
