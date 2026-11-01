# Raw SQ buffers are retained per observation, including repeated prepared IDs.
# The pooled reference transform is owned by preparation, never by this evaluator.
.msr_scale_point_ids <- function(inputs) {
  if (!is.null(inputs$binned)) {
    ids <- inputs$binned$point_ids
    n <- nrow(inputs$binned$vsum[[1L]])
  } else {
    ids <- rownames(inputs$kernel_dat); n <- length(inputs$d_list)
  }
  if (length(n) != 1L || is.null(ids) || length(ids) != n || anyNA(ids) ||
      anyDuplicated(as.character(ids)) || any(!nzchar(as.character(ids))))
    stop("Prepared point_ids must uniquely identify matching spatial rows.", call. = FALSE)
  as.character(ids)
}

.msr_scale_surface_names <- function(inputs) {
  s <- inputs$scale_vars
  if (is.null(s)) character() else s$covariate[
    s$type == "surface" & s$metric == "sq" & s$weighted & s$optimize]
}

# Shift exponents by the nearest retained cell, and center values before squaring.
# Saturation is below double-precision weight support, as in the limiting engine.
.msr_scale_surface_raw <- function(buffers, log_scale, kernel, ad = FALSE) {
  power <- if (kernel == "gaussian") 2 else 1
  cap <- if (ad) RTMB::Vectorize(function(x) {
    "if" <- RTMB::ADoverload("if")
    if (x > log(746)) RTMB::advector(log(746)) else x
  }) else function(x) pmin(x, log(746))
  values <- lapply(seq_along(buffers), function(i) {
    b <- buffers[[i]]
    if (any(!is.finite(b$d) | b$d < 0))
      stop("Weighted-SQ cell distances must be finite and nonnegative.", call. = FALSE)
    if (all(b$v == b$v[1L])) return(0 * log_scale[i])
    delta <- pmax(b$d^power - min(b$d)^power, 0)
    w <- exp(-exp(cap(log(delta) - power * log_scale[i] - if (power == 2) log(2) else 0)))
    v <- b$v - b$v[1L]
    mu <- sum(w * v) / sum(w)
    variance <- sum(w * (v - mu)^2) / sum(w)
    if (!ad) return(sqrt(variance))
    "if" <- RTMB::ADoverload("if")
    # Avoid differentiating sqrt(0) on an unselected branch. No positive
    # variance is floored and no epsilon is added to the likelihood.
    safe <- if (variance > 0) variance else RTMB::advector(1)
    if (variance > 0) sqrt(safe) else 0 * log_scale[i]
  })
  do.call(c, values)
}

.msr_scale_surface_size <- function(surface) {
  cells <- sum(vapply(surface$buffers, function(x) sum(vapply(x, function(b) length(b$v), 0L)), 0))
  list(observation_cells = cells, retained_bytes = as.numeric(utils::object.size(surface)),
    # A planning estimate, not a bound on the differentiation library's tape.
    estimated_AD_bytes = cells * 512)
}

.msr_scale_surface_ram <- function(spec) {
  design <- spec$response_design
  if (!length(spec$surface$buffers) && is.null(design)) return(invisible(NULL))
  memory <- .msr_scale_surface_size(spec$surface)
  total <- .msr_total_ram_bytes()
  estimate <- 4 * memory$retained_bytes + memory$estimated_AD_bytes + 1024^3
  if (!is.null(design)) {
    if (!identical(dim(design$multipliers), c(design$n, design$p)) ||
        design$n != length(spec$y) || design$p != length(spec$parameters$beta))
      stop("Expanded response dimensions do not match the differentiation payload.", call. = FALSE)
    # Additional planning allowance for retained numeric multipliers and their
    # AD response operations. This is not a bound on the library's peak tape.
    estimate <- estimate + 4 * as.numeric(utils::object.size(design)) +
      512 * as.double(design$n) * design$p
  }
  if (!is.finite(total) || total <= 0) {
    if (is.null(design)) stop("Cannot detect RAM for weighted-SQ differentiation; no tape was launched.", call. = FALSE)
    stop("Cannot detect RAM for expanded-response differentiation; no tape was launched.", call. = FALSE)
  }
  if (estimate > .75 * total) {
    if (is.null(design)) stop("Weighted-SQ raw buffers exceed the conservative differentiation RAM budget.", call. = FALSE)
    stop("Expanded response and raw buffers exceed the conservative differentiation RAM budget.", call. = FALSE)
  }
  invisible(estimate)
}


# Exact equal-cell limits of weighted SQ. The nearest limit includes every
# exactly tied valid minimum-distance cell; no finite bandwidth approximates it.
.msr_scale_surface_limit <- function(buffers, limit) {
  limit <- match.arg(limit, c("broad", "nearest"))
  if (!is.list(buffers) || !length(buffers))
    stop("Exact SQ limits require retained nonempty raw buffers.", call. = FALSE)
  vapply(buffers, function(b) {
    if (!is.numeric(b$v) || !is.numeric(b$d) || !length(b$v) ||
        length(b$v) != length(b$d) || any(!is.finite(b$v)) ||
        any(!is.finite(b$d) | b$d < 0))
      stop("Exact SQ limits require matching finite values and nonnegative distances.", call. = FALSE)
    v <- if (limit == "nearest") b$v[b$d == min(b$d)] else b$v
    v <- v - v[1L]
    result <- sqrt(mean((v - mean(v))^2))
    if (!is.finite(result)) stop("Exact SQ limit is not finite.", call. = FALSE)
    result
  }, 0)
}
