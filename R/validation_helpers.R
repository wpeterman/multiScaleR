validate_scalar_logical <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", arg), call. = FALSE)
  }
}

validate_scalar_numeric <- function(x,
                                    arg,
                                    lower = -Inf,
                                    upper = Inf,
                                    inclusive_lower = TRUE,
                                    inclusive_upper = TRUE,
                                    integerish = FALSE,
                                    positive = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x)) {
    stop(sprintf("`%s` must be a single finite numeric value.", arg), call. = FALSE)
  }

  if (integerish && x != as.integer(x)) {
    stop(sprintf("`%s` must be a whole number.", arg), call. = FALSE)
  }

  if (positive && x <= 0) {
    stop(sprintf("`%s` must be > 0.", arg), call. = FALSE)
  }

  lower_ok <- if (inclusive_lower) x >= lower else x > lower
  upper_ok <- if (inclusive_upper) x <= upper else x < upper

  if (!lower_ok || !upper_ok) {
    if (is.finite(lower) && is.finite(upper)) {
      comparator <- paste0(
        if (inclusive_lower) "[" else "(",
        lower, ", ", upper,
        if (inclusive_upper) "]" else ")"
      )
      stop(sprintf("`%s` must be within %s.", arg, comparator), call. = FALSE)
    } else if (is.finite(lower)) {
      comparator <- if (inclusive_lower) paste0(">= ", lower) else paste0("> ", lower)
      stop(sprintf("`%s` must be %s.", arg, comparator), call. = FALSE)
    } else {
      comparator <- if (inclusive_upper) paste0("<= ", upper) else paste0("< ", upper)
      stop(sprintf("`%s` must be %s.", arg, comparator), call. = FALSE)
    }
  }
}

validate_numeric_vector <- function(x,
                                    arg,
                                    length_ = NULL,
                                    positive = FALSE) {
  if (!is.numeric(x) || length(x) == 0 || anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must be a numeric vector with finite values.", arg), call. = FALSE)
  }

  if (!is.null(length_) && length(x) != length_) {
    stop(sprintf("`%s` must have length %s.", arg, length_), call. = FALSE)
  }

  if (positive && any(x <= 0)) {
    stop(sprintf("`%s` must contain only values > 0.", arg), call. = FALSE)
  }
}

validate_character_scalar <- function(x, arg) {
  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    stop(sprintf("`%s` must be a non-empty character string.", arg), call. = FALSE)
  }
}

validate_reserved_variable_names <- function(variable_names, context) {
  if ("x" %in% variable_names) {
    warning(
      paste0(
        "`x` is reserved for the x-coordinate column used internally during ",
        "raster extraction and cannot be used as a ", context, ". Rename the variable and rerun the analysis."
      ),
      call. = FALSE
    )
    stop(
      "Analysis stopped because `x` would confound with the internal x/y coordinate columns.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# Refuse a Bayesian fit in a function that needs a multiScale_optim() fit, and
# name the Bayesian tools that cover the same need. `fn` is the calling
# function's name, used in the message.
.msr_stop_bayes_fit <- function(x, fn = NULL) {
  if (!inherits(x, "multiScaleR_bayes")) return(invisible(NULL))
  fn <- if (is.null(fn)) "This function" else paste0("`", fn, "()`")
  stop(sprintf(paste0(
    "%s works on fits from `multiScale_optim()`, not on a Bayesian ",
    "`multiScaleR_bayes` fit from `multiScale_bayes()`. For a Bayesian fit, ",
    "use `posterior_scales()` to summarize the posterior scales, ",
    "`predict_scales()` to predict scales for new covariates or groups, and ",
    "`predict()` with new kernel inputs for posterior predictions. For a map, ",
    "pass a posterior scale to `kernel_scale.raster(sigma = )`."), fn),
    call. = FALSE)
}

validate_multiScaleR_input <- function(x, arg = "multiScaleR", fn = NULL) {
  .msr_stop_bayes_fit(x, fn)
  if (!inherits(x, c("multiScaleR", "multiScaleR_data"))) {
    stop(sprintf("`%s` must be a `multiScaleR` or `multiScaleR_data` object.", arg), call. = FALSE)
  }
}

validate_covariates_before_scale <- function(x,
                                             context = "calculated covariates",
                                             scale_vars = NULL) {
  x <- as.matrix(x)
  if (ncol(x) == 0) {
    return(invisible(TRUE))
  }

  for (j in seq_len(ncol(x))) {
    covariate <- colnames(x)[[j]]
    if (is.null(covariate) || is.na(covariate) || !nzchar(covariate)) {
      covariate <- paste0("column ", j)
    }

    values <- x[, j]
    non_missing <- !is.na(values)
    finite <- is.finite(values)

    if (!any(non_missing)) {
      stop(
        paste0(
          sprintf("Cannot scale %s: covariate `%s` is all NA.", context, covariate),
          .msr_covariate_scale_hint(covariate, scale_vars)
        ),
        call. = FALSE
      )
    }
    if (any(non_missing & !finite)) {
      stop(
        sprintf("Cannot scale %s: covariate `%s` contains non-finite values.", context, covariate),
        call. = FALSE
      )
    }
    if (!any(finite)) {
      stop(
        sprintf("Cannot scale %s: covariate `%s` has no finite values.", context, covariate),
        call. = FALSE
      )
    }

    finite_values <- values[finite]
    if (length(unique(finite_values)) < 2 || isTRUE(stats::sd(finite_values) == 0)) {
      stop(
        paste0(
          sprintf("Cannot scale %s: covariate `%s` has zero variance.", context, covariate),
          .msr_start_scale_hint(covariate, scale_vars, context),
          .msr_covariate_scale_hint(covariate, scale_vars)
        ),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

# kernel_prep() evaluates every estimated landscape or surface metric once at
# its starting scale. A metric that is constant there (for example `pr` when
# every buffer holds the same classes) may still vary at other scales, so name
# the starting scale and the `sigma =` argument that changes it. Other callers
# (optimization, fixed radii) keep the generic message.
.msr_start_scale_hint <- function(covariate, scale_vars = NULL, context = "") {
  start_contexts <- c("`kernel_prep()` covariates", "pooled group covariates")
  if (!context %in% start_contexts || is.null(scale_vars) ||
      !inherits(scale_vars, "multiScaleR_vars")) {
    return("")
  }
  spec <- scale_vars[scale_vars$covariate == covariate, , drop = FALSE]
  if (nrow(spec) != 1 || !spec$type %in% c("landscape", "surface") ||
      !isTRUE(spec$optimize)) {
    return("")
  }
  paste0(
    " Metric `", spec$metric, "` is constant across the sample points at the ",
    "starting scale (`max_D / 2` unless `sigma` is supplied), although it may ",
    "vary at other scales. Supply a starting scale at which it varies with ",
    "`sigma =` (in map units, one value per estimated scale)."
  )
}

.msr_covariate_scale_hint <- function(covariate, scale_vars = NULL) {
  if (is.null(scale_vars) || !inherits(scale_vars, "multiScaleR_vars")) {
    return("")
  }

  spec <- scale_vars[scale_vars$covariate == covariate, , drop = FALSE]
  if (nrow(spec) != 1 || !identical(spec$type, "landscape")) {
    return("")
  }

  metric <- spec$metric
  if (metric == "iji") {
    return(paste0(
      " Metric `iji` is fragile for scale optimization because it requires ",
      "at least three classes and enough between-class adjacencies in each ",
      "buffer; prescreen it for usable variation, increase the radius or ",
      "`max_D`, or use a more stable diversity or entropy metric."
    ))
  }

  if (metric %in% c("siei", "msiei")) {
    return(paste0(
      " Metric `", metric, "` can lose variation when local class evenness ",
      "is uniform or saturated; prescreen it across sample points before ",
      "optimizing its scale, increase the radius or `max_D`, or use the ",
      "corresponding diversity metric."
    ))
  }

  ""
}
