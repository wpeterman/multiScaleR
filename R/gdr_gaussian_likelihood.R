# Private fixed-point prerequisites. These functions do not admit fitted models.
.msr_gdr_annulus_program <- function() {
  if (!requireNamespace("RTMB", quietly = TRUE) ||
      as.character(getNamespaceVersion("RTMB")) != "2.0")
    stop("This private annulus derivative program requires qualified RTMB 2.0.")
  programs <- vector("list", 7L)
  for (order in 6:0) {
    # Only the numeric order and next moment are captured, never a caller/model.
    value_env <- list2env(list(n = order), parent = baseenv())
    value <- eval(quote(function(x) {
      if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0))
        stop("Annulus moments require finite nonnegative numeric arguments.")
      out <- numeric(length(x))
      small <- x <= .1
      if (any(small)) {
        z <- x[small]; term <- rep(1, length(z)); answer <- term/(n + 1)
        for (k in 1:18) {
          term <- term * (-z)/k
          answer <- answer + term/(n + k + 1)
        }
        out[small] <- answer
      }
      if (any(!small)) {
        z <- x[!small]
        out[!small] <- exp(lgamma(n + 1) + stats::pgamma(z, shape = n + 1, log.p = TRUE) -
          (n + 1) * log(z))
      }
      if (any(!is.finite(out)) || any(out < 0)) stop("Unrepresentable annulus moments.")
      out
    }), value_env)
    if (order == 6L) {
      derivative <- eval(quote(function(x, y, dy)
        stop("Annulus derivatives beyond the declared moment chain are unqualified.")), baseenv())
    } else {
      derivative_env <- list2env(list(higher = programs[[order + 2L]]), parent = baseenv())
      derivative <- eval(quote(function(x, y, dy) -higher(x) * dy), derivative_env)
    }
    # ADjoint retains promises: bind the current function values before looping.
    programs[[order + 1L]] <- do.call(RTMB::ADjoint, list(f = value, df = derivative))
  }
  list(evaluate = programs[[1L]], maximum_qualified_derivative = 4L,
    requires_serial_R_callbacks = TRUE, portable_objective = FALSE,
    fit_admission_established = FALSE)
}

# Validate an explicit native numeric layout, not a fitted-model approximation.
.msr_gdr_gaussian_data <- function(data, parameters) {
  fail <- function(message) stop(message, call. = FALSE)
  required <- c("Xlambda", "Xdist", "Xrem", "yd", "yr", "block", "observer",
    "periods", "breaks", "area", "K", "annulus_area_m2", "annulus_weights", "row_ids")
  if (!is.list(data) || is.null(names(data)) || anyNA(names(data)) ||
      anyDuplicated(names(data)) || !setequal(names(data), required))
    fail("Supply the complete named private GDR numeric layout.")
  real <- function(x) (is.double(x) || is.integer(x)) && !is.object(x) && all(is.finite(x))
  matrix_ok <- function(x) is.matrix(x) && real(x) && nrow(x) > 0L && ncol(x) > 0L
  if (!matrix_ok(data$yd) || !matrix_ok(data$yr) ||
      nrow(data$yd) != nrow(data$yr) ||
      any(data$yd < 0 | data$yd != floor(data$yd)) ||
      any(data$yr < 0 | data$yr != floor(data$yr)) ||
      any(rowSums(data$yd) != rowSums(data$yr)))
    fail("Use aligned complete integer distance/removal margins with equal totals.")
  M <- nrow(data$yd); J <- ncol(data$yr)
  if (!real(data$breaks) || !is.null(dim(data$breaks)) ||
      !identical(as.numeric(data$breaks), c(0, 25, 75)) || ncol(data$yd) != 2L)
    fail("This private native geometry requires complete c(0,25,75) meter point bins.")
  if (!real(data$K) || !is.null(dim(data$K)) || length(data$K) != 1L ||
      data$K < 0 || data$K != floor(data$K) || data$K > .Machine$integer.max ||
      any(rowSums(data$yd) > data$K))
    fail("A single finite integer cutoff must cover every observed total.")
  if (!real(data$periods) || !is.null(dim(data$periods)) || length(data$periods) != J ||
      any(data$periods <= 0 | data$periods != floor(data$periods)))
    fail("Removal durations must be positive integers aligned with the interval columns.")
  rows <- c(Xlambda = M, Xdist = M, Xrem = M * J)
  for (name in names(rows)) {
    X <- data[[name]]; columns <- colnames(X)
    if (!matrix_ok(X) || nrow(X) != rows[[name]] || is.null(columns) ||
        length(columns) != ncol(X) || anyNA(columns) || any(!nzchar(columns)) || anyDuplicated(columns))
      fail("Retain finite native design matrices, their original rows and unique column labels.")
  }
  ids <- data$row_ids
  if (!is.list(ids) || !identical(names(ids), c("survey", "removal_survey", "removal_interval")) ||
      !is.character(ids$survey) || length(ids$survey) != M || anyNA(ids$survey) ||
      any(!nzchar(ids$survey)) || anyDuplicated(ids$survey) ||
      !identical(ids$removal_survey, rep(ids$survey, each = J)) ||
      !real(ids$removal_interval) || !is.null(dim(ids$removal_interval)) ||
      !identical(as.numeric(ids$removal_interval), as.numeric(rep(seq_len(J), M))))
    fail("Retain unique survey identities and exact site-major removal interval ordering.")
  if (!real(data$area) || !is.null(dim(data$area)) || length(data$area) != M ||
      any(abs(data$area - pi * .075^2) >= 1e-15))
    fail("Expected circle abundance must use the original square-kilometer area.")
  annulus <- pi * diff(data$breaks^2)
  expected_area <- matrix(rep(annulus, each = M), M, 2L)
  expected_weights <- expected_area/(pi * 75^2)
  if (!matrix_ok(data$annulus_area_m2) || !matrix_ok(data$annulus_weights) ||
      !identical(dim(data$annulus_area_m2), c(M, 2L)) ||
      !identical(dim(data$annulus_weights), c(M, 2L)) ||
      any(abs(data$annulus_area_m2 - expected_area) >= 1e-10) ||
      any(abs(data$annulus_weights - expected_weights) >= 1e-12))
    fail("Native annulus areas and weights must describe the complete sampled circle.")
  blocks <- c("beta_lambda", "b_lambda", "lsigma_lambda", "beta_phi", "b_phi", "lsigma_phi",
    "beta_dist", "b_dist", "lsigma_dist", "beta_rem", "b_rem", "lsigma_rem", "beta_alpha", "beta_scale")
  if (!is.list(parameters) || is.null(names(parameters)) || anyNA(names(parameters)) ||
      anyDuplicated(names(parameters)) || !setequal(names(parameters), blocks) ||
      !all(vapply(parameters, function(x) real(x) && is.null(dim(x)), logical(1))))
    fail("Retain complete finite native parameter blocks, including empty inactive blocks.")
  if (length(parameters$beta_lambda) != ncol(data$Xlambda) ||
      length(parameters$beta_dist) != ncol(data$Xdist) ||
      length(parameters$beta_rem) != ncol(data$Xrem))
    fail("Regression coordinates must match every original native design column.")
  inactive <- c("beta_phi", "b_phi", "lsigma_phi", "b_rem", "lsigma_rem", "beta_alpha", "beta_scale")
  if (any(lengths(parameters[inactive])))
    fail("Annual availability, removal fields and non-Poisson/half-normal blocks are unqualified.")
  for (field in c("lambda", "dist")) {
    index <- data[[if (field == "lambda") "block" else "observer"]]
    states <- parameters[[paste0("b_", field)]]
    sd <- parameters[[paste0("lsigma_", field)]]
    if (!real(index) || !is.null(dim(index)) ||
        (length(states) == 0L && (length(index) || length(sd))) ||
        (length(states) > 0L && (length(index) != M || any(index < 1 | index != floor(index)) ||
          any(index > length(states)) || length(sd) != 1L ||
          !is.finite(exp(sd)) || exp(sd) <= 0)))
      fail("Each present iid intercept field needs its complete states, one log-SD and exact row map.")
  }
  invisible(data)
}

#' Private immutable observation context for normalized finite-K Gaussian GDR
#' @noRd
.msr_gdr_gaussian_context <- function(data, parameters, .lambda_loading = NULL) {
  .msr_gdr_gaussian_data(data, parameters)
  if (!is.null(.lambda_loading) && (!is.matrix(.lambda_loading) ||
      !is.numeric(.lambda_loading) || any(!is.finite(.lambda_loading)) ||
      !length(data$block) || !identical(dim(.lambda_loading),
        c(nrow(data$Xlambda), length(parameters$b_lambda)))))
    stop("Retain a finite abundance loading matrix in the exact survey/state layout.", call. = FALSE)
  abundance_dimension <- dim(data$Xlambda)
  abundance_columns <- colnames(data$Xlambda)
  abundance_rows <- rownames(data$Xlambda)
  # Keep invariant inputs separate from the abundance matrix that is refreshed.
  data <- data[setdiff(names(data), "Xlambda")]
  n <- rowSums(data$yd)
  cdf <- .msr_gdr_poisson_logcdf_program(data$K - n)$evaluate
  ratio <- .msr_gdr_annulus_program()$evaluate
  pd_positive <- which(data$yd > 0)
  pr_positive <- which(data$yr > 0)
  constant <- lgamma(n + 1) - rowSums(lgamma(data$yd + 1)) - rowSums(lgamma(data$yr + 1))
  M <- length(n)
  radius <- utils::tail(data$breaks, 1L)
  lower <- utils::head(data$breaks, -1L)/radius
  upper <- utils::tail(data$breaks, -1L)/radius
  width <- (upper - lower) * (upper + lower)
  block <- rep(names(parameters), lengths(parameters))
  within <- unlist(lapply(lengths(parameters), seq_len), use.names = FALSE)
  kind <- ifelse(startsWith(block, "b_"), "latent_state",
    ifelse(startsWith(block, "lsigma_"), "log_sd", "regression"))
  layout <- data.frame(position = seq_along(block), block = block, within_block = within, kind = kind)
  implementation <- lapply(list(.msr_gdr_gaussian_program, .msr_gdr_gaussian_context,
    .msr_gdr_gaussian_data, .msr_gdr_annulus_program,
    .msr_gdr_poisson_logcdf_program), function(fn)
      list(formals = deparse(formals(fn)), body = deparse(body(fn))))
  response_fingerprint <- digest::digest(list(data = data, lambda_loading = .lambda_loading,
    abundance_dimension = abundance_dimension, abundance_columns = abundance_columns,
    abundance_rows = abundance_rows, parameter_layout = layout, implementation = implementation,
    RTMB = as.character(getNamespaceVersion("RTMB"))), algo = "sha256")
  force(data); force(parameters); force(cdf); force(ratio); force(constant); force(.lambda_loading)
  prepare <- function(Xlambda, row_ids) {
    if (missing(row_ids) || !identical(row_ids, data$row_ids$survey) ||
        !is.matrix(Xlambda) || !identical(dim(Xlambda), abundance_dimension) ||
        !(is.double(Xlambda) || is.integer(Xlambda)) || is.object(Xlambda) ||
        any(!is.finite(Xlambda)) ||
        !identical(colnames(Xlambda), abundance_columns) ||
        !identical(rownames(Xlambda), abundance_rows))
      stop("Refresh the original abundance columns and row labels with exact survey identities.", call. = FALSE)
    # The complete observation layout was checked once. Only this numeric
    # matrix changes through the API; callers cannot update response fields.
    force(Xlambda); force(row_ids)
    objective <- function(p) {
      eta <- as.vector(Xlambda %*% p$beta_lambda)
      if (length(data$block)) eta <- eta + if (is.null(.lambda_loading))
        p$b_lambda[data$block] else as.vector(.lambda_loading %*% p$b_lambda)
      log_lambda <- eta + log(data$area)
      eta_dist <- as.vector(data$Xdist %*% p$beta_dist)
      if (length(data$observer)) eta_dist <- eta_dist + p$b_dist[data$observer]
      a <- .5 * (radius * exp(-eta_dist))^2
      pd <- lapply(seq_along(width), function(j)
        exp(-a * lower[j]^2) * width[j] * ratio(a * width[j]))
      pd <- do.call(cbind, pd)
      eta_rem <- RTMB::matrix(as.vector(data$Xrem %*% p$beta_rem), M, length(data$periods), byrow = TRUE)
      survived <- eta * 0
      pr <- vector("list", length(data$periods))
      for (j in seq_along(pr)) {
        current <- data$periods[j] * RTMB::logspace_add(eta * 0, eta_rem[, j])
        pr[[j]] <- exp(-survived) * (-expm1(-current))
        survived <- survived + current
      }
      pr <- do.call(cbind, pr)
      q <- RTMB::rowSums(pd) * (-expm1(-survived))
      # Keep one observed total, native margin constant and normalized fields.
      log_cutoff <- cdf(log_lambda + log1p(-q))
      ll <- sum(n * log_lambda - exp(log_lambda) * q + constant + log_cutoff) +
        sum(data$yd[pd_positive] * log(pd[pd_positive])) +
        sum(data$yr[pr_positive] * log(pr[pr_positive]))
      if (length(data$observer))
        ll <- ll + sum(RTMB::dnorm(p$b_dist, 0, exp(p$lsigma_dist), log = TRUE))
      if (length(data$block))
        ll <- ll + sum(RTMB::dnorm(p$b_lambda, 0, exp(p$lsigma_lambda), log = TRUE))
      -ll
    }
    list(evaluate = objective, parameter_map = layout,
      fixed_coordinate_count = sum(kind != "latent_state"),
      latent_coordinate_count = sum(kind == "latent_state"),
      unused_states = c(lambda = length(parameters$b_lambda) - length(unique(data$block)),
        dist = length(parameters$b_dist) - length(unique(data$observer))),
      fingerprint = digest::digest(list(response_context = response_fingerprint,
        Xlambda = Xlambda), algo = "sha256"),
      response_fingerprint = response_fingerprint, shared_observation_preparation = TRUE,
      fingerprint_is_not_full_provider_provenance = TRUE, format_version = 2L,
      target = "P_finite_K_native_margins_normalized_iid_Gaussian_joint",
      coordinate_count_is_not_AIC_admission = TRUE,
      requires_serial_R_callbacks = TRUE, portable_objective = FALSE,
      fit_admission_established = FALSE, uncertainty_evaluated = FALSE,
      prediction_recipe_qualified = FALSE,
      qualification = "private_prepared_numeric_derivative_prerequisite_only")
  }
  list(prepare = prepare, response_fingerprint = response_fingerprint,
    survey_count = M, abundance_columns = abundance_columns,
    shared_observation_preparation = TRUE, context_format_version = 1L,
    requires_serial_R_callbacks = TRUE, portable_objective = FALSE,
    fit_admission_established = FALSE, prediction_recipe_qualified = FALSE,
    qualification = "private_numeric_context_not_fitted_model_or_prediction_recipe")
}

#' Private one-program convenience interface
#' @noRd
.msr_gdr_gaussian_program <- function(data, parameters) {
  context <- .msr_gdr_gaussian_context(data, parameters)
  context$prepare(data$Xlambda, data$row_ids$survey)
}
