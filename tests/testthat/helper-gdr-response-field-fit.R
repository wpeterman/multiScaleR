# Reuse the published native mixed fixture. Its outcomes, seed and likelihood
# are unchanged; this adds an explicit fixed spatial covariance to abundance.
gdr_response_field_fit_fixture <- function(K = 400L) {
  f <- gdr_tmb_refit_fixture("mixed")
  native <- .msr_gdr_TMB_preparation(f$model)
  # parList's first argument contains FIXED coordinates, not the full joint
  # vector. Use its explicit full-vector argument and the already qualified
  # constructor map, since TMB reorders parameter blocks for its native tape.
  pars <- f$model@TMB$env$parList(par = f$model@TMB$env$last.par)
  random <- c("b_lambda", "b_dist")
  constructor_map <- gdr_start_values(f$model)$map
  for (i in seq_len(nrow(constructor_map)))
    pars[[constructor_map$block[i]]][constructor_map$within_block[i]] <-
      unname(f$constructor_args$starts[i])
  pars$b_lambda[] <- 0; pars$b_dist[] <- 0
  native$K <- as.integer(K)
  map <- function(Z) {
    z <- Matrix::summary(Z)
    stopifnot(nrow(z) == nrow(Z), all(z$x == 1), !anyDuplicated(z$i))
    as.integer(z$j[match(seq_len(nrow(Z)), z$i)])
  }
  ids <- paste0("mixed-survey-", seq_len(nrow(f$model@data@yDistance)))
  block <- map(native$Z_lambda)
  locations <- cbind(430000 + 40000 * (seq_along(pars$b_lambda) - 1), 4400000)
  coords <- locations[block, , drop = FALSE]; rownames(coords) <- ids
  data <- list(Xlambda = native$X_lambda, Xdist = native$X_dist, Xrem = native$X_rem,
    yd = unname(f$model@data@yDistance), yr = unname(f$model@data@yRemoval),
    block = block, observer = map(native$Z_dist), periods = c(1, 2, 1),
    breaks = c(0, 25, 50, 75), area = native$A, K = K,
    annulus_area_m2 = native$a, annulus_weights = native$u,
    row_ids = list(survey = ids, removal_survey = rep(ids, each = 3L),
      removal_interval = rep(1:3, length(ids))))
  list(data = data, pars = pars, native = native, random = random,
    coords = coords, locations = locations, crs = 32617, cold = f$starts$cold,
    constructor_map = constructor_map)
}

# Temporary small matrices remain in memory for reference comparisons only.
# Persist compact diagnostics for every start, including failures, never tapes.
gdr_response_field_fit_attempt <- function(f, context, engine, cold = FALSE) {
  record <- list(engine = engine, cold = cold, status = "started", warnings = character(),
    error = NULL, public_field_admission = FALSE, field_id = context$field_id)
  path <- paste0(gdr_test_audit_path(paste0("spatial-fit-", engine, "-")), ".rds")
  on.exit(.msr_gdr_atomic_audit(record, path))
  withCallingHandlers(tryCatch({
    pars <- f$pars
    if (cold) {
      for (i in seq_len(nrow(f$constructor_map)))
        pars[[f$constructor_map$block[i]]][f$constructor_map$within_block[i]] <-
          unname(f$cold[i])
    }
    if (engine == "RTMB") {
      program <- context$prepare(f$data$Xlambda, f$data$row_ids$survey)
      objective <- RTMB::MakeADFun(program$evaluate, pars, random = f$random, silent = TRUE)
    } else {
      native <- f$native
      L <- t(chol(exp(-as.matrix(dist(f$locations)) / context$field$range)))
      native$Z_lambda <- Matrix::Matrix(native$Z_lambda %*% L, sparse = TRUE)
      objective <- TMB::MakeADFun(native, pars, random = f$random,
        DLL = "unmarked_TMBExports", silent = TRUE)
    }
    fit <- stats::optim(objective$par, objective$fn, objective$gr, method = "BFGS",
      control = list(maxit = 2500L, reltol = 1e-10))
    record$optimizer <- .msr_gdr_optimizer_checkpoint(fit)
    value <- objective$fn(fit$par); score <- as.vector(objective$gr(fit$par))
    record$value <- value; record$maximum_score <- max(abs(score))
    report <- if (engine == "RTMB") RTMB::sdreport(objective, par.fixed = fit$par,
      getJointPrecision = TRUE) else TMB::sdreport(objective, par.fixed = fit$par,
      getJointPrecision = TRUE)
    objective$gr(fit$par)
    point <- objective$env$last.par
    precision <- as.matrix(objective$env$spHess(point, random = TRUE))
    covariance <- as.matrix(report$cov.fixed)
    joint <- as.matrix(report$jointPrecision)
    positive <- function(x) all(is.finite(x)) && isTRUE(tryCatch({chol(x); TRUE},
      error = function(e) FALSE))
    record$pdHess <- report$pdHess
    record$fixed_covariance_positive <- positive(covariance)
    record$joint_precision_positive <- positive(joint)
    record$current_latent_precision_positive <- positive(precision)
    record$coordinate_count <- length(fit$par)
    record$latent_count <- length(objective$env$random)
    record$status <- "returned"
    # sdreport joint precision follows original full parameter order. Retain
    # block/occurrence labels rather than assuming native and RTMB tape orders.
    list(record = record, fixed = fit$par, score = score, point = point,
      random = objective$env$random, precision = precision, covariance = covariance,
      joint_covariance = solve(joint), audit_file = path)
  }, error = function(e) {
    record$status <<- "failed"; record$error <<- conditionMessage(e); stop(e)
  }), warning = function(w) record$warnings <<- c(record$warnings, conditionMessage(w)))
}
