# Layout witnesses only. Native TMB fit_TMB is intercepted before optimization;
# joint tapes fix every coordinate and are never evaluated or checkpointed.
gdr_start_frame <- function(primary) {
  n <- 4L
  yd <- matrix(0, n, 2L * primary)
  yr <- matrix(0, n, 3L * primary)
  yd[, 1L] <- yr[, 1L] <- 1
  if (primary == 2L) yd[, 3L] <- yr[, 4L] <- 1
  site <- data.frame(hab = seq(-.3, .3, length.out = n),
    cover = factor(c("None", "A", "None", "A")),
    Observer = factor(c("X", "X", "Y", "Y")))
  contrasts(site$cover) <- stats::contr.sum(2)
  frame <- unmarked::unmarkedFrameGDR(yd, yr, numPrimary = primary,
    siteCovs = site, yearlySiteCovs = data.frame(survey_wind = seq(-.2, .2, length.out = n * primary)),
    dist.breaks = c(0, 25, 75), unitsIn = "m", period.lengths = c(1, 2, 1))
  expand_gdr_intervals(frame, data.frame(survey = paste0("s", seq_len(n * primary))),
    "survey", "minutes", site_covariates = c(rem_hab = "hab", rem_observer = "Observer"),
    yearly_covariates = c(rem_wind = "survey_wind"))$data
}

gdr_start_forms <- function(random = FALSE) list(lambda = ~hab * cover,
  phi = ~survey_wind, dist = if (random) ~hab + (1 | Observer) else ~hab,
  rem = if (random) ~rem_hab * factor(interval) + I(rem_wind^2) + (1 | rem_observer) else
    ~rem_hab * factor(interval) + I(rem_wind^2))

gdr_start_C_points <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    frame <- gdr_start_frame(2L)
    forms <- gdr_start_forms()
    design <- get("getDesign", asNamespace("unmarked"))(frame, forms)
    cache <<- lapply(c("uniform", "hazard"), function(key) {
      mixture <- if (key == "uniform") "P" else "NB"
      counts <- c(lambda = ncol(design$X_lambda), alpha = as.integer(mixture != "P"),
        phi = ncol(design$X_phi), dist = if (key == "uniform") 0L else ncol(design$X_dist),
        scale = as.integer(key == "hazard"), rem = ncol(design$X_rem))
      values <- lapply(names(counts), function(name) {
        x <- seq_len(counts[[name]]) * .01
        if (length(x)) x[1L] <- switch(name, lambda = log(8), alpha = log(2),
          phi = qlogis(.6), dist = log(40), scale = log(2), rem = qlogis(.3))
        x
      })
      unmarked::gdistremoval(lambdaformula = forms$lambda, phiformula = forms$phi,
        distanceformula = forms$dist, removalformula = forms$rem, data = frame,
        keyfun = key, mixture = mixture, output = "abund", K = 40, threads = 1,
        starts = unlist(values, use.names = FALSE), engine = "C", se = FALSE,
        control = list(maxit = 0))
    })
    names(cache) <<- c("uniform", "hazard")
    cache
  }
})

gdr_start_TMB_witnesses <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    cache <<- lapply(c("single_NB", "uniform_P", "random_NB"), function(kind) {
      frame <- gdr_start_frame(if (kind == "single_NB") 1L else 2L)
      forms <- gdr_start_forms(kind == "random_NB")
      key <- if (kind == "uniform_P") "uniform" else "halfnorm"
      mixture <- if (kind == "uniform_P") "P" else "NB"
      original <- get("fit_TMB", asNamespace("unmarked"))
      captured <- testthat::with_mocked_bindings(
        tryCatch(unmarked::gdistremoval(lambdaformula = forms$lambda, phiformula = forms$phi,
          distanceformula = forms$dist, removalformula = forms$rem, data = frame,
          keyfun = key, mixture = mixture, output = "abund", K = 40, engine = "TMB"),
          msr_test_GDR_layout = function(e) e),
        fit_TMB = function(model, data, params, random, starts, method, ...)
          stop(structure(list(message = "native input layout only; optimizer not entered", call = NULL,
            model_name = model, data = data, params = params, random = random),
            class = c("msr_test_GDR_layout", "error", "condition"))), .package = "unmarked")
      stopifnot(inherits(captured, "msr_test_GDR_layout"),
        identical(get("fit_TMB", asNamespace("unmarked")), original))
      for (name in names(captured$params)) {
        x <- captured$params[[name]]
        if (!length(x)) next
        x[] <- seq_along(x) * .01
        if (startsWith(name, "beta_")) x[1L] <- switch(sub("beta_", "", name),
          lambda = log(8), phi = qlogis(.6), dist = log(40), rem = qlogis(.3), alpha = log(2), scale = log(2))
        if (startsWith(name, "lsigma_")) x[] <- log(.4 + seq_along(x) / 10)
        captured$params[[name]] <- x
      }
      # random=NULL creates a joint fixed-coordinate tape, not a Laplace fit.
      joint <- TMB::MakeADFun(data = c(list(model = captured$model_name), captured$data),
        parameters = captured$params, random = NULL, silent = TRUE, DLL = "unmarked_TMBExports")
      fixed <- joint$par[!names(joint$par) %in% captured$random]
      constructor <- rep(names(captured$params)[!names(captured$params) %in% captured$random],
        lengths(captured$params[!names(captured$params) %in% captured$random]))
      design <- get("getDesign", asNamespace("unmarked"))(frame, forms)
      # This S4 object is explicitly a wrapper metadata fixture, NOT a fit.
      metadata <- gdr_guard_fixture()$NB
      metadata@data <- frame
      metadata@formlist <- forms
      metadata@mixture <- mixture
      metadata@keyfun <- key
      active <- c("lambda", if (mixture != "P") "alpha", if (frame@numPrimary > 1) "phi",
        if (key != "uniform") "dist", "rem")
      components <- stats::setNames(lapply(active, function(name) {
        beta <- captured$params[[paste0("beta_", name)]]
        modes <- captured$params[[paste0("b_", name)]]
        sd <- captured$params[[paste0("lsigma_", name)]]
        labels <- if (name == "alpha") "alpha" else colnames(design[[paste0("X_", name)]])
        values <- c(stats::setNames(beta, labels),
          if (length(modes)) stats::setNames(modes, paste0("mode", seq_along(modes))))
        table <- switch(name, lambda = unmarked::siteCovs(frame), phi = unmarked::yearlySiteCovs(frame),
          dist = unmarked::siteCovs(frame), rem = unmarked::obsCovs(frame))
        random <- if (length(sd)) list(estimates = sd,
          names = get("sigma_names", asNamespace("unmarked"))(forms[[name]], table),
          covMat = matrix(NA_real_, length(sd), length(sd))) else list()
        methods::new("unmarkedEstimate", name = name, short.name = name, estimates = values,
          fixed = seq_along(beta), covMat = matrix(NA_real_, length(values), length(values)),
          invlink = "exp", invlinkGrad = "exp", randomVarInfo = random)
      }), active)
      metadata@estimates <- methods::new("unmarkedEstimateList", estimates = components)
      metadata@opt <- list(par = fixed, value = Inf, convergence = 1L)
      metadata@nllFun <- function(...) stop("must not evaluate saved likelihood")
      metadata@TMB <- list(starts_order = constructor,
        fn = function(...) stop("must not evaluate saved tape"))
      list(kind = kind, model = metadata, params = captured$params, random = captured$random,
        joint_blocks = names(joint$par), fixed = fixed, constructor = constructor,
        qualifier = "native_input_and_joint_coordinate_witness_plus_S4_metadata_fixture_not_fit")
    })
    names(cache) <<- c("single_NB", "uniform_P", "random_NB")
    cache
  }
})

test_that("actual C points retain all active constructor coordinates without callbacks", {
  for (model in c(gdr_guard_fixture(), gdr_start_C_points())) {
    before <- gdr_guard_state(model)
    snapshot <- gdr_start_values(model)
    expect_identical(snapshot$starts, unname(model@opt$par))
    expect_identical(snapshot$fitted_values, model@opt$par)
    expect_identical(snapshot$coordinate_count, length(model@opt$par))
    expect_identical(snapshot$provider$engine, "C")
    expect_identical(snapshot$map$constructor_position, seq_along(model@opt$par))
    expect_identical(snapshot$map$fitted_position, seq_along(model@opt$par))
    expect_false(snapshot$fit_eligibility_established)
    expect_false(snapshot$likelihood_evaluated)
    expect_false(snapshot$uncertainty_evaluated)
    expect_identical(gdr_guard_state(model), before)
    expect_identical(unserialize(serialize(snapshot, NULL, version = 2)), snapshot)
    restored <- numeric(snapshot$coordinate_count)
    restored[snapshot$map$fitted_position] <- snapshot$starts
    expect_identical(restored, unname(model@opt$par))
  }
  uniform <- gdr_start_values(gdr_start_C_points()$uniform)
  expect_false(any(uniform$map$component == "dist"))
  expect_true(any(uniform$map$component == "phi"))
  hazard <- gdr_start_values(gdr_start_C_points()$hazard)
  expect_true(any(hazard$map$component == "scale"))
  expect_true(any(hazard$map$within_block > 1L))
})

test_that("actual TMB joint order maps to intercepted constructor input order", {
  for (witness in gdr_start_TMB_witnesses()) {
    model <- witness$model
    before <- gdr_guard_state(model)
    snapshot <- gdr_start_values(model)
    expected <- unlist(witness$params[!names(witness$params) %in% witness$random], use.names = FALSE)
    expect_identical(snapshot$starts, unname(expected))
    expect_identical(snapshot$fitted_values, witness$fixed)
    expect_identical(snapshot$map$block, witness$constructor)
    restored <- numeric(length(snapshot$fitted_values))
    restored[snapshot$map$fitted_position] <- snapshot$starts
    expect_identical(restored, unname(witness$fixed))
    expect_false(snapshot$fit_eligibility_established)
    expect_false(snapshot$likelihood_evaluated)
    expect_false(snapshot$uncertainty_evaluated)
    expect_identical(snapshot$optimizer$value, Inf)
    expect_identical(snapshot$optimizer$convergence, 1L)
    expect_identical(gdr_guard_state(model), before)
    expect_identical(unserialize(serialize(snapshot, NULL, version = 2)), snapshot)
    plain <- function(x) {
      if (is.environment(x) || is.function(x) || isS4(x) || typeof(x) == "externalptr") return(FALSE)
      if (is.list(x)) return(all(vapply(x, plain, logical(1))))
      TRUE
    }
    expect_true(plain(snapshot))
  }
  examples <- gdr_start_TMB_witnesses()
  nb <- gdr_start_values(examples$single_NB$model)
  expect_false(identical(nb$map$constructor_position, nb$map$fitted_position))
  expect_true(any(nb$map$block == "beta_alpha"))
  uniform <- gdr_start_values(examples$uniform_P$model)
  expect_false(any(uniform$map$component %in% c("dist", "alpha")))
  expect_true(any(uniform$map$component == "phi"))
  random <- gdr_start_values(examples$random_NB$model)
  expect_setequal(random$map$component[random$map$kind == "log_sd"], c("dist", "rem"))
  expect_true(all(random$random_modes_omitted[c("dist", "rem")] > 0L))
  expect_false(any(startsWith(random$map$block, "b_")))
  cat("\nGDR_START_METADATA C_points=5 TMB_layouts=3 fit_admission=0 latent_optimization=0\n")
})

test_that("inconsistent or incomplete records fail without repairing native values", {
  cmodel <- gdr_guard_fixture()$NB
  before <- gdr_guard_state(cmodel)
  bad <- cmodel
  bad@opt$par[1] <- NA_real_
  expect_error(gdr_start_values(bad), "finite")
  bad <- cmodel
  bad@estimates@estimates$lambda@estimates[1] <- .1
  expect_error(gdr_start_values(bad), "disagree")
  bad <- cmodel
  bad@estimates@estimates$lambda@fixed <- c(2, 1)
  expect_error(gdr_start_values(bad), "noncanonical")
  bad <- cmodel
  names(bad@estimates@estimates$lambda@estimates) <- NULL
  expect_error(gdr_start_values(bad), "labels")
  bad <- cmodel
  bad@estimates@estimates$extra <- bad@estimates@estimates$lambda
  expect_error(gdr_start_values(bad), "complete active")
  bad <- cmodel
  bad@formlist$lambda <- ~hab + offset(exposure)
  expect_error(gdr_start_values(bad), "offsets are not qualified")
  expect_error(gdr_start_values(list()), "native GDR coordinate")
  expect_identical(gdr_guard_state(cmodel), before)
  tmb <- gdr_start_TMB_witnesses()$random_NB$model
  for (blocks in list(NULL, "beta_lambda", c(tmb@TMB$starts_order, "b_rem"))) {
    bad <- tmb
    bad@TMB$starts_order <- blocks
    expect_error(gdr_start_values(bad), "identities|bijection|random modes")
  }
  bad <- tmb
  names(bad@opt$par) <- NULL
  expect_error(gdr_start_values(bad), "identities")
  bad <- tmb
  bad@estimates@estimates$dist@randomVarInfo$estimates[1] <- Inf
  expect_error(gdr_start_values(bad), "disagree")
  bad <- tmb
  bad@estimates@estimates$rem@randomVarInfo$names <- NULL
  expect_error(gdr_start_values(bad), "variance coordinates")
  bad <- tmb
  bad@mixture <- "ZIP"
  expect_error(gdr_start_values(bad), "TMB ZIP")
  bad <- cmodel
  bad@nllFun <- function(...) stop("must not evaluate")
  expect_identical(gdr_start_values(bad)$starts, unname(cmodel@opt$par))
})
