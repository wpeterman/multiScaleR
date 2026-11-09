gdr_gaussian_fixture <- function(block, observer, K, unused = TRUE) {
  skip_if_not_installed("unmarked")
  skip_if(as.character(packageVersion("unmarked")) != "1.5.2", "Native layout needs qualification")
  n <- c(0, 1, 2, 0, 1, 2, 0, 1)
  yd <- cbind(pmin(n, 1), pmax(n - 1, 0))
  yr <- cbind(pmin(n, 1), pmax(n - 1, 0), 0)
  sc <- data.frame(h = c(-.1, .1, .3, .5, .65, .8, 1, 1.2),
    cover = factor(rep(c("A", "B"), each = 4)),
    Observer = factor(rep(c("X", "Y"), 4)),
    block = factor(rep(c("a", "b"), 4), levels = c("a", "b", "unused")))
  contrasts(sc$cover) <- contr.sum(2)
  frame <- unmarked::unmarkedFrameGDR(yd, yr, siteCovs = sc,
    dist.breaks = c(0, 25, 75), unitsIn = "m", period.lengths = c(1, 2, 1))
  ids <- paste0("toy", seq_len(8))
  frame <- expand_gdr_intervals(frame, data.frame(survey = ids), "survey", "minutes",
    site_covariates = c(rem_hab = "h"))$data
  forms <- list(lambda = ~cover * splines::ns(h, knots = c(.2, .6), Boundary.knots = c(0, 1)),
    phi = ~1, dist = if (observer) ~h + (1 | Observer) else ~h,
    rem = ~rem_hab + I(interval == 1))
  if (block) forms$lambda <- update(forms$lambda, ~. + (1 | block))
  original <- get("fit_TMB", asNamespace("unmarked"))
  constructor <- unmarked::gdistremoval
  local <- new.env(parent = environment(constructor))
  local$fit_TMB <- function(model, data, params, random, starts, method, ...)
    stop(structure(list(message = "native numeric preparation only", call = NULL,
      model = model, data = data, parameters = params, random = random),
      class = c("msr_test_gaussian_GDR", "error", "condition")))
  environment(constructor) <- local
  captured <- tryCatch(constructor(lambdaformula = forms$lambda, phiformula = forms$phi,
    distanceformula = forms$dist, removalformula = forms$rem, data = frame,
    K = K, keyfun = "halfnorm", mixture = "P", output = "density", unitsOut = "kmsq",
    engine = "TMB", threads = 1), msr_test_gaussian_GDR = identity)
  stopifnot(inherits(captured, "msr_test_gaussian_GDR"),
    identical(get("fit_TMB", asNamespace("unmarked")), original))
  if (block && unused) {
    captured$data$Z_lambda <- cbind(captured$data$Z_lambda,
      Matrix::Matrix(0, 8L, 1L, sparse = TRUE))
    captured$data$n_grouplevels_lambda <- captured$data$n_grouplevels_lambda + 1L
    captured$parameters$b_lambda <- c(captured$parameters$b_lambda, 0)
  }
  gd <- get("getDesign", asNamespace("unmarked"))(frame, forms)
  map <- function(Z) {
    if (!ncol(Z)) return(integer())
    z <- Matrix::summary(Z)
    stopifnot(nrow(z) == 8L, all(z$x == 1), !anyDuplicated(z$i))
    as.integer(z$j[match(seq_len(8L), z$i)])
  }
  pars <- captured$parameters
  pars$beta_lambda <- seq_len(ncol(gd$X_lambda)) * .01
  pars$beta_lambda[1L] <- log(8/(pi * .075^2))
  pars$beta_dist <- c(log(40), .1)
  pars$beta_rem <- c(qlogis(.25), .1, .2)
  pars$b_lambda <- if (block) c(.1, -.2, if (unused) .3) else numeric()
  pars$b_dist <- if (observer) c(.08, -.05) else numeric()
  pars$lsigma_lambda <- if (block) log(.3) else numeric()
  pars$lsigma_dist <- if (observer) log(.4) else numeric()
  data <- list(Xlambda = gd$X_lambda, Xdist = gd$X_dist, Xrem = gd$X_rem,
    yd = unname(frame@yDistance), yr = unname(frame@yRemoval),
    block = map(gd$Z_lambda), observer = map(gd$Z_dist),
    periods = c(1, 2, 1), breaks = c(0, 25, 75), area = captured$data$A, K = K,
    annulus_area_m2 = captured$data$a, annulus_weights = captured$data$u,
    row_ids = list(survey = ids, removal_survey = rep(ids, each = 3L),
      removal_interval = rep(1:3, 8L)))
  list(data = data, pars = pars, captured = captured)
}

gdr_gaussian_require_RTMB <- function() {
  skip_if_not_installed("RTMB")
  skip_if(as.character(packageVersion("RTMB")) != "2.0", "Private derivative layout needs qualification")
}

gdr_gaussian_numeric <- function(data, p, uniform = FALSE) {
  eta <- as.vector(data$Xlambda %*% p$beta_lambda)
  if (length(data$block)) eta <- eta + p$b_lambda[data$block]
  eta_dist <- as.vector(data$Xdist %*% p$beta_dist)
  if (length(data$observer)) eta_dist <- eta_dist + p$b_dist[data$observer]
  pd <- .msr_gdr_exact_distance_bins(exp(eta_dist), data$breaks,
    if (uniform) "uniform" else "halfnorm")
  removal <- matrix(plogis(as.vector(data$Xrem %*% p$beta_rem)), 8L, 3L, byrow = TRUE)
  pr <- .msr_gdr_exact_removal_bins(removal, data$periods)
  ll <- .msr_gdr_finite_poisson_score_bins(exp(eta) * data$area,
    pd, pr, data$yd, data$yr, data$K)$logLik
  if (length(data$block)) ll <- ll + sum(dnorm(p$b_lambda, 0, exp(p$lsigma_lambda), log = TRUE))
  if (length(data$observer)) ll <- ll + sum(dnorm(p$b_dist, 0, exp(p$lsigma_dist), log = TRUE))
  -ll
}

test_that("annulus callbacks retain moment orders and update across the numeric branch", {
  gdr_gaussian_require_RTMB()
  ratio <- .msr_gdr_annulus_program()
  tape <- RTMB::MakeTape(function(x) sum(ratio$evaluate(x)), 0)
  for (x in c(0, 1e-12, .05, .099, .101, 1, 10, 1e4)) {
    actual <- c(as.numeric(tape(x)), vapply(1:4, function(order) {
      derivative <- tape
      if (order > 1L) for (j in seq_len(order - 1L)) derivative <- derivative$jacfun()
      as.numeric(derivative$jacobian(x))
    }, numeric(1L)))
    expected <- vapply(0:4, function(order) {
      moment <- if (x == 1e4) factorial(order)/x^(order + 1) else
        integrate(function(t) t^order * exp(-x * t), 0, 1, rel.tol = 2e-12, abs.tol = 0)$value
      (-1)^order * moment
    }, numeric(1L))
    expect_lt(max(abs(actual - expected)), 1e-10)
  }
  expect_false(ratio$portable_objective)
  expect_false(ratio$fit_admission_established)
})

test_that("native finite-K joint targets and derivatives survive exact field reductions", {
  gdr_gaussian_require_RTMB()
  key <- function(x) paste(x, ave(seq_along(x), x, FUN = seq_along), sep = "#")
  for (block in c(FALSE, TRUE)) for (observer in c(FALSE, TRUE)) for (K in c(2L, 40L)) {
    f <- gdr_gaussian_fixture(block, observer, K)
    program <- .msr_gdr_gaussian_program(f$data, f$pars)
    candidate <- RTMB::MakeADFun(program$evaluate, f$pars, random = NULL, silent = TRUE)
    native <- TMB::MakeADFun(data = c(list(model = f$captured$model), f$captured$data),
      parameters = f$pars, random = NULL, DLL = "unmarked_TMBExports", silent = TRUE)
    index <- match(key(names(candidate$par)), key(names(native$par)))
    expect_false(anyNA(index))
    expect_lt(abs(candidate$fn(candidate$par) - native$fn(native$par)), 1e-9)
    expect_lt(abs(candidate$fn(candidate$par) - gdr_gaussian_numeric(f$data, f$pars)), 1e-9)
    expect_lt(max(abs(as.vector(candidate$gr(candidate$par)) - as.vector(native$gr(native$par))[index])), 1e-8)
    expect_lt(max(abs(candidate$he(candidate$par) - native$he(native$par)[index, index, drop = FALSE])), 1e-8)
    expect_identical(program$fixed_coordinate_count, 13L + as.integer(block) + as.integer(observer))
    expect_identical(program$latent_coordinate_count, 3L * as.integer(block) + 2L * as.integer(observer))
    expect_false(program$fit_admission_established)
    expect_false(program$uncertainty_evaluated)
    rm(candidate, native, program); gc(verbose = FALSE)
  }
})

test_that("unused Gaussian coordinates normalize without changing fixed-coordinate counts", {
  gdr_gaussian_require_RTMB()
  values <- lapply(c(TRUE, FALSE), function(unused) {
    f <- gdr_gaussian_fixture(TRUE, TRUE, 2L, unused)
    program <- .msr_gdr_gaussian_program(f$data, f$pars)
    f$pars$b_lambda[] <- f$pars$b_dist[] <- 0
    obj <- RTMB::MakeADFun(program$evaluate, f$pars,
      random = f$captured$random, silent = TRUE)
    list(value = obj$fn(obj$par), score = as.vector(obj$gr(obj$par)),
      count = program$fixed_coordinate_count, unused = program$unused_states,
      fingerprint = program$fingerprint)
  })
  expect_lt(abs(values[[1L]]$value - values[[2L]]$value), 1e-10)
  expect_lt(max(abs(values[[1L]]$score - values[[2L]]$score)), 1e-10)
  expect_identical(values[[1L]]$count, values[[2L]]$count)
  expect_identical(unname(values[[1L]]$unused['lambda']), 1L)
  expect_identical(unname(values[[2L]]$unused['lambda']), 0L)
  expect_false(identical(values[[1L]]$fingerprint, values[[2L]]$fingerprint))
})

test_that("the stable distance limit stays finite without granting uncertainty", {
  gdr_gaussian_require_RTMB()
  f <- gdr_gaussian_fixture(FALSE, FALSE, 40L)
  f$pars$beta_dist <- c(40, 0)
  program <- .msr_gdr_gaussian_program(f$data, f$pars)
  obj <- RTMB::MakeADFun(program$evaluate, f$pars, random = NULL, silent = TRUE)
  expect_lt(abs(obj$fn(obj$par) - gdr_gaussian_numeric(f$data, f$pars, uniform = TRUE)), 1e-9)
  expect_lt(max(abs(obj$gr(obj$par)[names(obj$par) == 'beta_dist'])), 1e-10)
  expect_false(program$uncertainty_evaluated)
  expect_true(program$coordinate_count_is_not_AIC_admission)
})

test_that("prepared layouts reject changed rows, geometry and unsupported parameter blocks", {
  f <- gdr_gaussian_fixture(TRUE, TRUE, 2L)
  check <- function(change, message) {
    changed <- change(f)
    expect_error(.msr_gdr_gaussian_data(changed$data, changed$pars), message)
  }
  check(function(x) { x$data$K <- 1L; x }, "cutoff")
  check(function(x) { x$data$yd[1, 1] <- 1; x }, "equal totals")
  check(function(x) { x$data$periods[2] <- 0; x }, "durations")
  check(function(x) { x$data$Xrem <- x$data$Xrem[-1, , drop = FALSE]; x }, "design matrices")
  check(function(x) { x$data$Xlambda[1, 1] <- NA; x }, "design matrices")
  check(function(x) { x$data$row_ids$survey[2] <- x$data$row_ids$survey[1]; x }, "identities")
  check(function(x) { x$data$row_ids$removal_survey <- rev(x$data$row_ids$removal_survey); x }, "ordering")
  check(function(x) { x$data$breaks <- c(0, 30, 75); x }, "geometry")
  check(function(x) { x$data$area[1] <- x$data$area[1] * 1e6; x }, "square-kilometer")
  check(function(x) { x$data$annulus_weights[1, 1] <- .2; x }, "annulus")
  check(function(x) { x$data$block[1] <- 0; x }, "row map")
  check(function(x) { x$pars$lsigma_lambda <- numeric(); x }, "log-SD")
  check(function(x) { x$pars$beta_lambda <- x$pars$beta_lambda[-1]; x }, "design column")
  check(function(x) { x$pars$beta_phi <- 0; x }, "Annual availability")
  check(function(x) { x$pars$b_rem <- 0; x }, "removal fields")
})

test_that("data values and identities invalidate a derivative context fingerprint", {
  gdr_gaussian_require_RTMB()
  f <- gdr_gaussian_fixture(TRUE, TRUE, 2L)
  original <- .msr_gdr_gaussian_program(f$data, f$pars)
  changed <- f$data
  changed$Xlambda[1, 2] <- changed$Xlambda[1, 2] + .01
  expect_false(identical(original$fingerprint, .msr_gdr_gaussian_program(changed, f$pars)$fingerprint))
  changed <- f$data
  changed$row_ids$survey <- paste0('new', changed$row_ids$survey)
  changed$row_ids$removal_survey <- rep(changed$row_ids$survey, each = 3L)
  expect_false(identical(original$fingerprint, .msr_gdr_gaussian_program(changed, f$pars)$fingerprint))
  expect_true(original$fingerprint_is_not_full_provider_provenance)
  cat("\nGDR_GAUSSIAN_JOINT prerequisites=1 native_P_HN_only=1 fitting_admission=0\n")
})
