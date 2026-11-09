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
