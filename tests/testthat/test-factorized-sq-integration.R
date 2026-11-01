test_that("one-group raw-SQ quadrature matches independent normalized values and full derivatives", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE) {
        s <- private_sq_fixture(kernel, family, response = response)$spec
        if (family == "binomial") s$trials[] <- 4
        before <- s; g <- .msr_factorized_group(s, 2L)
        p <- private_sq_parameters(g); p$beta[3] <- .12
        p$scale_z <- p$response_z <- numeric()
        expect_null(g$priors)
        ordinary <- .msr_factorized_group_integral(g, FALSE)
        traced <- .msr_factorized_group_integral(g, TRUE)
        obj <- RTMB::MakeADFun(function(q) {
          part <- traced(q)
          -log(part$value) - part$shift
        }, p, silent = TRUE)
        native <- function(x) sq_group_native_integral(g, private_sq_unpack(x, p))
        value <- ordinary(p)
        expect_equal(-log(value$value) - value$shift, native(obj$par), tolerance = 1e-6)
        expect_equal(as.numeric(obj$fn(obj$par)), native(obj$par), tolerance = 1e-6)
        expect_equal(as.numeric(obj$gr(obj$par)), private_sq_fd_gradient(native, obj$par), tolerance = 2e-5)
        expect_equal(unname(obj$he(obj$par)), private_sq_fd_hessian(native, obj$par), tolerance = 2e-5)
        expect_identical(g$center, s$center)
        expect_identical(g$spread, s$spread)
        expect_identical(s, before)
        expect_error(.msr_factorized_objective(s, TRUE), "Weighted SQ")
        expect_error(.msr_factorized_adfun(s), "Weighted SQ")
        expect_error(.msr_factorized_score_adfun(s), "Weighted SQ")
        expect_error(.msr_factorized_conditional_cache(s), "Weighted SQ")
      }
})

test_that("group quadrature rejects complete or malformed models before density evaluation", {
  s <- private_sq_fixture(response = TRUE)$spec
  expect_error(.msr_factorized_group_integral(s), "one recoded")
  g <- .msr_factorized_group(s, 1L)
  invalid <- g; invalid$scales[[1L]]$random <- list()
  expect_error(.msr_factorized_group_integral(invalid), "one scale grouping")
  invalid <- g; invalid$scales[[1L]]$random[[1L]]$gaussian_block <- "unsupported"
  expect_error(.msr_factorized_group_integral(invalid), "iid scale group")
  invalid <- g; invalid$rows <- s$rows
  expect_error(.msr_factorized_group_integral(invalid), "observation metadata")
  invalid <- g; invalid$observation_ids <- s$observation_ids
  expect_error(.msr_factorized_group_integral(invalid), "observation metadata")
  invalid <- g; invalid$data <- s$data
  expect_error(.msr_factorized_group_integral(invalid), "observation metadata")
  invalid <- g; invalid$family <- "poisson"
  expect_error(.msr_factorized_group_integral(invalid), "response structure")
  invalid <- g; invalid$response$random[[1L]]$gaussian_block <- "unsupported"
  expect_error(.msr_factorized_group_integral(invalid), "response structure")
  expect_true(is.function(.msr_factorized_group_integral(g)))
})

test_that("raw-SQ group quadrature enforces scientific RAM admission before tape entry", {
  skip_if_not_installed("RTMB")
  s <- private_sq_fixture()$spec; g <- .msr_factorized_group(s, 1L)
  launched <- 0L
  constructor <- function() {
    fn <- .msr_factorized_group_integral(g, TRUE)
    RTMB::MakeADFun(function(p) fn(p)$value, g$parameters, silent = TRUE)
  }
  for (ram in c(NA_real_, 1024^3)) with_mocked_bindings({
    with_mocked_bindings({
      expect_error(constructor(), if (is.na(ram)) "Cannot detect RAM" else "RAM budget")
    }, MakeADFun = function(...) { launched <<- launched + 1L; stop("Tape reached") }, .package = "RTMB")
    expect_true(is.function(.msr_factorized_group_integral(g, FALSE)))
  }, .msr_total_ram_bytes = function() ram, .package = "multiScaleR")
  expect_equal(launched, 0L)
  expect_error(.msr_factorized_validate(s), "Weighted SQ")
})
