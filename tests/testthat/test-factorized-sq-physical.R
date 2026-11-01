test_that("physical SQ coordinates preserve normalized values and every nuisance derivative", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE)
        for (sd in c(.3, 2)) {
          s <- private_sq_fixture(kernel, family, response = response)$spec
          if (family == "binomial") s$trials[] <- 4
          before <- s; g <- .msr_factorized_group(s, 2L)
          p <- private_sq_parameters(g); p$beta[3] <- .12
          p$log_scale_sd[] <- log(sd); p$scale_z <- p$response_z <- numeric()
          normal <- .msr_factorized_group_integral(g)
          physical <- .msr_factorized_group_scale_integral(g)
          normal_ad <- .msr_factorized_group_integral(g, TRUE)
          physical_ad <- .msr_factorized_group_scale_integral(g, TRUE)
          make <- function(fn) RTMB::MakeADFun(function(q) {
            part <- fn(q); -log(part$value) - part$shift
          }, p, silent = TRUE)
          a <- make(normal_ad); b <- make(physical_ad)
          native <- function(x) sq_group_native_integral(g, private_sq_unpack(x, p))
          nv <- normal(p); pv <- physical(p)
          expect_equal(-log(pv$value) - pv$shift, native(b$par), tolerance = 1e-6)
          expect_equal(-log(pv$value) - pv$shift, -log(nv$value) - nv$shift, tolerance = 1e-6)
          expect_equal(as.numeric(b$fn(b$par)), native(b$par), tolerance = 1e-6)
          expect_equal(as.numeric(b$gr(b$par)), private_sq_fd_gradient(native, b$par), tolerance = 2e-5)
          expect_equal(unname(b$he(b$par)), private_sq_fd_hessian(native, b$par), tolerance = 2e-5)
          expect_equal(as.numeric(b$gr(b$par)), as.numeric(a$gr(a$par)), tolerance = 2e-5)
          expect_equal(unname(b$he(b$par)), unname(a$he(a$par)), tolerance = 2e-5)
          expect_identical(g$center, s$center)
          expect_identical(g$spread, s$spread)
          expect_identical(s, before)
          expect_error(.msr_factorized_scale_objective(s), "Weighted SQ")
        }
})

test_that("physical coordinate quadrature preserves the original physical units and response covariance", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp")) {
    f <- private_sq_fixture(kernel, response = TRUE)
    s <- f$spec; g <- .msr_factorized_group(s, 2L)
    moved <- f$inputs
    moved$d_list <- lapply(moved$d_list, function(d) d * moved$unit_conv)
    moved$unit_conv <- 1
    other <- .msr_scale_prepare(f$formula, moved, f$data, f$scales, "gaussian", "key", reference_scale = 250)
    gm <- .msr_factorized_group(other, 2L)
    p <- private_sq_parameters(g); p$beta[3] <- .12
    p$log_scale_sd[] <- log(2); p$scale_z <- p$response_z <- numeric()
    q <- p; q$scale_beta <- p$scale_beta + log(s$unit_conv)
    a <- .msr_factorized_group_scale_integral(g, TRUE)
    b <- .msr_factorized_group_scale_integral(gm, TRUE)
    make <- function(fn, parameters) RTMB::MakeADFun(function(z) {
      part <- fn(z); -log(part$value) - part$shift
    }, parameters, silent = TRUE)
    oa <- make(a, p); ob <- make(b, q)
    expect_equal(gm$center, g$center, tolerance = 1e-12)
    expect_equal(gm$spread, g$spread, tolerance = 1e-12)
    expect_identical(gm$observation_ids, g$observation_ids)
    # Prepared numeric columns can round differently after changing units.
    # All grouping/formula/prediction metadata must remain exactly identical.
    response_metadata <- gm$response; response_metadata$frame <- g$response$frame
    expect_identical(response_metadata, g$response)
    expect_equal(gm$response$frame, g$response$frame, tolerance = 1e-12)
    expect_equal(as.numeric(ob$fn(ob$par)), as.numeric(oa$fn(oa$par)), tolerance = 1e-6)
    expect_equal(as.numeric(ob$gr(ob$par)), as.numeric(oa$gr(oa$par)), tolerance = 2e-5)
    expect_equal(unname(ob$he(ob$par)), unname(oa$he(oa$par)), tolerance = 2e-5)
  }
})

test_that("physical component retains strict group and pre-tape RAM failure behavior", {
  skip_if_not_installed("RTMB")
  s <- private_sq_fixture()$spec; g <- .msr_factorized_group(s, 1L)
  expect_error(.msr_factorized_group_scale_integral(s), "one recoded")
  bad <- g; bad$observation_ids <- s$observation_ids
  expect_error(.msr_factorized_group_scale_integral(bad), "observation metadata")
  launched <- 0L
  constructor <- function() {
    fn <- .msr_factorized_group_scale_integral(g, TRUE)
    RTMB::MakeADFun(function(q) fn(q)$value, g$parameters, silent = TRUE)
  }
  for (ram in c(NA_real_, 1024^3)) with_mocked_bindings({
    with_mocked_bindings({
      expect_error(constructor(), if (is.na(ram)) "Cannot detect RAM" else "RAM budget")
    }, MakeADFun = function(...) { launched <<- launched + 1L; stop("Tape reached") }, .package = "RTMB")
    expect_true(is.function(.msr_factorized_group_scale_integral(g)))
  }, .msr_total_ram_bytes = function() ram, .package = "multiScaleR")
  expect_equal(launched, 0L)
  fixed <- .msr_factorized_kernel_limit(s, "nearest")
  expect_error(.msr_factorized_group_scale_integral(fixed), "one scale grouping")
  expect_equal(.msr_factorized_scale_objective(fixed)(fixed$parameters),
    .msr_factorized_objective(fixed)(fixed$parameters), tolerance = 1e-12)
  expect_error(.msr_factorized_scale_objective(s), "Weighted SQ")
})
