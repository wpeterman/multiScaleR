factorized_refinement_fixture <- function() {
  list(covariates = "x", scales = list(list(random = list(list()))),
    response = list(random = list()),
    parameters = list(beta = c(1, 2), scale_beta = matrix(c(.3, .4), 1),
      scale_z = 0, log_scale_sd = -1, response_z = numeric(),
      log_response_sd = numeric(), log_dispersion = 0))
}

test_that("bounded integration refinement retries the identical point", {
  spec <- factorized_refinement_fixture()
  constructions <- list()
  factory <- function(spec, parameters = spec$parameters, rel.tol = 1e-8, step = 1) {
    constructions[[length(constructions) + 1L]] <<- list(step = step, parameters = parameters)
    if (step == 4) stop("injected construction failure")
    x <- unlist(parameters, use.names = FALSE); e <- new.env()
    e$last.par <- x
    list(par = x, env = e, marker = step,
      fn = function(x) sum(x^2),
      gr = function(x) if (step == 2) rep(NaN, length(x)) else 2*x,
      he = function(x) diag(2, length(x)))
  }
  local_mocked_bindings(.msr_factorized_adfun = factory, .package = "multiScaleR")
  obj <- multiScaleR:::.msr_factorized_optimizer(spec, step = 4)
  x <- seq_along(obj$par) / 10
  expect_equal(obj$fn(x), sum(x^2))
  expect_identical(obj$env$integration_current_step, 2L)
  generation <- obj$env$integration_generation
  expect_equal(obj$gr(x), 2*x)
  expect_identical(obj$env$integration_current_step, 1L)
  expect_gt(obj$env$integration_generation, generation)
  expect_equal(obj$he(x), diag(2, length(x)))
  expect_length(obj$env$integration_events, 2L)
  expect_identical(vapply(obj$env$integration_events, `[[`, "", "phase"), c("fn", "gr"))
  for (event in obj$env$integration_events) expect_identical(event$parameters, x)
  restored <- obj$env$parList(x)
  expect_equal(unlist(restored, use.names = FALSE), x)
  expect_identical(dim(restored$scale_beta), c(1L, 2L))
  expect_length(restored$scale_z, 0L)
  expect_length(constructions, 3L)
  expect_equal(unlist(constructions[[3L]]$parameters, use.names = FALSE), x)
  expect_error(obj$fn(x[-1L]), "Invalid factorized parameter")
  expect_error(obj$gr(replace(x, 1, Inf)), "Invalid factorized parameter")
})

test_that("value and curvature failures refine without inventing derivatives", {
  spec <- factorized_refinement_fixture()
  for (phase in c("fn", "he")) {
    factory <- local({
      failed_phase <- phase
      function(spec, parameters = spec$parameters, rel.tol = 1e-8, step = 1) {
        x <- unlist(parameters, use.names = FALSE); e <- new.env(); e$last.par <- x
        list(par = x, env = e,
          fn = function(x) if (step == 2 && failed_phase == "fn") Inf else sum(x^2),
          gr = function(x) 2*x,
          he = function(x) if (step == 2 && failed_phase == "he") stop("injected Hessian failure") else diag(2, length(x)))
      }
    })
    local_mocked_bindings(.msr_factorized_adfun = factory, .package = "multiScaleR")
    obj <- multiScaleR:::.msr_factorized_optimizer(spec, step = 2)
    x <- as.numeric(obj$par)
    value <- obj[[phase]](x)
    expect_equal(value, if (phase == "fn") sum(x^2) else diag(2, length(x)))
    expect_identical(obj$env$integration_current_step, 1L)
    expect_length(obj$env$integration_events, 1L)
    expect_identical(obj$env$integration_events[[1L]]$phase, phase)
    expect_identical(obj$env$integration_events[[1L]]$parameters, x)
  }
})

test_that("terminal failures remain failures and default integration is direct", {
  spec <- factorized_refinement_fixture()
  local_mocked_bindings(.msr_factorized_adfun = function(...) stop("injected terminal failure"),
    .package = "multiScaleR")
  obj <- multiScaleR:::.msr_factorized_optimizer(spec, step = 8)
  expect_error(obj$fn(obj$par), "injected terminal failure")
  expect_length(obj$env$integration_events, 4L)
  expect_identical(vapply(obj$env$integration_events, `[[`, 0L, "from"), c(8L, 4L, 2L, 1L))
  expect_true(is.na(obj$env$integration_events[[4L]]$to))
  expect_error(multiScaleR:::.msr_factorized_optimizer(spec, step = 1), "injected terminal failure")
  for (bad in list(0, 3, 2.5, 16, NA_real_, "4", c(1, 2)))
    expect_error(multiScaleR:::.msr_factorized_optimizer(spec, step = bad), "Invalid optimizing")
})
