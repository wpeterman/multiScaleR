# A default refit search that never leaves its starting values must say so.
# Two causes share one signature (L-BFGS-B stops after one evaluation with a
# zero finite-difference gradient): every optimized scale is a hard radius, so
# the likelihood changes in steps, or every refit at the start fails and
# returns the penalty value.

running_example <- function(kernel) {
  land <- terra::rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))
  inputs <- kernel_prep(pts = surv_pts, raster_stack = land, max_D = 1700,
                        kernel = kernel, verbose = FALSE)
  dat <- data.frame(landscape_counts, inputs$kernel_dat)
  mod <- glm(counts ~ site + land1 + land2, family = poisson(), data = dat)
  list(inputs = inputs, mod = mod)
}

test_that("a fixed-kernel single start that cannot move warns and is recorded", {
  ex <- running_example("fixed")
  # Since 0.7.42 a fixed kernel screens by default, so request the single
  # start explicitly to exercise the unscreened stall.
  expect_warning(
    fit <- multiScale_optim(ex$mod, ex$inputs, n_cores = 1, verbose = FALSE,
                            start_strategy = "single"),
    "hard radius.*start_strategy = \"screen\""
  )
  stall <- diagnostics(fit)$start_stall
  expect_true(stall$triggered)
  expect_identical(stall$reason, "hard_radius")
  expect_false(stall$screened)
  expect_equal(stall$function_evaluations, 1L)
  # The reported radii are the starting values.
  expect_equal(fit$scale_est$Mean, unname(stall$start))

  # summary() and print() repeat the explanation.
  expect_output(print(summary(fit)), "hard radius")
  expect_output(print(fit), "hard radius")
})

test_that("a screened fixed-kernel start that cannot move is described as screened", {
  ex <- running_example("fixed")
  set.seed(1)
  fit <- withCallingHandlers(
    multiScale_optim(ex$mod, ex$inputs, n_cores = 1, verbose = FALSE,
                     start_strategy = "screen"),
    warning = function(w) {
      expect_match(conditionMessage(w), "screened start")
      invokeRestart("muffleWarning")
    })
  stall <- diagnostics(fit)$start_stall
  # The screened start improves on the default start, then stays put.
  expect_identical(stall$reason, "hard_radius")
  expect_true(stall$screened)
  expect_equal(fit$scale_est$Mean, unname(stall$start))
})

test_that("a fixed kernel screens starting values by default", {
  ex <- running_example("fixed")
  set.seed(1)
  screened <- suppressWarnings(
    multiScale_optim(ex$mod, ex$inputs, n_cores = 1, verbose = FALSE))
  set.seed(1)
  explicit <- suppressWarnings(
    multiScale_optim(ex$mod, ex$inputs, n_cores = 1, verbose = FALSE,
                     start_strategy = "screen"))
  single <- suppressWarnings(
    multiScale_optim(ex$mod, ex$inputs, n_cores = 1, verbose = FALSE,
                     start_strategy = "single"))
  # The default equals an explicit screen and fits better than one start.
  expect_true(diagnostics(screened)$start_stall$screened)
  expect_equal(screened$scale_est, explicit$scale_est)
  expect_gt(as.numeric(logLik(screened$opt_mod)),
            as.numeric(logLik(single$opt_mod)))
})

test_that("smooth kernels keep the single default start", {
  ex <- running_example("gaussian")
  fit <- multiScale_optim(ex$mod, ex$inputs, n_cores = 1, verbose = FALSE)
  single <- multiScale_optim(ex$mod, ex$inputs, n_cores = 1, verbose = FALSE,
                             start_strategy = "single")
  expect_null(diagnostics(fit)$start_stall)
  expect_equal(fit$scale_est, single$scale_est)
  expect_equal(fit$optim_results$counts, single$optim_results$counts)
})

test_that("a start where every refit fails warns with the refit remedy", {
  ex <- running_example("gaussian")
  # Stand-in for a glmmTMB model with a non-positive-definite Hessian, whose
  # log-likelihood is NA: logLik.glm() returns NA when the AIC is missing.
  no_loglik <- function(model, data, context) {
    m <- stats::update(model, data = data)
    m$aic <- NA_real_
    m
  }
  expect_warning(
    fit <- multiScale_optim(ex$mod, ex$inputs, n_cores = 1, verbose = FALSE,
                            refit_fn = no_loglik),
    "no refitted model gave a finite log-likelihood.*explicit `start` values"
  )
  stall <- diagnostics(fit)$start_stall
  expect_identical(stall$reason, "refit_failed")
  expect_equal(fit$scale_est$Mean, unname(stall$start))
})

test_that("fits that move from their start are unchanged and not flagged", {
  ex <- running_example("gaussian")
  expect_no_warning(
    fit <- multiScale_optim(ex$mod, ex$inputs, n_cores = 1, verbose = FALSE)
  )
  expect_null(diagnostics(fit)$start_stall)
  expect_gt(fit$optim_results$counts[[1L]], 1L)
  expect_output(print(summary(fit)), "Optimized Scale of Effect|Fitted Model")
})

test_that("a user start on a bound is not flagged as a stall", {
  # A zero projected gradient at a bound is a legitimate stop.
  diag <- multiScaleR:::.msr_start_stall_diagnostic(
    opt_results = list(par = c(0.5, 0.01), value = 10, counts = c(`function` = 1L, gradient = 1L)),
    start = c(0.5, 0.01), lower = c(0.01, 0.01), upper = c(1, 1),
    unit_conv = 1000, covariates = c("a", "b"), hard_only = TRUE, screened = FALSE)
  expect_null(diag)
  moved <- multiScaleR:::.msr_start_stall_diagnostic(
    opt_results = list(par = c(0.5, 0.4), value = 10, counts = c(`function` = 1L, gradient = 1L)),
    start = c(0.5, 0.3), lower = c(0.01, 0.01), upper = c(1, 1),
    unit_conv = 1000, covariates = c("a", "b"), hard_only = TRUE, screened = FALSE)
  expect_null(moved)
})

test_that("a default start returned after a failed line search is flagged", {
  # Seen with an irregular refit objective: three evaluations, no step taken.
  stalled <- list(par = c(0.244), value = 524,
                  counts = c(`function` = 3L, gradient = 3L))
  diag <- multiScaleR:::.msr_start_stall_diagnostic(
    opt_results = stalled, start = 0.244, lower = 0.01, upper = 1,
    unit_conv = 1700, covariates = "land1", hard_only = FALSE,
    screened = FALSE, default_start = TRUE)
  expect_identical(diag$reason, "no_progress")
  expect_equal(unname(diag$start), 0.244 * 1700)
  expect_match(diag$message, "found no improvement")

  # A user-supplied or screened start may already be an optimum.
  expect_null(multiScaleR:::.msr_start_stall_diagnostic(
    opt_results = stalled, start = 0.244, lower = 0.01, upper = 1,
    unit_conv = 1700, covariates = "land1", hard_only = FALSE,
    screened = FALSE, default_start = FALSE))
})
