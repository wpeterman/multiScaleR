# Skipped on CRAN to keep the package check within CRAN's time limit. CI sets
# NOT_CRAN and runs this file.
skip_on_cran()

test_that('custom refitters distinguish objectives from returned models', {
  f <- make_unmarked_fixture()
  flags <- logical()
  refit <- function(model,data,context) {
    flags <<- c(flags,context$return_model)
    if(identical(context$return_model,FALSE)) unmarked::update(model,data=data,se=FALSE)
    else unmarked::update(model,data=data)
  }
  ctx <- build_opt_context(f$mod,f$kernel_inputs$raw_cov,refit_fn=refit)
  par <- c(40,60)/f$kernel_inputs$unit_conv
  run <- function(return_model=NULL,context=ctx) kernel_scale_fn(par,
    f$kernel_inputs$d_list,f$kernel_inputs$raw_cov,f$kernel_inputs$kernel,
    f$mod,mod_return=return_model,opt_context=context)
  scalar <- run()
  result <- run(TRUE)
  scalar_again <- run()
  expect_identical(flags,c(FALSE,TRUE,FALSE))
  expect_null(ctx$return_model)
  baseline <- ctx; baseline$refit_fn <- NULL
  reference <- run(TRUE,baseline)
  expect_equal(scalar,-as.numeric(unmarked::logLik(reference$mod)),tolerance=1e-7)
  expect_equal(scalar_again,scalar,tolerance=1e-7)
  expect_equal(unmarked::coef(result$mod),unmarked::coef(reference$mod),tolerance=1e-7)
  expect_equal(unmarked::vcov(result$mod),unmarked::vcov(reference$mod),tolerance=1e-7)
  expect_identical(result$mod@sitesRemoved,reference$mod@sitesRemoved)
  # A model originally fitted without uncertainty must remain that way.
  lean_ctx <- ctx; lean_ctx$fitted_mod <- unmarked::update(f$mod,se=FALSE)
  lean <- run(TRUE,lean_ctx)
  lean_ref <- lean_ctx; lean_ref$refit_fn <- NULL
  expect_error(unmarked::vcov(lean$mod))
  expect_error(unmarked::vcov(run(TRUE,lean_ref)$mod))
})

test_that("default unmarked objective refits skip only the Hessian", {
  f <- make_unmarked_fixture()
  ctx <- build_opt_context(f$mod, f$kernel_inputs$raw_cov)
  # Refit on the fitted frame so the comparison isolates the se setting.
  lean <- .refit_model(f$mod, f$umf, utils::modifyList(ctx, list(return_model = FALSE)))
  full <- .refit_model(f$mod, f$umf, utils::modifyList(ctx, list(return_model = TRUE)))
  expect_false(lean@call$se)
  expect_error(unmarked::vcov(lean))
  expect_true(all(is.finite(unmarked::vcov(full))))
  expect_identical(as.numeric(unmarked::logLik(lean)), as.numeric(unmarked::logLik(full)))
  expect_identical(unmarked::coef(lean), unmarked::coef(full))
  # Objective values through kernel_scale_fn() match the returned model, and
  # the returned model keeps its coefficient covariance.
  par <- c(40, 60) / f$kernel_inputs$unit_conv
  run <- function(return_model = NULL) kernel_scale_fn(par, f$kernel_inputs$d_list,
    f$kernel_inputs$raw_cov, f$kernel_inputs$kernel, f$mod,
    mod_return = return_model, opt_context = ctx)
  returned <- run(TRUE)
  expect_equal(run(), -as.numeric(unmarked::logLik(returned$mod)), tolerance = 1e-10)
  expect_true(all(is.finite(unmarked::vcov(returned$mod))))
  expect_null(returned$mod@call$se)
  # Only unmarked fitting functions with an `se` argument are changed.
  expect_true(.msr_unmarked_accepts_se(f$mod))
  bare <- f$mod; bare@call[[1L]] <- as.name("pcount")
  expect_true(.msr_unmarked_accepts_se(bare))
  wrapped <- f$mod; wrapped@call[[1L]] <- as.name("my_pcount")
  expect_false(.msr_unmarked_accepts_se(wrapped))
  other <- f$mod; other@call[[1L]] <- quote(mypkg::pcount)
  expect_false(.msr_unmarked_accepts_se(other))
  no_se <- f$mod; no_se@call[[1L]] <- quote(unmarked::IDS)
  expect_false(.msr_unmarked_accepts_se(no_se))
})
