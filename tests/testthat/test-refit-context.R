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
