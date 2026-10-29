test_that("NB2 continuation resets and audits cold independently", {
  f <- profile_glm_fixture()
  original <- .msr_profile_contract(f$model, f$inputs)$objective
  make <- function(mode) .msr_joint_nb2_objective(original$engine, original$y,
    X_fixed=original$X_fixed, offset=original$offset, weights=original$weights,
    coefficient_start=mode)
  cold <- make("cold"); warm <- make("warm")
  p <- c(log(.25), log1p(.5)); q <- p + c(.002, 0)
  warm$reset_start()
  expect_equal(warm$evaluate(p)$value, cold$evaluate(p)$value, tolerance=1e-10)
  a <- warm$evaluate(q); b <- cold$evaluate(q)
  expect_equal(a$value,b$value,tolerance=1e-10)
  expect_equal(a$gradient,b$gradient,tolerance=1e-7)
  before <- warm$coefficient_audit()$counts
  expect_gt(unname(before['warm']),0)
  audited <- warm$cold_evaluate(q)
  expect_equal(audited$value,b$value,tolerance=1e-10)
  expect_equal(unname(warm$coefficient_audit()$counts['cold']),unname(before['cold'])+1)
  warm$reset_start(); warm$evaluate(q)
  expect_equal(unname(warm$coefficient_audit()$counts['cold']),unname(before['cold'])+2)
  expect_null(cold$reset_start)
  expect_identical(.msr_profile_control(list(),TRUE)$coefficient_start,"cold")
  expect_identical(.msr_profile_control(list(coefficient_start="warm"))$coefficient_start,"warm")
  solver <- .msr_nb2_warm_coefficients()
  solver$state$beta <- rep(NA_real_, ncol(a$design))
  recovered <- solver$solve(a$design,original$y,original$offset,original$weights,
    a$alpha,.msr_nb2_normalizer(original$y,original$weights))
  expect_equal(unname(solver$state$counts['fallback']),1)
  expect_equal(recovered$terms$value,a$value,tolerance=1e-7)
  expect_length(solver$state$failures,1)
  fallback <- .msr_nb2_cold_fallback
  scope <- new.env(parent=environment(fallback))
  scope$.msr_nb2_coefficients <- function(...)stop("injected cold failure")
  environment(fallback) <- scope
  recovered <- fallback(a$design,original$y,original$offset,original$weights,
    a$alpha,.msr_nb2_normalizer(original$y,original$weights))
  expect_equal(recovered$terms$value,a$value,tolerance=1e-7)
  expect_identical(recovered$initialization_audit$reference_error,"injected cold failure")
  expect_lte(recovered$initialization_audit$score,1e-7)
})

test_that("a disagreeing cold audit cannot produce an eligible warm fit", {
  resets <- 0L
  obj <- list(fn=function(x)sum(x^2),gr=function(x)2*x,
    reset_start=function()resets<<-resets+1L,
    cold_evaluate=function(x)list(value=sum(x^2)+1,gradient=2*x))
  z <- .msr_multistart(obj,matrix(c(.2,.4),ncol=1),-1,1)
  expect_null(z$best)
  expect_false(any(z$records$eligible))
  expect_equal(resets,2)
  expect_true(all(grepl("cold coefficient audit",z$records$message)))
})

test_that("warm reduced searches retain starts and cold terminal checks", {
  f <- profile_glm_fixture()
  base <- .msr_profile_contract(f$model,f$inputs)$objective
  make <- function(mode) .msr_joint_nb2_objective(base$engine,base$y,
    X_fixed=base$X_fixed,offset=base$offset,weights=base$weights,coefficient_start=mode)
  starts <- rbind(c(log(.1),0),c(log(.4),log1p(.5)))
  lower <- c(log(.05),0); upper <- c(log(.8),log1p(1e4))
  control <- .msr_profile_control(list(),TRUE)
  cold <- make("cold"); warm <- make("warm")
  a <- .msr_inference_fit(cold,starts,lower,upper,control,fixed_scale=log(.25))
  b <- .msr_inference_fit(warm,starts,lower,upper,control,fixed_scale=log(.25))
  expect_equal(b$value,a$value,tolerance=1e-7)
  expect_identical(nrow(b$records),nrow(a$records))
  expect_gt(unname(warm$coefficient_audit()$counts['warm']),0)
  expect_gte(unname(warm$coefficient_audit()$counts['resets']),4)
  expect_true(all(b$records$eligible))
  expect_true(all(is.finite(b$records$cold_gap) & b$records$cold_gap<=1e-5))
})

test_that("opt-in continuation preserves public prediction and uncertainty", {
  # Reuse the package's existing regular NB2 regression fixture and seed.
  f <- profile_glm_fixture()
  fit <- function(control) suppressWarnings(multiScale_optim(f$model,f$inputs,
    backend="joint_nb2",n_cores=1,verbose=FALSE,profile_control=control))
  cold <- fit(list())
  warm <- fit(list(coefficient_start="warm",kernel_engine="stream"))
  expect_equal(warm$optim_results$value,cold$optim_results$value,tolerance=1e-7)
  expect_identical(warm$profile_fit$starts,cold$profile_fit$starts)
  expect_identical(warm$uncertainty$available,cold$uncertainty$available)
  expect_equal(vcov(warm),vcov(cold),tolerance=1e-5)
  expect_equal(predict(warm,se.fit=TRUE),predict(cold,se.fit=TRUE),tolerance=1e-5)
  expect_equal(predict(warm,newdata=f$data,kernel_inputs=f$inputs,se.fit=TRUE),
    predict(cold,newdata=f$data,kernel_inputs=f$inputs,se.fit=TRUE),tolerance=1e-5)
  expect_true(is.list(warm$profile_fit$coefficient_audit))
  expect_null(cold$profile_fit$coefficient_audit)
})

test_that("finite-theta continuation preserves profile coordinates and cold metadata", {
  f <- profile_glm_fixture()
  a <- .msr_profile_contract(f$model,f$inputs)$objective
  b <- .msr_profile_nb2_continuation(a)
  expect_false(isTRUE(b$boundary_nb2))
  expect_identical(b$parameter_names,a$parameter_names)
  for (theta in c(.1,2,1e5)) {
    p <- c(log(.25),log(theta))
    expect_equal(b$evaluate(p)$value,a$evaluate(p)$value,tolerance=1e-7)
    expect_equal(b$gr(p),a$gr(p),tolerance=1e-5)
    expect_equal(b$cold_evaluate(p)$fit,a$evaluate(p)$fit)
  }
  a$family <- 'gaussian'
  expect_error(.msr_profile_nb2_continuation(a),'NB2')
})

test_that("profile NB2 continuation retains public uncertainty and clone contracts", {
  f <- profile_glm_fixture()
  fit <- function(control) suppressWarnings(multiScale_optim(f$model,f$inputs,
    backend='profile_glm',n_cores=1,verbose=FALSE,profile_control=control))
  a <- fit(list())
  b <- fit(list(coefficient_start='warm',kernel_engine='stream'))
  expect_equal(b$optim_results$value,a$optim_results$value,tolerance=1e-7)
  expect_identical(b$profile_fit$starts,a$profile_fit$starts)
  expect_identical(b$uncertainty$available,a$uncertainty$available)
  expect_equal(vcov(b),vcov(a),tolerance=1e-5)
  expect_equal(predict(b,se.fit=TRUE),predict(a,se.fit=TRUE),tolerance=1e-5)
  expect_equal(predict(b,newdata=f$data,kernel_inputs=f$inputs,se.fit=TRUE),
    predict(a,newdata=f$data,kernel_inputs=f$inputs,se.fit=TRUE),tolerance=1e-5)
  expect_gt(unname(b$profile_fit$coefficient_audit$counts['warm']),0)
  setup <- .msr_inference_setup(b,2L)
  expect_identical(setup$obj$coefficient_start,'warm')
})

test_that("warm objectives remain functional after serialization without output cycles", {
  f <- profile_glm_fixture()
  base <- .msr_profile_contract(f$model,f$inputs)$objective
  joint <- .msr_joint_nb2_objective(base$engine,base$y,base$family,
    base$X_fixed,base$offset,base$weights,coefficient_start='warm')
  profile <- .msr_profile_nb2_continuation(base)
  for(obj in list(joint,profile)) {
    p <- c(log(.25),if(isTRUE(obj$boundary_nb2))log1p(.5) else log(2))
    before <- obj$cold_evaluate(p)
    restored <- unserialize(serialize(obj,NULL))
    expect_false(exists('obj',environment(restored$evaluate),inherits=FALSE))
    restored$reset_start()
    expect_equal(restored$evaluate(p)$value,before$value,tolerance=1e-7)
    expect_equal(restored$cold_evaluate(p)$gradient,before$gradient,tolerance=1e-7)
    expect_true(is.list(restored$coefficient_audit()))
  }
})
