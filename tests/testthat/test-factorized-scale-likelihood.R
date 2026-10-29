test_that("factorized kernel densities agree and remain differentiable in extreme tails", {
  skip_if_not_installed("RTMB")
  for(kernel in c("gaussian","exp")) for(family in c("gaussian","poisson","binomial","nbinom2")) {
    f <- bayesian_scale_fixture(n_site=4,n_year=3,kernel=kernel)
    if(family=="binomial")f$data$count <- as.integer(f$data$count>2)
    s <- multiScaleR:::.msr_marginal_prepare(count~habitat+temp,f$inputs,f$data,
      list(habitat=~1+(1|year)),family,"key",200)
    p <- s$parameters;p$beta <- c(.2,-.1,.4);p$scale_z[] <- c(-.2,.1,.3)
    expect_equal(unname(multiScaleR:::.msr_factorized_response_nll(s,p)),
      unname(multiScaleR:::.msr_marginal_response_nll(s,p)),tolerance=1e-10)
    g <- multiScaleR:::.msr_factorized_group(s,2L);one <- p;one$scale_z <- p$scale_z[2]
    rows <- which(s$scales[[1]]$random[[1]]$index==2)
    expect_equal(unname(multiScaleR:::.msr_factorized_response_nll(g,one)),
      unname(multiScaleR:::.msr_factorized_response_nll(s,p)[rows]),tolerance=1e-10)
    for(z in c(-1e6,1e6)) {
      p$scale_z[] <- z
      obj <- RTMB::MakeADFun(function(q)sum(multiScaleR:::.msr_factorized_response_nll(s,q,TRUE)),p,silent=TRUE)
      expect_true(all(is.finite(c(obj$fn(obj$par),obj$gr(obj$par),obj$he(obj$par)))))
    }
  }
})

test_that("ordinary binomial likelihoods retain finite extreme-logit values", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site=4,n_year=3)
  f$data$count <- as.integer(f$data$count>2)
  s <- multiScaleR:::.msr_marginal_prepare(count~habitat+temp,f$inputs,f$data,
    list(habitat=~1+(1|year)),"binomial","key",200)
  s$trials[] <- 10
  p <- s$parameters
  for(eta in c(-1000,1000)) {
    p$beta <- c(eta,0,0)
    value <- sum(multiScaleR:::.msr_factorized_response_nll(s,p))
    obj <- RTMB::MakeADFun(function(q)sum(multiScaleR:::.msr_factorized_response_nll(s,q,TRUE)),p,silent=TRUE)
    expect_true(is.finite(value))
    expect_equal(unname(value),unname(obj$fn(obj$par)),tolerance=1e-10)
  }
})

test_that("concentrated Gaussian integrals retain finite audited derivatives", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site=20,n_year=5,seed=8321)
  s <- multiScaleR:::.msr_marginal_prepare(count~habitat+temp,f$inputs,f$data,
    list(habitat=~1+(1|year)),"gaussian","key",200)
  p <- s$parameters;p$beta <- c(.5,.2,2)
  p$scale_z[] <- c(-1.3,-.6,.1,.7,1.1);p$log_scale_sd[] <- log(.4)
  f$data$count <- multiScaleR:::.msr_bayes_components(s,p)$eta+.15*sin(seq_len(nrow(f$data)))
  s$y <- f$data$count;p$log_dispersion[] <- log(.106)
  p$scale_z <- numeric()
  obj <- multiScaleR:::.msr_factorized_adfun(s,p)
  reference <- function(theta)multiScaleR:::.msr_factorized_objective(s,FALSE,1e-10,.5)(obj$env$parList(x=theta))
  gradient <- as.numeric(obj$gr(obj$par))
  expect_true(all(is.finite(c(obj$fn(obj$par),gradient,obj$he(obj$par)))))
  expect_lt(abs(obj$fn(obj$par)-reference(obj$par)),1e-7)
  h <- 1e-6
  fd <- vapply(seq_along(obj$par),function(j) {
    hi<-lo<-obj$par;hi[j]<-hi[j]+h;lo[j]<-lo[j]-h
    (reference(hi)-reference(lo))/(2*h)
  },0)
  expect_lt(max(abs(gradient-fd)),1e-4)
})

test_that("a second integration coordinate resolves very large annual SDs", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site=4,n_year=3)
  s <- multiScaleR:::.msr_marginal_prepare(count~habitat+temp,f$inputs,f$data,
    list(habitat=~1+(1|year)),"nbinom2","key",200)
  p <- s$parameters;p$beta <- c(.2,-.1,.4)
  p$log_scale_sd[] <- log(1000);p$scale_z <- numeric()
  obj <- RTMB::MakeADFun(multiScaleR:::.msr_factorized_scale_objective(s,TRUE),p,silent=TRUE)
  reference <- function(theta)multiScaleR:::.msr_factorized_scale_objective(s,FALSE)(obj$env$parList(x=theta))
  gradient <- as.numeric(obj$gr(obj$par))
  expect_true(all(is.finite(gradient)))
  expect_lt(abs(obj$fn(obj$par)-reference(obj$par)),1e-7)
  h <- 1e-5
  fd <- vapply(seq_along(obj$par),function(j) {
    hi<-lo<-obj$par;hi[j]<-hi[j]+h;lo[j]<-lo[j]-h
    (reference(hi)-reference(lo))/(2*h)
  },0)
  expect_lt(max(abs(gradient-fd)),2e-5)
})

test_that("adaptive factorized integration has independent values and derivatives", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site=4,n_year=3)
  s <- multiScaleR:::.msr_marginal_prepare(count~habitat+temp,f$inputs,f$data,
    list(habitat=~1+(1|year)),"gaussian","key",200)
  p <- s$parameters;p$beta <- c(.5,.1,.6);p$log_dispersion[] <- log(1)
  p$log_scale_sd[] <- log(.2)
  obj <- multiScaleR:::.msr_factorized_adfun(s,p)
  reference <- function(x) marginal_annual_integral(s,obj$env$parList(x=x),rel.tol=1e-10)
  expect_lt(abs(obj$fn(obj$par)-reference(obj$par)),1e-7)
  h <- 1e-5
  fd <- vapply(seq_along(obj$par),function(j) {
    hi <- lo <- obj$par;hi[j] <- hi[j]+h;lo[j] <- lo[j]-h
    (reference(hi)-reference(lo))/(2*h)
  },0)
  expect_lt(max(abs(as.numeric(obj$gr(obj$par))-fd)),2e-5)
  H <- obj$he(obj$par)
  coarse <- stats::optimHess(obj$par,obj$fn,obj$gr)
  fine <- stats::optimHess(obj$par,obj$fn,obj$gr,
    control=list(ndeps=rep(1e-4,length(obj$par))))
  expect_lt(max(abs(H-fine)),1e-4)
  expect_lt(max(abs(H-fine)),max(abs(H-coarse))/10)
  cold <- multiScaleR:::.msr_factorized_adfun(s,p,1e-10,.5)
  expect_lt(abs(obj$fn(obj$par)-cold$fn(obj$par)),1e-7)
  expect_lt(max(abs(obj$gr(obj$par)-cold$gr(obj$par))),1e-5)
  expect_false(any(names(obj$par)=="scale_z"))
  p$log_scale_sd[] <- log(1e-7)
  near <- multiScaleR:::.msr_factorized_adfun(s,p)
  reduced <- multiScaleR:::.msr_marginal_reduce(s,"scale")
  # The reduced parameter list must also remove the positive-SD entry.
  zp <- p;zp$scale_z <- zp$log_scale_sd <- numeric()
  zero <- multiScaleR:::.msr_factorized_adfun(reduced,zp)
  expect_lt(abs(near$fn(near$par)-zero$fn(zero$par)),1e-7)
})

test_that("factorized fitting preserves exact reduced models and rejects coupled years", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site=12,n_year=4,annual=FALSE)
  f$data$temp <- rep(seq(-1,1,length.out=12),4)
  s <- multiScaleR:::.msr_marginal_prepare(count~habitat+temp,f$inputs,f$data,
    list(habitat=~1+(1|year)),"gaussian","key",200)
  p <- s$parameters;p$beta <- c(.5,.2,2)
  f$data$count <- multiScaleR:::.msr_bayes_components(s,p)$eta+rep(.15*sin(1:12),4)
  s <- multiScaleR:::.msr_marginal_prepare(count~habitat+temp,f$inputs,f$data,
    list(habitat=~1+(1|year)),"gaussian","key",200)
  s <- multiScaleR:::.msr_marginal_reduce(s,"scale")
  before <- .Random.seed
  fit <- multiScaleR:::.msr_factorized_fit(s,list(starts=2L))
  expect_identical(.Random.seed,before)
  expect_true(fit$diagnostics$passed)
  expect_s3_class(fit,"msr_factorized_prototype")
  expect_false(fit$uncertainty$available)
  expect_true(all(!fit$annual$interval_available))
  expect_equal(fit$annual$conditional_mean_log,fit$annual$log_estimate)
  expect_length(fit$attempts,2L)
  expect_true(all(vapply(fit$attempts,`[[`,TRUE,"eligible")))
  coupled <- multiScaleR:::.msr_marginal_prepare(count~habitat+(1|site),f$inputs,f$data,
    list(habitat=~1+(1|year)),"gaussian","key",200)
  expect_error(multiScaleR:::.msr_factorized_adfun(coupled),"coupling")
  expect_error(multiScaleR:::.msr_factorized_fit(s,list(prior_sd=1)),"control")

  # A failed derived summary must not discard the completed optimization.
  testthat::local_mocked_bindings(.msr_factorized_predictions=function(...)stop("forced summary failure"),
    .package="multiScaleR")
  failed_summary <- multiScaleR:::.msr_factorized_fit(s,list(starts=1L))
  expect_false(failed_summary$diagnostics$passed)
  expect_true(failed_summary$diagnostics$fit_passed_before_prediction)
  expect_length(failed_summary$attempts,1L)
  expect_true(length(failed_summary$parameters)>0)
  expect_match(failed_summary$prediction_error,"forced summary failure")

  # Preserve the fitted point if an independent numerical audit itself errors.
  original <- multiScaleR:::.msr_factorized_objective
  testthat::local_mocked_bindings(.msr_factorized_objective=function(spec,ad=FALSE,...) {
    if(!ad)stop("forced ordinary audit failure")
    original(spec,ad,...)
  },.package="multiScaleR")
  failed_audit <- multiScaleR:::.msr_factorized_fit(s,list(starts=1L))
  a <- failed_audit$attempts[[1L]]
  expect_false(failed_audit$diagnostics$passed)
  expect_true(is.finite(a$objective))
  expect_true(length(a$parameters)>0)
  expect_match(a$error,"forced ordinary audit failure")
  expect_identical(a$error_stage,"ordinary integration")
})
