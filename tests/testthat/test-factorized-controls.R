test_that("optimizing subdivisions preserve factorized derivatives across families", {
  skip_if_not_installed("RTMB")
  for(kernel in c("gaussian","exp"))for(family in c("gaussian","poisson","binomial","nbinom2")) {
    f<-bayesian_scale_fixture(n_site=4,n_year=3,kernel=kernel,seed=8321)
    if(family=="binomial")f$data$count<-as.integer(f$data$count>2)
    s<-multiScaleR:::.msr_marginal_prepare(count~habitat+temp,f$inputs,f$data,
      list(habitat=~1+(1|year)),family,"key",200)
    p<-s$parameters;p$beta<-c(.5,.1,.4);p$log_scale_sd[]<-log(.2)
    if(length(p$log_dispersion))p$log_dispersion[]<-if(family=="gaussian")0 else log(3)
    a<-multiScaleR:::.msr_factorized_adfun(s,p,step=1)
    b<-multiScaleR:::.msr_factorized_adfun(s,p,step=4)
    expect_lt(abs(a$fn(a$par)-b$fn(b$par)),1e-6)
    expect_lt(max(abs(a$gr(a$par)-b$gr(b$par))),1e-4)
    H<-a$he(a$par)
    expect_lt(max(abs(H-b$he(b$par)))/max(1,max(abs(H))),1e-5)
  }
})

test_that("factorized subdivision control is explicit and preserves defaults", {
  skip_if_not_installed("RTMB")
  f<-bayesian_scale_fixture(n_site=8,n_year=3,annual=FALSE,seed=8321)
  s<-multiScaleR:::.msr_marginal_prepare(count~habitat+temp,f$inputs,f$data,
    list(habitat=~1+(1|year)),"gaussian","key",200)
  p<-s$parameters;p$beta<-c(.5,.2,2);p$scale_z[]<-0
  s$y<-multiScaleR:::.msr_bayes_components(s,p)$eta+.15*sin(seq_along(s$y))
  s<-multiScaleR:::.msr_marginal_reduce(s,"scale")
  before<-.Random.seed
  a<-multiScaleR:::.msr_factorized_fit(s,list(starts=1L))
  b<-multiScaleR:::.msr_factorized_fit(s,list(starts=1L,integration_step=4L))
  expect_identical(.Random.seed,before)
  expect_identical(a$control$integration_step,1L)
  expect_identical(b$control$integration_step,4L)
  expect_true(a$diagnostics$passed&&b$diagnostics$passed)
  expect_equal(a$log_likelihood,b$log_likelihood,tolerance=1e-10)
  expect_equal(a$annual,b$annual,tolerance=1e-10)
  expect_identical(b$attempts[[1L]]$integration_step,4L)
  expect_false(b$uncertainty$available)
  for(bad in c(0,2.5,3,16))expect_error(
    multiScaleR:::.msr_factorized_fit(s,list(integration_step=bad)),"integration_step")
  expect_error(multiScaleR:::.msr_factorized_fit(s,list(optimizer_engine="ad")),"control")
})
