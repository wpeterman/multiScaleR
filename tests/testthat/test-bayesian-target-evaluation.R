testthat::test_that("importance targets preserve every finite AD value across families and kernels", {
  testthat::skip_if_not_installed("RTMB")
  for(family in c("gaussian","poisson","binomial","nbinom2"))for(kernel in c("gaussian","exp")) {
    f<-bayesian_scale_fixture(n_site=4,n_year=3,kernel=kernel)
    if(family=="binomial")f$data$count<-as.numeric(seq_len(nrow(f$data))%%2)
    prepared<-multiScaleR::multiScale_bayes(count~habitat+temp+(1|site),f$inputs,f$data,
      scales=list(habitat=~1+(1|year)),family=family,row_id="key",method="prepare")
    s<-prepared$specification
    obj<-RTMB::MakeADFun(multiScaleR:::.msr_bayes_objective(s,TRUE),s$parameters,silent=TRUE)
    samples<-rbind(obj$par,obj$par+.02*sin(seq_along(obj$par)),obj$par-.03*cos(seq_along(obj$par)))
    original<-vapply(seq_len(nrow(samples)),function(i)-obj$fn(samples[i,]),0)
    value<-multiScaleR:::.msr_bayes_importance_target(samples,obj,s)
    testthat::expect_identical(value$log_density,original)
    testthat::expect_equal(unname(value$diagnostics["ad_invalid"]),0)
    testthat::expect_equal(unname(value$diagnostics["evaluations"]),3)
  }
})
testthat::test_that("overflowing Poisson means retain valid zero posterior mass", {
  testthat::skip_if_not_installed("RTMB")
  f<-bayesian_scale_fixture(n_site=3,n_year=2)
  f$data$count<-rep(c(0,1),length.out=nrow(f$data))
  fit<-multiScaleR::multiScale_bayes(count~habitat+(1|site)+(1|year),f$inputs,f$data,
    scales=list(habitat=~1+(1|year)),family="poisson",row_id="key",method="prepare")
  s<-fit$specification;p<-s$parameters;p$beta[]<-0;p$response_z[]<-0
  term<-s$response$random[["year"]]
  p$log_response_sd[term$sd_index]<-log(200)
  p$response_z[term$z_index]<-c(-5,5)
  vector<-unlist(p,use.names=FALSE)
  obj<-RTMB::MakeADFun(multiScaleR:::.msr_bayes_objective(s,TRUE),s$parameters,silent=TRUE)
  testthat::expect_true(all(is.finite(vector)))
  testthat::expect_true(is.na(-obj$fn(vector)))
  testthat::expect_identical(-multiScaleR:::.msr_bayes_objective(s)(p),-Inf)
  component<-multiScaleR:::.msr_bayes_components(s,p)
  testthat::expect_true(all(is.finite(component$eta)))
  testthat::expect_true(any(is.infinite(exp(component$eta))))
  value<-multiScaleR:::.msr_bayes_importance_target(matrix(vector,nrow=1),obj,s)
  testthat::expect_identical(value$log_density,-Inf)
  testthat::expect_equal(unname(value$diagnostics[c("ad_invalid","recovered_zero","unresolved")]),c(1,1,0))
  testthat::expect_length(value$log_density,1)
})
testthat::test_that("ordinary rechecks preserve RNG and leave unresolved targets rejected", {
  f<-bayesian_scale_fixture(n_site=3,n_year=2)
  fit<-multiScaleR::multiScale_bayes(count~habitat+(1|site),f$inputs,f$data,
    scales=list(habitat=~1+(1|year)),family="gaussian",row_id="key",method="prepare")
  s<-fit$specification;vector<-unlist(s$parameters,use.names=FALSE);samples<-matrix(vector,nrow=1)
  set.seed(71027);before<-.Random.seed
  value<-multiScaleR:::.msr_bayes_importance_target(samples,list(fn=function(p)NaN),s)
  testthat::expect_identical(.Random.seed,before)
  testthat::expect_equal(value$log_density,-multiScaleR:::.msr_bayes_objective(s)(s$parameters),tolerance=0)
  testthat::expect_equal(unname(value$diagnostics[c("recovered_finite","unresolved")]),c(1,0))
  zero<-multiScaleR:::.msr_bayes_importance_target(samples,list(fn=function(p)Inf),s)
  testthat::expect_identical(zero$log_density,-Inf)
  testthat::expect_equal(unname(zero$diagnostics["ad_invalid"]),0)
  samples[1,1]<-NaN
  bad<-suppressWarnings(multiScaleR:::.msr_bayes_importance_target(samples,list(fn=function(p)NaN),s))
  testthat::expect_true(is.na(bad$log_density))
  testthat::expect_equal(unname(bad$diagnostics["unresolved"]),1)
  testthat::expect_error(multiScaleR:::.msr_bayes_importance_checks(samples,bad$log_density,s),"Invalid posterior importance ratios")
  testthat::expect_error(multiScaleR:::.msr_bayes_importance_checks(samples,Inf,s),"Invalid posterior importance ratios")
})
