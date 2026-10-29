test_that("Gaussian response intercept integration equals a normalized covariance density", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site=4,n_year=3,seed=8321)
  f$data <- f$data[-c(2,7),,drop=FALSE]
  s <- multiScaleR:::.msr_marginal_prepare(count~habitat+temp+(1|year)+offset(log(effort)),
    f$inputs,f$data,list(habitat=~1+(1|year)),"gaussian","key",200)
  p <- s$parameters;p$beta <- c(.2,-.1,.4);p$scale_z[] <- c(-.2,.1,.3)
  p$log_scale_sd[] <- log(.2);p$log_response_sd[] <- log(.35)
  p$log_dispersion[] <- log(.7);p$response_z <- numeric()
  reference <- function(spec,q) {
    eta <- multiScaleR:::.msr_bayes_components(spec,q,include_response=FALSE)$eta
    Z <- model.matrix(~0+factor(spec$response$random[[1]]$index))
    covariance <- diag(exp(2*q$log_dispersion),length(spec$y))+
      exp(2*q$log_response_sd)*tcrossprod(Z)
    residual <- spec$y-eta
    length(residual)/2*log(2*pi)+sum(log(diag(chol(covariance))))+
      drop(crossprod(residual,solve(covariance,residual)))/2
  }
  value <- multiScaleR:::.msr_factorized_response_nll(s,p)
  obj <- RTMB::MakeADFun(function(q)multiScaleR:::.msr_factorized_response_nll(s,q,TRUE),p,silent=TRUE)
  expect_equal(value,reference(s,p),tolerance=1e-10)
  expect_equal(unname(obj$fn(obj$par)),reference(s,p),tolerance=1e-10)
  h <- 1e-5
  fd <- vapply(seq_along(obj$par),function(j) {
    hi <- lo <- obj$par;hi[j] <- hi[j]+h;lo[j] <- lo[j]-h
    (reference(s,obj$env$parList(x=hi))-reference(s,obj$env$parList(x=lo)))/(2*h)
  },0)
  expect_lt(max(abs(as.numeric(obj$gr(obj$par))-fd)),2e-5)
  expect_true(all(is.finite(obj$he(obj$par))))
  small <- p;small$log_response_sd[] <- log(1e-8)
  reduced <- multiScaleR:::.msr_marginal_reduce(s,"response")
  rp <- p;rp$log_response_sd <- rp$response_z <- numeric()
  expect_equal(multiScaleR:::.msr_factorized_response_nll(s,small),
    sum(multiScaleR:::.msr_factorized_response_nll(reduced,rp)),tolerance=1e-7)
  # Subsetting remaps the response levels, even with unbalanced observations.
  pieces <- vapply(1:3,function(i) {
    g <- multiScaleR:::.msr_factorized_group(s,i);q <- p;q$scale_z <- p$scale_z[i]
    multiScaleR:::.msr_factorized_response_nll(g,q)
  },0)
  expect_equal(sum(pieces),value,tolerance=1e-10)
})

test_that("Gaussian nested response partitions preserve independent annual integration", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site=4,n_year=3,seed=8321)
  f$data$year_region <- interaction(f$data$year,rep(c("a","a","b","b"),3))
  s <- multiScaleR:::.msr_marginal_prepare(count~habitat+temp+(1|year_region),f$inputs,f$data,
    list(habitat=~1+(1|year)),"gaussian","key",200)
  p <- s$parameters;p$beta <- c(.2,-.1,.4);p$log_dispersion[] <- log(.7)
  p$log_response_sd[] <- log(.35);p$log_scale_sd[] <- log(.2)
  obj <- multiScaleR:::.msr_factorized_adfun(s,p)
  ordinary <- multiScaleR:::.msr_factorized_objective(s,FALSE,1e-10,.5)
  expect_lt(abs(obj$fn(obj$par)-ordinary(obj$env$parList(x=obj$par))),1e-7)
  expect_true(all(is.finite(obj$gr(obj$par))))
  expect_false(any(names(obj$par)%in%c("scale_z","response_z")))
  # Each response subgroup must be wholly inside one integrated scale group.
  coupled <- multiScaleR:::.msr_marginal_prepare(count~habitat+temp+(1|site),f$inputs,f$data,
    list(habitat=~1+(1|year)),"gaussian","key",200)
  expect_error(multiScaleR:::.msr_factorized_adfun(coupled),"coupling")
  shared <- multiScaleR:::.msr_marginal_reduce(coupled,"scale")
  expect_silent(multiScaleR:::.msr_factorized_validate(shared))
  shared_obj <- multiScaleR:::.msr_factorized_adfun(shared)
  expect_true(all(is.finite(c(shared_obj$fn(shared_obj$par),shared_obj$gr(shared_obj$par)))))
  count <- multiScaleR:::.msr_marginal_prepare(count~habitat+temp+(1|year),f$inputs,f$data,
    list(habitat=~1+(1|year)),"nbinom2","key",200)
  expect_error(multiScaleR:::.msr_factorized_adfun(count),"Gaussian")
})

test_that("known Gaussian variance nonidentification is retained before fitting", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site=4,n_year=3)
  s <- multiScaleR:::.msr_marginal_prepare(count~habitat+temp+(1|key),f$inputs,f$data,
    list(habitat=~1+(1|year)),"gaussian","key",200)
  # Fixed scales remove the coupling restriction and expose variance confounding.
  s <- multiScaleR:::.msr_marginal_reduce(s,"scale")
  fit <- multiScaleR:::.msr_factorized_fit(s,list(starts=1L))
  expect_false(fit$diagnostics$passed)
  expect_match(fit$diagnostics$reason,"singleton")
  expect_length(fit$attempts,0L)
  f$data$effort[] <- 1
  f$data$count <- rep(1:3,each=4)+log(f$data$effort)
  s <- multiScaleR:::.msr_marginal_prepare(count~habitat+temp+(1|year)+offset(log(effort)),
    f$inputs,f$data,list(habitat=~1+(1|year)),"gaussian","key",200)
  fit <- multiScaleR:::.msr_factorized_fit(s,list(starts=1L))
  expect_false(fit$diagnostics$passed)
  expect_match(fit$diagnostics$reason,"unbounded")
  expect_length(fit$attempts,0L)
})
