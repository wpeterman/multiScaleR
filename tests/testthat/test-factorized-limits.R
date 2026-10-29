test_that("kernel limit likelihoods retain normalized densities and offsets", {
  skip_if_not_installed("RTMB")
  for(kernel in c("gaussian","exp"))for(family in c("gaussian","poisson","binomial","nbinom2")) {
    f<-bayesian_scale_fixture(n_site=4,n_year=3,kernel=kernel,seed=8321)
    if(family=="binomial")f$data$count<-as.integer(f$data$count>2)
    s<-multiScaleR:::.msr_marginal_prepare(count~habitat+temp,f$inputs,f$data,
      list(habitat=~1+(1|year)),family,"key",200)
    s$offset<-log(f$data$effort)
    if(family=="binomial")s$trials[]<-10
    p<-s$parameters;p$beta<-c(.5,.1,.4)
    if(length(p$log_dispersion))p$log_dispersion[]<-if(family=="gaussian")0 else log(3)
    for(limit in c("broad","nearest"))for(poisson in if(family=="nbinom2")c(FALSE,TRUE) else FALSE) {
      g<-multiScaleR:::.msr_factorized_limit_spec(s,limit,poisson)
      q<-p;q$scale_beta<-q$scale_z<-q$log_scale_sd<-numeric()
      if(poisson)q$log_dispersion<-numeric()
      keep<-if(limit=="nearest")s$delta[[1L]]==0 else s$delta[[1L]]>=0
      k<-rowSums(s$V[[1L]]*keep)/rowSums(s$C[[1L]]*keep)
      k<-(k-s$center[[1L]])/s$spread[[1L]]
      eta<-drop(s$X_fixed%*%p$beta[seq_len(ncol(s$X_fixed))])+s$offset+
        p$beta[ncol(s$X_fixed)+1L]*k
      ref<--sum(switch(g$family,
        gaussian=dnorm(s$y,eta,exp(q$log_dispersion),log=TRUE),
        poisson=dpois(s$y,exp(eta),log=TRUE),
        binomial=dbinom(s$y,s$trials,plogis(eta),log=TRUE),
        nbinom2=dnbinom(s$y,mu=exp(eta),size=exp(q$log_dispersion),log=TRUE)))
      obj<-multiScaleR:::.msr_factorized_adfun(g,q)
      expect_lt(abs(obj$fn(obj$par)-ref),1e-8)
      expect_lt(abs(multiScaleR:::.msr_factorized_objective(g,FALSE)(q)-ref),1e-8)
      expect_length(g$parameters$scale_beta,0L)
      expect_length(g$scales[[1L]]$random,0L)
      if(g$family=="poisson") {
        native<-glm.fit(cbind(s$X_fixed,k),s$y,offset=s$offset,family=poisson())
        q$beta<-unname(native$coefficients)
        target<--sum(dpois(s$y,native$fitted.values,log=TRUE))
        expect_lt(abs(multiScaleR:::.msr_factorized_objective(g,FALSE)(q)-target),1e-8)
      }
    }
  }
})

test_that("Poisson and kernel limits reject a finite-scale identification claim", {
  skip_if_not_installed("RTMB")
  f<-bayesian_scale_fixture(n_site=6,n_year=3,seed=8321)
  s<-multiScaleR:::.msr_marginal_prepare(count~habitat+temp,f$inputs,f$data,
    list(habitat=~1+(1|year)),"nbinom2","key",200)
  s<-multiScaleR:::.msr_marginal_reduce(s,"scale")
  g<-multiScaleR:::.msr_factorized_limit_spec(s,"broad",TRUE)
  k<-(rowSums(s$V[[1L]])/rowSums(s$C[[1L]])-s$center[[1L]])/s$spread[[1L]]
  native<-glm.fit(cbind(s$X_fixed,k),s$y,offset=s$offset,family=poisson())
  p<-s$parameters;p$beta<-unname(native$coefficients);p$scale_beta[]<-20
  p$log_dispersion[]<-log(1e10)
  candidate<-multiScaleR:::.msr_factorized_objective(s,FALSE)(p)
  audit<-multiScaleR:::.msr_factorized_limit_audit(s,p,candidate)
  expect_true(audit$completed)
  expect_true("broad-poisson"%in%audit$witnesses)
  expect_true("fitted-scales-poisson"%in%audit$witnesses)
  expect_length(audit$references$`broad-poisson`$attempts,3L)
  expect_error(multiScaleR:::.msr_factorized_limit_spec(s,"invalid"),"kernel limit")
})
test_that("fixed-kernel witnesses also support scales without an annual SD", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 6, n_year = 2)
  for (scale_formula in list(~1, ~0 + year)) {
    s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp + offset(log(effort)),
      f$inputs, f$data, list(habitat = scale_formula), "gaussian", "key", 200)
    p <- s$parameters; p$beta <- c(.2, -.1, .4)
    for (kernel in c("broad", "nearest")) {
      limited <- multiScaleR:::.msr_factorized_limit_spec(s, kernel)
      keep <- if (kernel == "nearest") s$delta[[1L]] == 0 else
        matrix(TRUE, nrow(s$V[[1L]]), ncol(s$V[[1L]]))
      weighted <- rowSums(s$V[[1L]] * keep) / rowSums(s$C[[1L]] * keep)
      eta <- drop(s$X_fixed %*% p$beta[seq_len(ncol(s$X_fixed))]) + s$offset +
        p$beta[ncol(s$X_fixed) + 1L] * (weighted - s$center[[1L]]) / s$spread[[1L]]
      expected <- -sum(dnorm(s$y, eta, exp(p$log_dispersion), log = TRUE))
      p$scale_beta <- p$scale_z <- p$log_scale_sd <- numeric()
      expect_equal(multiScaleR:::.msr_factorized_objective(limited, FALSE)(p), expected,
        tolerance = 1e-10)
    }
  }
})

