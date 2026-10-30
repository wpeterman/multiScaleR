bayesian_sq_fixture <- function(kernel = "gaussian", only = FALSE, n_site = 4) {
  f <- bayesian_scale_fixture(n_site = n_site, n_year = 3, annual = FALSE, kernel = kernel)
  b <- f$inputs$binned
  f$inputs$d_list <- rep(list(b$bin_rep), length(b$point_ids))
  f$inputs$raw_cov <- lapply(seq_along(b$point_ids), function(i)
    data.frame(elev = sin(b$bin_rep * (i + 1)) + i*b$bin_rep^2))
  f$inputs$raw_cov[[1]]$elev[] <- 2
  f$inputs$raw_cov[[2]]$elev[2] <- NA_real_
  f$inputs$kernel_dat <- data.frame(row.names = b$point_ids)
  f$inputs$scale_vars <- if (only) msr_vars(rough = surface_var("elev", "sq", weighted = TRUE)) else
    msr_vars(habitat = kernel_var("habitat"), rough = surface_var("elev", "sq", weighted = TRUE))
  if (only) f$inputs$binned <- NULL
  covs <- if (only) "rough" else c("habitat", "rough")
  f$formula <- as.formula(paste("count ~", paste(covs, collapse = "+"),
    "+ temp + (1|site) + offset(log(effort))"))
  f$scales <- setNames(rep(list(~1+(1|year)), length(covs)), covs)
  f
}

test_that("Bayesian SQ preserves fixed reference normalization and proper targets", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp")) for (only in c(FALSE, TRUE))
    for (family in c("gaussian", "poisson", "nbinom2", "binomial")) {
      f <- bayesian_sq_fixture(kernel, only)
      if (family == "binomial") f$data$count <- rep(0:1, 6)
      fit <- multiScale_bayes(f$formula, f$inputs, f$data, f$scales, family,
        row_id = "key", method = "prepare", priors = list(scale_median = 250))
      s <- fit$specification; p <- s$parameters
      expect_identical(fit$inference_metadata$target, "posterior")
      expect_equal(s$center, colMeans(.msr_eval_kernel_engine(.msr_prepare_surface_engine(
        f$inputs, s$covariates, s$rows, FALSE), rep(log(.25), length(s$covariates)))$raw))
      expect_identical(s$priors, .msr_bayes_priors(list(scale_median = 250), s$covariates, f$inputs))
      p$beta[] <- .2; p$scale_z[] <- seq_along(p$scale_z)/10
      p$response_z[] <- .1
      z <- .msr_bayes_components(s, p)
      raw <- vapply(seq_along(s$surface$buffers$rough), function(i) {
        b <- s$surface$buffers$rough[[i]]; d <- b$d
        sigma <- z$scales$rough[i]/s$unit_conv
        w <- if (kernel == "gaussian") exp(-(d^2-min(d)^2)/(2*sigma^2)) else exp(-(d-min(d))/sigma)
        v <- b$v-b$v[1]; sqrt(sum(w*(v-sum(w*v)/sum(w))^2)/sum(w))
      }, 0)
      expect_equal(z$kernel$rough, (raw-s$center["rough"])/s$spread["rough"], ignore_attr = TRUE)
      obj <- RTMB::MakeADFun(.msr_bayes_objective(s, TRUE), p, silent = TRUE)
      expect_equal(as.numeric(obj$fn(obj$par)), .msr_bayes_objective(s)(p), tolerance = 1e-9)
      ll <- switch(family, gaussian = dnorm(s$y,z$eta,exp(p$log_dispersion),log=TRUE),
        poisson = dpois(s$y,exp(z$eta),log=TRUE),
        binomial = dbinom(s$y,s$trials,plogis(z$eta),log=TRUE),
        nbinom2 = dnbinom(s$y,mu=exp(z$eta),size=exp(p$log_dispersion),log=TRUE))
      pr <- s$priors
      prior <- sum(dnorm(p$beta,0,pr$beta_sd,log=TRUE)) +
        sum(dnorm(p$scale_z,log=TRUE)) + sum(dnorm(p$response_z,log=TRUE)) +
        sum(dnorm(exp(p$log_scale_sd),0,pr$scale_group_sd,log=TRUE)+log(2)+p$log_scale_sd) +
        sum(dnorm(exp(p$log_response_sd),0,pr$response_group_sd,log=TRUE)+log(2)+p$log_response_sd) +
        sum(dnorm(p$log_dispersion,pr$dispersion_log_mean,pr$dispersion_log_sd,log=TRUE))
      for (sm in s$scales) prior <- prior + sum(dnorm(p$scale_beta[sm$beta_index],
        sm$prior_mean,sm$prior_sd,log=TRUE))
      expect_equal(.msr_bayes_objective(s)(p), -sum(ll)-prior, tolerance = 1e-10)
      score <- vapply(seq_along(obj$par), function(j) {
        a <- b <- obj$par; a[j] <- a[j]+1e-5; b[j] <- b[j]-1e-5
        (obj$fn(a)-obj$fn(b))/2e-5
      }, 0)
      expect_equal(as.numeric(obj$gr(obj$par)), score, tolerance = 2e-5)
      expect_true(all(is.finite(obj$he(obj$par))))
      expect_equal(s$surface$memory$observation_cells,
        sum(vapply(s$surface$buffers, function(x) sum(vapply(x, function(b) length(b$v), 0L)), 0)))
      expect_error(.msr_scale_check_capability(s, "marginal_likelihood", "factorized"), "Weighted SQ")
    }
})

test_that("SQ public Laplace and marginal fits retain their posterior diagnostics", {
  skip_if_not_installed("RTMB")
  f <- bayesian_sq_fixture(only = TRUE, n_site = 8)
  form <- count ~ rough + temp
  prep <- multiScale_bayes(form,f$inputs,f$data,list(rough=~1+(1|year)),
    "gaussian",row_id="key",method="prepare")
  p <- prep$parameters; p$beta <- c(.5,.3,1.5); p$scale_z <- c(-.4,0,.5)
  p$log_scale_sd <- log(.3)
  f$data$count <- .msr_bayes_components(prep$specification,p)$eta + .3*sin(1:24)
  state <- .Random.seed
  for (method in c("laplace","marginal")) {
    fit <- multiScale_bayes(form,f$inputs,f$data,list(rough=~1+(1|year)),
      "gaussian",row_id="key",method=method,laplace_draws=8,seed=913)
    expect_identical(.Random.seed,state)
    expect_null(fit$sampler)
    expect_identical(fit$inference_metadata$target,"posterior")
    expect_identical(fit$inference_metadata$uncertainty$intervals,"marginal_credible")
    expect_equal(nrow(predict(fit,draws=3)),24)
    expect_equal(dim(simulate(fit,nsim=2,seed=11)),c(24,2))
    expect_true(all(is.finite(as.matrix(posterior_scales(fit)[c("mean","median")]))))
    if(method=="laplace") expect_lt(fit$diagnostics$gradient_max,1e-4) else
      expect_true(fit$diagnostics$grid_converged)
  }
})

test_that("SQ Gaussian response integration agrees with a normalized covariance density", {
  skip_if_not_installed("RTMB")
  f <- bayesian_sq_fixture(only = TRUE)
  fit <- multiScale_bayes(f$formula, f$inputs, f$data, f$scales,
    "gaussian", row_id = "key", method = "prepare")
  s <- fit$specification; p <- s$parameters; p$beta[] <- .3
  eta <- .msr_bayes_components(s,p,include_response=FALSE)$eta
  group <- s$response$random[[1]]$index
  covariance <- diag(exp(2*p$log_dispersion),length(eta)) +
    outer(group,group,"==")*exp(2*p$log_response_sd)
  ch <- chol(covariance); residual <- s$y-eta
  ref <- .5*(length(eta)*log(2*pi)+2*sum(log(diag(ch)))+
    sum(backsolve(ch,residual,transpose=TRUE)^2))
  expect_equal(.msr_factorized_response_nll(s,p),ref,tolerance=1e-10)
  obj <- RTMB::MakeADFun(function(p) .msr_factorized_response_nll(s,p,TRUE),p,silent=TRUE)
  expect_equal(as.numeric(obj$fn(obj$par)),ref,tolerance=1e-10)
  expect_true(all(is.finite(obj$he(obj$par))))
})

test_that("SQ predictions and proposals share frozen reference maps and raw sources", {
  for (only in c(FALSE, TRUE)) {
    f <- bayesian_sq_fixture(only = only)
    fit <- fake_bayesian_draws(multiScale_bayes(f$formula, f$inputs, f$data,
      f$scales, "gaussian", row_id = "key", method = "prepare"), 5)
    fit$draws$beta[,] <- .3
    before <- fit$specification
    old <- predict(fit, summary = FALSE)
    ids <- c(8, 2, 2, 5)
    new <- predict(fit, f$data[ids, ], kernel_inputs = f$inputs, summary = FALSE)
    expect_equal(as.numeric(new), as.numeric(old[,ids]), tolerance = 1e-12)
    expect_error(prepare_prediction(fit), "Preparation requires")
    restored <- unserialize(serialize(fit, NULL))
    expect_equal(predict(restored, f$data[ids, ], kernel_inputs = f$inputs, summary = FALSE), new)
    expect_identical(fit$specification, before)
    moved <- f$inputs; order <- rev(seq_along(moved$d_list))
    moved$d_list <- moved$d_list[order]; moved$raw_cov <- moved$raw_cov[order]
    moved$kernel_dat <- moved$kernel_dat[order,,drop=FALSE]
    if (!is.null(moved$binned)) {
      moved$binned$point_ids <- moved$binned$point_ids[order]
      moved$binned$vsum <- lapply(moved$binned$vsum,function(x)x[order,,drop=FALSE])
      moved$binned$csum <- lapply(moved$binned$csum,function(x)x[order,,drop=FALSE])
    }
    expect_equal(predict(fit,f$data[ids,],kernel_inputs=moved,summary=FALSE),new)
    units <- f$inputs; units$unit_conv <- units$unit_conv/2
    units$d_list <- lapply(units$d_list, function(d) d*2)
    if (!is.null(units$binned)) units$binned$bin_rep <- units$binned$bin_rep*2
    expect_equal(predict(fit, f$data[ids, ], kernel_inputs = units, summary = FALSE), new)
    changed <- f$inputs; changed$scale_vars$source[changed$scale_vars$type == "surface"] <- "wrong"
    expect_error(predict(fit, f$data, kernel_inputs = changed), "specifications")
    missing <- f$inputs; missing$raw_cov <- NULL
    expect_error(predict(fit, f$data, kernel_inputs = missing), "raw buffers")
    blocks <- fit$draws
    expect_equal(unname(.msr_bayes_batch_eta(before, blocks)), unname(old))
    unseen <- f$data[ids,]; unseen$year <- factor(rep("new",length(ids)))
    expect_error(predict(fit,unseen,kernel_inputs=f$inputs),"New scale group")
    tab <- predict(fit,unseen,kernel_inputs=f$inputs,new_scale_levels="zero")
    expect_true(all(tab$new_scale_level))
    expect_equal(dim(simulate(fit,nsim=2,seed=51,newdata=unseen,
      kernel_inputs=f$inputs,new_scale_levels="zero")),c(length(ids),2))
    group <- .msr_factorized_group(before, 1)
    expect_length(group$surface$buffers$rough, sum(before$scales[[1]]$random[[1]]$index == 1))
  }
})

test_that("SQ differentiation fails before tape construction when RAM is unresolved", {
  f <- bayesian_sq_fixture(only = TRUE)
  fit <- multiScale_bayes(f$formula, f$inputs, f$data, f$scales,
    row_id = "key", method = "prepare")
  local_mocked_bindings(.msr_total_ram_bytes = function() NA_real_, .package = "multiScaleR")
  expect_error(.msr_bayes_objective(fit$specification, TRUE), "Cannot detect RAM")
  local_mocked_bindings(.msr_total_ram_bytes = function() 1024^3, .package = "multiScaleR")
  expect_error(.msr_bayes_objective(fit$specification, TRUE), "RAM budget")
})

test_that("SQ reference preparation retains finite nearest ties at tiny prior medians", {
  skip_if_not_installed("RTMB")
  f <- bayesian_sq_fixture(only=TRUE)
  for(i in 2:4) {
    f$inputs$d_list[[i]][2] <- f$inputs$d_list[[i]][1]
    f$inputs$raw_cov[[i]]$elev[1:2] <- c(0,i)
  }
  fit <- multiScale_bayes(f$formula,f$inputs,f$data,f$scales,"gaussian",
    row_id="key",priors=list(scale_median=1e-300),method="prepare")
  s <- fit$specification
  expected <- c(0,1,1.5,2)[s$rows]
  expect_equal(s$center["rough"],mean(expected),ignore_attr=TRUE)
  expect_equal(s$spread["rough"],sd(expected),ignore_attr=TRUE)
  obj <- RTMB::MakeADFun(.msr_bayes_objective(s,TRUE),s$parameters,silent=TRUE)
  expect_true(all(is.finite(c(obj$fn(obj$par),obj$gr(obj$par),obj$he(obj$par)))))
})

test_that("SQ raw evaluation retains broad, nearest and constant-buffer limits", {
  skip_if_not_installed("RTMB")
  buffers <- list(list(d = c(1, 1, 2), v = c(2, 4, 10)),
    list(d = c(1, 2, 3), v = c(2, 4, 8)), list(d = 1:3, v = c(5,5,5)))
  for (kernel in c("gaussian", "exp")) {
    invalid <- buffers; invalid[[1]]$d[1] <- -1
    expect_error(.msr_scale_surface_raw(invalid,rep(0,3),kernel),"nonnegative")
    expect_equal(.msr_scale_surface_raw(buffers, rep(1000,3), kernel),
      vapply(buffers, function(b) sqrt(mean((b$v-mean(b$v))^2)), 0))
    expect_identical(.msr_scale_surface_raw(buffers, rep(-1000,3), kernel), c(1,0,0))
    for (value in c(-1000, -4, 1, 1000)) {
      shifted <- lapply(buffers, function(b) { b$v <- b$v+2^40; b })
      expect_equal(.msr_scale_surface_raw(shifted,rep(value,3),kernel),
        .msr_scale_surface_raw(buffers,rep(value,3),kernel),tolerance=1e-14)
      obj <- RTMB::MakeADFun(function(p) sum(.msr_scale_surface_raw(buffers,
        rep(p$x,3), kernel, TRUE)), list(x = value), silent = TRUE)
      expect_equal(as.numeric(obj$fn(obj$par)), sum(.msr_scale_surface_raw(buffers, rep(value,3), kernel)))
      expect_true(all(is.finite(c(obj$gr(obj$par), obj$he(obj$par)))))
    }
  }
})
