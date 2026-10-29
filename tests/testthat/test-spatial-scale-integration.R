test_that("latent Gaussian integration retains normalization and an uninformative scale limit", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site=8,n_year=3)
  angle <- 2*pi*(as.integer(f$data$site)-1)/8
  f$data$east <- cos(angle);f$data$north <- sin(angle)
  object <- multiScaleR:::.msr_spatial_scale_prepare(count~habitat+temp,f$inputs,f$data,
    list(habitat=~1+(1|year)),family="gaussian",row_id="key",covariate="habitat",
    coordinates=c("east","north"),knots=matrix(c(-.5,0,.5,0,.5,-.5),3,2),range=.7,
    target="marginal_likelihood",reference_scale=200)
  object$parameters$beta <- c(.2,-.1,0)
  object$parameters$log_scale_sd[] <- log(c(.3,.5))
  obj <- multiScaleR:::.msr_spatial_scale_adfun(object,integrate_latent=TRUE)
  expected <- sum(multiScaleR:::.msr_marginal_response_nll(object$specification,object$parameters))
  expect_equal(as.numeric(obj$fn(obj$par)),expected,tolerance=1e-8)
  gradient <- as.numeric(obj$gr(obj$par))
  expect_equal(gradient[names(obj$par)=="log_scale_sd"],c(0,0),tolerance=1e-8)
  # The zero response coefficient removes likelihood information about the
  # ecological scale field. Flat outer SD curvature is retained, not repaired.
  hi <- obj$par;lo <- obj$par
  at <- which(names(hi)=="log_scale_sd")
  hi[at] <- hi[at]+.5;lo[at] <- lo[at]-.5
  expect_equal(obj$fn(hi),obj$fn(lo),tolerance=1e-8)
})
