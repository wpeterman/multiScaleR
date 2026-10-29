spatial_scale_fixture <- function(kernel="gaussian",family="gaussian",target="posterior") {
  f <- bayesian_scale_fixture(n_site=8,n_year=3,kernel=kernel)
  angle <- 2*pi*(as.integer(f$data$site)-1)/8
  f$data$east <- cos(angle); f$data$north <- sin(angle)
  if(family=="binomial") f$data$count <- as.integer(f$data$count>2)
  args <- list(formula=count~habitat+temp,inputs=f$inputs,data=f$data,
    scales=list(habitat=~1+(1|year)),family=family,row_id="key",covariate="habitat",
    coordinates=c("east","north"),knots=matrix(c(-.5,0,.5,0,.5,-.5),3,2),range=.7,target=target)
  if(target=="posterior") {args$field_sd_prior <- .4;args$priors <- list(scale_median=200)} else args$reference_scale <- 200
  list(f=f,args=args,object=do.call(multiScaleR:::.msr_spatial_scale_prepare,args))
}

test_that("spatial basis is centered, normalized and frozen for prediction", {
  x <- spatial_scale_fixture(); out <- x$object; spec <- out$specification
  term <- spec$scales$habitat$random$.spatial_field; B <- term$gaussian_block$design
  expect_equal(colMeans(B),rep(0,ncol(B)),tolerance=1e-12)
  expect_equal(mean(rowSums(B^2)),1,tolerance=1e-12)
  expect_identical(spec$priors$scale_group_sd,c(.5,.4))
  expect_equal(length(spec$parameters$scale_z),6)
  p <- out$parameters; p$scale_z[] <- c(-.2,0,.3,.4,-.5,.1)
  p$log_scale_sd[] <- log(c(.3,.5))
  all <- multiScaleR:::.msr_bayes_components(spec,p)$scales$habitat
  predicted <- multiScaleR:::.msr_spatial_scale_predict(out,x$f$data,parameters=p)
  expect_equal(predicted,all,tolerance=1e-12)
  reverse <- rev(seq_len(nrow(x$f$data)))
  expect_equal(multiScaleR:::.msr_spatial_scale_predict(out,x$f$data[reverse,],parameters=p),all[reverse],tolerance=1e-12)
  subset <- x$f$data[c(1,5,12),]; subset$east <- subset$east + .15
  joint <- multiScaleR:::.msr_spatial_scale_predict(out,subset,parameters=p)
  separate <- vapply(seq_len(nrow(subset)),function(i) multiScaleR:::.msr_spatial_scale_predict(out,subset[i,,drop=FALSE],parameters=p),0)
  expect_equal(unname(joint),unname(separate),tolerance=1e-12)
  expect_false(inherits(out,"multiScaleR_bayes"))
  expect_error(posterior_scales(out),"multiScaleR_bayes")
  expect_error(multiScaleR:::.msr_bayes_prediction_scales(spec,x$f$data,"error",TRUE),"dedicated private")
  for(pair in list(c("posterior","marginal_laplace"),c("marginal_likelihood","factorized"),c("marginal_likelihood","laplace")))
    expect_error(multiScaleR:::.msr_scale_check_capability(spec,pair[1],pair[2]),"Gaussian projection")
})

test_that("spatial prediction propagates joint vectors and coherent new annual effects", {
  x <- spatial_scale_fixture(); out <- x$object
  out$draws <- lapply(out$parameters,function(v) matrix(rep(v,each=4),nrow=4))
  out$draws$scale_z[,4:6] <- matrix(seq(-.3,.4,length.out=12),4)
  predicted <- multiScaleR:::.msr_spatial_scale_predict(out,x$f$data,scale="log")
  expect_equal(dim(predicted),c(4,nrow(x$f$data)))
  for(i in 1:4) expect_equal(predicted[i,],multiScaleR:::.msr_spatial_scale_predict(out,x$f$data,
    parameters=multiScaleR:::.msr_bayes_draw(out,i),scale="log"),tolerance=1e-12)
  future <- x$f$data[1:2,]; future$year <- "new"
  expect_error(multiScaleR:::.msr_spatial_scale_predict(out,future),"new")
  p <- out$parameters; p$scale_z[4:6] <- 0
  future$east[] <- 0;future$north[] <- 0
  set.seed(123)
  sampled <- multiScaleR:::.msr_spatial_scale_predict(out,future,parameters=p,new_levels="sample")
  expect_identical(unname(sampled[1]),unname(sampled[2]))
  args <- x$args;args$field_sd_prior <- NULL
  expect_error(do.call(multiScaleR:::.msr_spatial_scale_prepare,args),"field_sd_prior")
  args <- x$args;args$knots <- rbind(args$knots,args$knots[1,])
  expect_error(do.call(multiScaleR:::.msr_spatial_scale_prepare,args),"distinct")
  args <- x$args;args$scales <- list(habitat=~1+temp+(1|year))
  expect_error(do.call(multiScaleR:::.msr_spatial_scale_prepare,args),"intercept-only fixed scale")
})

test_that("spatial target adapters reject mislabeled population priors", {
  object <- spatial_scale_fixture()$object
  bad <- object;bad$target <- "unknown"
  expect_error(multiScaleR:::.msr_spatial_scale_objective(bad),"Unknown spatial prototype target")
  bad <- object;bad$target <- "marginal_likelihood"
  expect_error(multiScaleR:::.msr_spatial_scale_objective(bad),"population-prior specification")
  bad <- object;bad$specification$priors <- NULL
  expect_error(multiScaleR:::.msr_spatial_scale_objective(bad),"population-prior specification")
})

test_that("spatial targets preserve population priors and derivatives for both kernels and families", {
  skip_if_not_installed("RTMB")
  for(kernel in c("gaussian","exp")) for(family in c("gaussian","poisson","binomial","nbinom2")) {
    b <- spatial_scale_fixture(kernel,family,"posterior")$object
    m <- spatial_scale_fixture(kernel,family,"marginal_likelihood")$object
    p <- b$parameters; p$beta <- c(.2,-.1,.5); p$scale_z[] <- seq(-.2,.3,length.out=6)
    p$log_scale_sd[] <- log(c(.3,.4))
    sm <- b$specification$scales$habitat; pr <- b$specification$priors
    penalty <- -sum(dnorm(p$beta,0,pr$beta_sd,log=TRUE)) -
      sum(dnorm(p$scale_beta,sm$prior_mean,sm$prior_sd,log=TRUE)) -
      sum(dnorm(exp(p$log_scale_sd),0,pr$scale_group_sd,log=TRUE)+log(2)+p$log_scale_sd) -
      sum(dnorm(p$log_dispersion,pr$dispersion_log_mean,pr$dispersion_log_sd,log=TRUE))
    expect_equal(multiScaleR:::.msr_spatial_scale_objective(b)(p)-multiScaleR:::.msr_spatial_scale_objective(m)(p),penalty,tolerance=1e-10)
    m$parameters <- p
    ad <- multiScaleR:::.msr_spatial_scale_adfun(m)
    ordinary <- multiScaleR:::.msr_spatial_scale_objective(m)
    expect_equal(as.numeric(ad$fn(ad$par)),ordinary(p),tolerance=1e-10)
    step <- 1e-5
    gradient <- vapply(seq_along(ad$par),function(j){hi<-lo<-ad$par;hi[j]<-hi[j]+step;lo[j]<-lo[j]-step;(ad$fn(hi)-ad$fn(lo))/(2*step)},0)
    expect_equal(as.numeric(ad$gr(ad$par)),gradient,tolerance=2e-5)
    expect_equal(unname(ad$he(ad$par)),unname(optimHess(ad$par,ad$fn,ad$gr)),tolerance=2e-5)
    p$beta[3] <- 0
    p2 <- p;p2$scale_z[4:6] <- c(3,-2,1)
    expect_equal(multiScaleR:::.msr_marginal_response_nll(m$specification,p),multiScaleR:::.msr_marginal_response_nll(m$specification,p2),tolerance=1e-12)
  }
})
