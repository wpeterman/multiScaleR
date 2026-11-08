# Existing bounded software fixtures only. No recovery/sampler/seed search.
test_that("fixed spline matrices follow native crossings, columns and row layouts", {
  dat <- data.frame(y = 1:12, h = c(-1, 0, .2, .5, .8, 1, 2, .2-1e-8,
    .2+1e-8, .5-1e-8, .5+1e-8, .3), k = seq(.05, .95, length.out = 12),
    z = seq(.3, .9, length.out = 12), group = factor(rep(c("a", "b", "c"), 4)))
  forms <- list(
    y ~ splines::ns(h, knots = c(.2,.5), Boundary.knots = c(0,1))*group + z,
    y ~ 0 + splines::ns(h/2, knots = c(.1,.25), Boundary.knots = c(0,.5), intercept = TRUE):group + k,
    y ~ splines::bs(h, knots = c(.2,.5), Boundary.knots = c(0,1), warn.outside = FALSE)*I(k^2) + z,
    y ~ splines::ns(h, knots = c(.2,.5), Boundary.knots = c(0,1)):
      splines::bs(k, knots = c(.3,.7), Boundary.knots = c(0,1)) + group,
    y ~ splines::ns(h, knots = NULL, Boundary.knots = c(0,1)) + I(h^2):group,
    y ~ splines::bs(h, knots = NULL, Boundary.knots = c(0,1), intercept = TRUE, warn.outside = FALSE) + z)
  contrasts(dat$group) <- contr.sum(3)
  for (form in forms) {
    tt <- terms(form); frame <- model.frame(tt, dat)
    design <- .msr_compile_response_design(tt, frame, c("h", "k"))
    kernel <- as.matrix(dat[c("h", "k")]); X <- model.matrix(tt, frame)
    beta <- seq(-.2, .2, length.out = ncol(X))
    expect_equal(.msr_response_design_matrix(design, kernel), X, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(.msr_response_design_eta(design, kernel, beta), drop(X %*% beta), tolerance = 1e-12, ignore_attr = TRUE)
    sensitivity <- .msr_response_kernel_sensitivity(design, kernel, beta)
    for (j in seq_len(ncol(kernel))) {
      hi <- lo <- kernel;hi[,j] <- hi[,j]+1e-6;lo[,j] <- lo[,j]-1e-6
      numerical <- drop((.msr_response_design_matrix(design, hi)-
        .msr_response_design_matrix(design, lo)) %*% beta)/2e-6
      expect_lt(max(abs(sensitivity[,j]-numerical)), 1e-6)
    }
    next_kernel <- kernel[c(5,1,8,8,12,4,11,3,7,2,9,6),,drop=FALSE]
    updated <- .msr_response_update_frame(design, frame, next_kernel)
    nd <- dat;nd$h <- next_kernel[,1];nd$k <- next_kernel[,2]
    expect_equal(model.matrix(tt, updated), model.matrix(tt, model.frame(tt, nd)),
      tolerance = 1e-12, ignore_attr = TRUE)
    for (nm in design$expression_frame_columns) if (is.matrix(frame[[nm]]))
      expect_identical(attributes(updated[[nm]]), attributes(frame[[nm]]))
    selected <- unique(c(ncol(X), 1L, min(3L,ncol(X))))
    subset <- .msr_compile_response_design(tt, frame, c("h","k"), columns = colnames(X)[selected])
    expect_equal(.msr_response_design_matrix(subset, kernel), X[,selected,drop=FALSE], tolerance = 1e-12, ignore_attr = TRUE)
    rows <- c(5L,1L,5L)
    sliced <- .msr_response_design_rows(design, rows)
    expect_equal(.msr_response_design_matrix(sliced, kernel[rows,,drop=FALSE]), X[rows,,drop=FALSE], tolerance = 1e-12, ignore_attr = TRUE)
    for (program in design$expression_programs) if (!is.null(program$basis_recipe)) {
      expect_setequal(ls(environment(program$value)),c("r","f"))
      expect_identical(parent.env(environment(program$value)),environment(.msr_spline_columns))
      expect_true(all(vapply(program$derivative,function(d)
        setequal(ls(environment(d)),c("r","f","d")),logical(1))))
    }
    saved <- unserialize(serialize(design, NULL, version = 2))
    expect_equal(.msr_response_design_matrix(saved, kernel), X, tolerance = 1e-12, ignore_attr = TRUE)
  }
})

test_that("fixed spline AD tapes retain scores and curvature when crossing knots", {
  skip_if_not_installed("RTMB")
  original_config <- RTMB::TapeConfig()
  dat <- data.frame(y = 1:7, h = c(-.1,0,.2,.5,.8,1,1.2),
    k = c(.1,.6,.4,.7,.2,.5,.9), group = factor(c("a","b","a","b","a","b","a")))
  for (type in c("ns", "bs")) {
    form <- as.formula(paste0("y ~ splines::",type,
      "(h, knots=c(.2,.5), Boundary.knots=c(0,1)",
      if(type=="bs") ", warn.outside=FALSE" else "", ")*group + I(h*k)"))
    tt <- terms(form);frame <- model.frame(tt,dat)
    design <- .msr_compile_response_design(tt,frame,c("h","k"))
    kernel <- as.matrix(dat[c("h","k")]);beta <- seq(-.2,.2,length.out=design$p)
    obj <- RTMB::MakeADFun(function(p) sum(.msr_response_design_eta(design,
      p$kernel,p$beta)^2),list(kernel=kernel,beta=beta),silent=TRUE)
    numerical <- vapply(seq_along(obj$par),function(j) {
      hi <- lo <- obj$par;hi[j] <- hi[j]+1e-5;lo[j] <- lo[j]-1e-5
      (obj$fn(hi)-obj$fn(lo))/2e-5
    },0)
    expect_lt(max(abs(as.numeric(obj$gr(obj$par))-numerical)),2e-5)
    H <- vapply(seq_along(obj$par),function(j) {
      # Exact knots have continuous curvature but a jump in the next derivative.
      # Cancel the O(h) central-difference error rather than loosen the gate.
      difference <- function(h) {
        hi <- lo <- obj$par;hi[j] <- hi[j]+h;lo[j] <- lo[j]-h
        (as.numeric(obj$gr(hi))-as.numeric(obj$gr(lo)))/(2*h)
      }
      2*difference(5e-6)-difference(1e-5)
    },numeric(length(obj$par)))
    expect_lt(max(abs(obj$he(obj$par)-H)),1e-4)
    expect_gt(max(abs(obj$he(obj$par))),.01)
    next_dat <- dat;next_dat$h <- rev(dat$h);at <- obj$par
    at[seq_along(kernel)] <- as.vector(as.matrix(next_dat[c("h","k")]))
    expected <- drop(model.matrix(form,next_dat)%*%beta)
    expect_equal(as.numeric(obj$fn(at)),sum(expected^2),tolerance=1e-10)
    expect_identical(RTMB::TapeConfig(),original_config)
  }
})

test_that("fixed spline profiles retain native family targets and kernel derivatives", {
  for (family in c("gaussian","poisson","binomial","nbinom2")) {
    f <- formula_integration_fixture(family,two=TRUE)
    response <- formula(f$model)[[2L]]
    form <- as.formula(call("~",response,quote(
      splines::ns(hab,knots=c(-.4,.6),Boundary.knots=c(-5,5))*group +
      splines::bs(cover,knots=c(-.5,.5),Boundary.knots=c(-5,5)):hab +
      z + offset(log(exposure)))))
    f$model <- glm(form,f$data,family=f$model$family)
    contract <- .msr_profile_contract(f$model,f$inputs);obj <- contract$objective
    p <- c(log(c(.3,.4)),if(family=="nbinom2") log(2))
    result <- obj$evaluate(p);nd <- f$data
    kd <- .msr_eval_kernel_engine(obj$engine,p[1:2])
    nd$hab <- kd$value[,"hab"];nd$cover <- kd$value[,"cover"]
    native <- glm(form,nd,family=f$model$family,control=glm.control(epsilon=1e-11,maxit=100))
    expect_lt(abs(result$value+as.numeric(logLik(native))),1e-7)
    rebuilt <- .msr_profile_conditional_model(contract,result,f$model,p[1:2])
    expect_equal(model.matrix(rebuilt),model.matrix(native),tolerance=1e-10,ignore_attr=TRUE)
    numerical <- vapply(seq_along(p),function(j) {
      hi <- lo <- p;hi[j] <- hi[j]+1e-5;lo[j] <- lo[j]-1e-5
      (obj$fn(hi)-obj$fn(lo))/2e-5
    },0)
    expect_lt(max(abs(as.numeric(obj$gr(p))-numerical)),2e-5)
    joint <- .msr_profile_joint(obj,p)
    expect_true(is.finite(joint$fn(joint$par)))
  }
})

test_that("fixed spline prepared weighted-SQ targets preserve priors and prediction recipes", {
  skip_if_not_installed("RTMB")
  for (surface in c(FALSE,TRUE)) for (family in c("gaussian","poisson","binomial","nbinom2")) {
    f <- prepared_formula_fixture(surface=surface)
    if(family=="binomial") f$data$count <- as.integer(f$data$count>2)
    # Eight independent sites cannot identify two rich bases simultaneously.
    # Keep this fixture/seed and use one nonlinear source, including SQ itself.
    form <- if(surface)
      count ~ splines::ns(other,knots=c(-.4,.6),Boundary.knots=c(-5,5))*year + habitat + poly(temp,2) + offset(log(effort)) else
      count ~ splines::ns(habitat,knots=c(-.4,.6),Boundary.knots=c(-5,5))*year + other + poly(temp,2) + offset(log(effort))
    fit <- multiScale_bayes(form,f$inputs,f$data,f$scales,family,row_id="key",method="prepare",
      priors=list(scale_median=250))
    s <- fit$specification;p <- s$parameters
    p$beta[] <- seq(-.2,.2,length.out=length(p$beta));p$scale_z[] <- seq(-.1,.2,length.out=length(p$scale_z))
    native <- prepared_formula_native_eta(s,p,f$data)
    expect_equal(.msr_bayes_components(s,p)$eta,native$eta,tolerance=1e-12,ignore_attr=TRUE)
    expect_identical(s$beta_names,colnames(native$X))
    expect_identical(s$priors,.msr_bayes_priors(list(scale_median=250),s$covariates,f$inputs))
    eta <- native$eta
    likelihood <- switch(family,gaussian=dnorm(s$y,eta,exp(p$log_dispersion),log=TRUE),
      poisson=dpois(s$y,exp(eta),log=TRUE),binomial=dbinom(s$y,s$trials,plogis(eta),log=TRUE),
      nbinom2=dnbinom(s$y,mu=exp(eta),size=exp(p$log_dispersion),log=TRUE))
    expect_equal(.msr_bayes_objective(s)(p),-sum(likelihood)-prepared_formula_prior(s,p),tolerance=1e-10)
    obj <- RTMB::MakeADFun(.msr_bayes_objective(s,TRUE),p,silent=TRUE)
    expect_equal(as.numeric(obj$fn(obj$par)),.msr_bayes_objective(s)(p),tolerance=1e-10)
    numerical <- vapply(seq_along(obj$par),function(j) {
      hi <- lo <- obj$par;hi[j] <- hi[j]+1e-5;lo[j] <- lo[j]-1e-5
      (obj$fn(hi)-obj$fn(lo))/2e-5
    },0)
    expect_lt(max(abs(as.numeric(obj$gr(obj$par))-numerical)),2e-5)
    if(surface) expect_error(.msr_scale_check_capability(s,"marginal_likelihood","factorized"),"Weighted SQ")
    if(family=="poisson") {
      fit$parameters <- p;fit <- fake_bayesian_draws(fit,n=4)
      original <- predict(fit,summary=FALSE);ids <- c(7,1,13,7)
      projected <- predict(fit,newdata=f$data[ids,],kernel_inputs=f$inputs,summary=FALSE)
      expect_equal(as.numeric(projected),as.numeric(original[,ids,drop=FALSE]),tolerance=1e-12)
      single <- predict(fit,newdata=f$data[13,,drop=FALSE],kernel_inputs=f$inputs,summary=FALSE)
      expect_equal(as.numeric(single),as.numeric(original[,13,drop=FALSE]),tolerance=1e-12)
      grDevices::pdf(tempfile(fileext=".pdf"))
      plots <- tryCatch(plot_marginal_effects(fit,length.out=5),finally=grDevices::dev.off())
      expect_setequal(names(plots),c("habitat","other","year","temp"))
      expect_true(all(vapply(plots,function(p) all(is.finite(p$data$fit)),logical(1))))
      expect_identical(fit$specification[c("center","spread")],s[c("center","spread")])
    }
  }
})

test_that("unidentifiable combined spline preparation retains its full-rank guard", {
  f <- prepared_formula_fixture()
  form <- count ~ splines::ns(habitat,knots=c(-.4,.6),Boundary.knots=c(-5,5))*year +
    splines::bs(other,knots=c(-.5,.5),Boundary.knots=c(-5,5)) + poly(temp,2) + offset(log(effort))
  expect_error(multiScale_bayes(form,f$inputs,f$data,f$scales,"poisson",row_id="key",method="prepare"),"full-rank")
})

test_that("weighted-SQ spline profile scores retain scale and shape derivatives", {
  for(kernel in c("gaussian","exp")) {
    f <- formula_weighted_fixture(kernel=kernel)
    form <- y ~ splines::bs(rough,knots=c(-.4,.6),Boundary.knots=c(-5,5),warn.outside=FALSE)*group + hab + z + offset(log(exposure))
    model <- glm(form,f$data,family=poisson())
    contract <- .msr_profile_contract(model,f$inputs);obj <- contract$objective
    p <- c(log(c(.3,.4)),if(kernel=="expow") log(1.8))
    kd <- .msr_eval_kernel_engine(obj$engine,p);nd <- f$data
    nd$hab <- kd$value[,"hab"];nd$rough <- kd$value[,"rough"]
    native <- glm(form,nd,family=poisson(),control=glm.control(epsilon=1e-11,maxit=100))
    expect_lt(abs(obj$fn(p)+as.numeric(logLik(native))),1e-7)
    numerical <- vapply(seq_along(p),function(j) {
      hi <- lo <- p;hi[j] <- hi[j]+1e-5;lo[j] <- lo[j]-1e-5
      (obj$fn(hi)-obj$fn(lo))/2e-5
    },0)
    expect_lt(max(abs(as.numeric(obj$gr(p))-numerical)),2e-5)
    expect_setequal(contract$objective$response_design$covariates,c("rough","hab"))
  }
})

test_that("dynamic spline policy rejects adaptive and ambiguous recipes", {
  validate <- function(form) .msr_validate_kernel_expressions(terms(form),"h")
  expect_error(validate(y~splines::ns(h,df=3)),"explicit fixed")
  expect_error(validate(y~splines::ns(h,knots=c(.2,.5))),"explicit fixed")
  expect_error(validate(y~splines::ns(h,df=3,knots=c(.2,.5),Boundary.knots=c(0,1))),"adaptive df")
  expect_error(validate(y~splines::ns(h,knots=c(.2,.2),Boundary.knots=c(0,1))),"distinct")
  expect_error(validate(y~splines::ns(h,knots=c(0,.5),Boundary.knots=c(0,1))),"strictly inside")
  expect_error(validate(y~splines::ns(h,knots=c(.2,.5),Boundary.knots=c(1,0))),"increasing")
  expect_error(validate(y~splines::bs(h,knots=c(.2,.5),Boundary.knots=c(0,1),degree=2)),"cubic")
  adaptive <- c(.2,.5)
  expect_error(validate(y~splines::ns(h,knots=adaptive,Boundary.knots=c(0,1))),"literal fixed")
  env <- new.env(parent=baseenv());env$ns <- function(...) 1
  form <- as.formula("y~ns(h,knots=c(.2,.5),Boundary.knots=c(0,1))",env=env)
  expect_error(validate(form),"original splines")
  expect_no_error(validate(y~splines::ns(h,knots=c(.2,.5),Boundary.knots=c(0,1))))
  program <- .msr_kernel_expression(quote(splines::ns(h,knots=c(.2,.5),Boundary.knots=c(0,1))),"h",environment())
  f <- formula_integration_fixture()
  engine <- .msr_prepare_kernel_engine(f$inputs$binned,"hab")
  normalized <- lapply(log(c(.15,.55)),function(kp) {
    kd <- .msr_eval_kernel_engine(engine,kp)
    raw_knots <- kd$center[1]+kd$spread[1]*c(.2,.5)
    raw_bounds <- kd$center[1]+kd$spread[1]*c(0,1)
    native <- splines::ns(kd$raw[,1],knots=raw_knots,Boundary.knots=raw_bounds)
    list(modeled=do.call(cbind,program$value(list(kd$value[,1]))),
      raw=matrix(as.numeric(native),nrow=nrow(kd$raw)),physical_knots=raw_knots)
  })
  # A fixed knot in candidate-standardized units maps to changing raw cover.
  # The native raw-unit oracle transforms knots/boundaries by that same map.
  expect_equal(lapply(normalized,`[[`,"modeled"),lapply(normalized,`[[`,"raw"),tolerance=1e-12)
  expect_gt(max(abs(normalized[[1]]$physical_knots-normalized[[2]]$physical_knots)),1e-8)
  cat("FIXED_SPLINE_NATIVE_AD recipes=fixed_modeled_input guards=retained\n")
})


test_that("public fixed spline fits preserve coefficient SEs, K and frozen curves", {
  f <- formula_integration_fixture()
  model <- glm(y ~ splines::ns(hab,knots=c(-.4,.6),Boundary.knots=c(-5,5)) +
    group + z + offset(log(exposure)),f$data,family=poisson())
  fit <- multiScale_optim(model,f$inputs,backend="profile_glm",n_cores=1,verbose=FALSE,
    profile_control=list(starts=matrix(c(.08,.25,.6),ncol=1)))
  expect_equal(nrow(fit$profile_fit$start_results),3L)
  expect_true(any(fit$profile_fit$start_results$eligible))
  expect_equal(unname(predict(fit)),unname(predict(fit$opt_mod,type="response")),tolerance=1e-8)
  conditional <- predict(fit,type="link",se.fit=TRUE,uncertainty="conditional")
  native <- predict(fit$opt_mod,type="link",se.fit=TRUE)
  expect_equal(conditional$fit,native$fit,tolerance=1e-8)
  expect_equal(conditional$se.fit,native$se.fit,tolerance=1e-6)
  expect_equal(.msr_selection_parameter_count(fit),as.numeric(attr(logLik(fit$opt_mod),"df"))+1L)
  expect_identical(predict(prepare_prediction(fit)),predict(fit))
  order <- rev(seq_len(nrow(f$data)))
  expect_equal(unname(predict(fit,f$data[order,],f$inputs)),rev(unname(predict(fit))),tolerance=1e-10)
  grDevices::pdf(tempfile(fileext=".pdf"))
  plots <- tryCatch(suppressMessages(plot_marginal_effects(fit,length.out=5)),finally=grDevices::dev.off())
  expect_setequal(names(plots),c("hab","group","z"))
  expect_true(all(vapply(plots,function(p) all(is.finite(p$data$fit)),logical(1))))
})

test_that("fixed splines retain native mixed count and choice-set targets", {
  skip_if_not_installed("RTMB");skip_if_not_installed("glmmTMB");skip_if_not_installed("survival")
  f <- readRDS(test_path("fixtures","joint-mixed-existing.rds"))
  dat <- f$data;dat$site <- factor(dat$site);dat$year <- factor(dat$year)
  engine <- .msr_prepare_kernel_engine(f$inputs$binned,"habitat");kp <- log(.23)
  dat$habitat <- .msr_eval_kernel_engine(engine,kp)$value[,1L]
  for(family in c("poisson","nbinom2")) {
    model <- glmmTMB::glmmTMB(count ~ splines::ns(habitat,knots=c(-.4,.6),Boundary.knots=c(-5,5)) +
      year + poly(temp,2) + (1|site) + offset(log(effort)),dat,
      family=if(family=="poisson") poisson() else glmmTMB::nbinom2())
    args <- .msr_native_mixed_design(model,engine,"habitat",kp,rownames(dat))
    obj <- do.call(.msr_mixed_rtmb,c(args,list(fixed_scale=TRUE)))
    expect_lt(abs(obj$fn(obj$par)+as.numeric(logLik(model))),1e-5)
    at <- obj$par+.02
    numerical <- vapply(seq_along(at),function(j) {
      hi <- lo <- at;hi[j] <- hi[j]+1e-5;lo[j] <- lo[j]-1e-5
      (obj$fn(hi)-obj$fn(lo))/2e-5
    },0)
    expect_lt(max(abs(as.numeric(obj$gr(at))-numerical)),2e-5)
  }
  f <- profile_glm_fixture();dat <- f$data
  dat$case_ <- rep(c(1,0,0),length.out=nrow(dat));dat$choice <- rep(seq_len(nrow(dat)/3),each=3)
  form <- as.formula("case_ ~ splines::ns(hab,knots=c(-.4,.6),Boundary.knots=c(-5,5)) + group + z + offset(log(exposure)) + strata(choice)",env=asNamespace("survival"))
  coxph <- survival::coxph
  model <- survival::clogit(form,dat,model=TRUE,x=TRUE,y=TRUE)
  contract <- .msr_clogit_contract(model,f$inputs);obj <- .msr_clogit_objective(contract)
  native <- .msr_clogit_native(contract,log(.3))
  at <- c(log(.3),unname(coef(native)[contract$original_columns]))
  expect_lt(abs(obj$fn(at)+as.numeric(logLik(native))),1e-10)
  numerical <- vapply(seq_along(at),function(j) {
    hi <- lo <- at;hi[j] <- hi[j]+1e-5;lo[j] <- lo[j]-1e-5
    (obj$fn(hi)-obj$fn(lo))/2e-5
  },0)
  expect_lt(max(abs(as.numeric(obj$gr(at))-numerical)),2e-5)
})


test_that("fixed spline private groups and endpoints preserve native Gaussian targets", {
  f <- prepared_formula_fixture(two=FALSE)
  form <- count ~ splines::ns(habitat,knots=c(-.4,.6),Boundary.knots=c(-5,5))*year +
    poly(temp,2) + offset(log(effort))
  s <- .msr_marginal_prepare(form,f$inputs,f$data,f$scales,"gaussian",row_id="key",reference_scale=250)
  p <- s$parameters;p$beta[] <- seq(-.2,.2,length.out=length(p$beta))
  p$scale_z[] <- seq(-.1,.1,length.out=length(p$scale_z))
  eta <- prepared_formula_native_eta(s,p,f$data)$eta
  expected <- -dnorm(s$y,eta,exp(p$log_dispersion),log=TRUE)
  expect_equal(.msr_factorized_response_nll(s,p),expected,tolerance=1e-10,ignore_attr=TRUE)
  term <- s$scales[[1]]$random[[1]]
  for(j in seq_along(term$levels)) {
    group <- .msr_factorized_group(s,j);ids <- which(term$index==j)
    q <- p;q$scale_z <- p$scale_z[j]
    expect_equal(.msr_factorized_response_nll(group,q),expected[ids],tolerance=1e-10,ignore_attr=TRUE)
  }
  for(kind in c("zero","broad","nearest")) {
    endpoint <- if(kind=="zero") .msr_marginal_reduce(s,"scale") else .msr_factorized_kernel_limit(s,kind)
    q <- p;q$scale_z <- q$log_scale_sd <- numeric()
    if(kind!="zero") q$scale_beta <- numeric()
    reference <- prepared_formula_native_eta(endpoint,q,f$data)$eta
    expect_equal(sum(.msr_factorized_response_nll(endpoint,q)),
      -sum(dnorm(endpoint$y,reference,exp(q$log_dispersion),log=TRUE)),tolerance=1e-10)
  }
})

test_that("fixed spline Gaussian coefficient proposals retain normalized prior factors", {
  skip_if_not_installed("RTMB")
  f <- prepared_formula_fixture(two=FALSE)
  fit <- multiScale_bayes(count ~ splines::ns(habitat,knots=c(-.4,.6),Boundary.knots=c(-5,5))*year +
    poly(temp,2) + offset(log(effort)),f$inputs,f$data,f$scales,"gaussian",row_id="key",method="prepare")
  s <- fit$specification;p <- s$parameters
  obj <- RTMB::MakeADFun(.msr_bayes_objective(s,TRUE),p,silent=TRUE)
  keep <- which(names(obj$par)!="beta");outer <- matrix(rep(obj$par[keep],each=3),nrow=3)
  previous_seed <- get0(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  set.seed(5811)
  candidate <- tryCatch(.msr_bayes_beta_candidates(outer,keep,s,p$beta),finally={
    if(is.null(previous_seed)) rm(".Random.seed",envir=.GlobalEnv) else
      assign(".Random.seed",previous_seed,envir=.GlobalEnv)
  })
  difference <- vapply(seq_len(3),function(i) -obj$fn(candidate$samples[i,]),0)-candidate$logq
  expect_lt(diff(range(difference)),1e-8)
  expect_equal(candidate$fallbacks,0)
})
