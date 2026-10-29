test_that("fitted prediction and simulation retain their original calculations and RNG", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f, "gaussian"), 7)
  fit$draws$beta[,] <- matrix(seq_len(length(fit$draws$beta)) / 30, 7)
  before <- fit
  reference <- t(vapply(1:7, function(i)
    multiScaleR:::.msr_bayes_components(fit$specification,
      multiScaleR:::.msr_bayes_draw(fit, i))$eta, numeric(nrow(f$data))))
  dimnames(reference) <- NULL
  set.seed(110); state <- .Random.seed
  expect_identical(predict(fit, summary = FALSE), reference)
  expect_identical(.Random.seed, state)
  set.seed(25); selected <- sample(1:7, 5, replace = TRUE)
  expected <- sapply(selected, function(i) {
    p <- multiScaleR:::.msr_bayes_draw(fit, i)
    eta <- multiScaleR:::.msr_bayes_components(fit$specification, p)$eta
    rnorm(length(eta), eta, exp(p$log_dispersion))
  })
  state <- .Random.seed
  actual <- simulate(fit, nsim = 5, seed = 25)
  expect_identical(as.numeric(as.matrix(actual)), as.numeric(expected))
  expect_identical(.Random.seed, state)
  expect_identical(fit, before)
})

test_that("new rows reuse fitted maps, normalization, offsets and explicit IDs", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f, "gaussian"))
  fit$draws$beta[,ncol(fit$draws$beta)] <- 1
  old <- predict(fit, summary = FALSE)
  rows <- c(8, 2, 2, 5)
  nd <- f$data[rows, ]; nd$count <- NULL
  new <- predict(fit, nd, kernel_inputs = f$inputs, row_id = "key", summary = FALSE)
  expect_equal(as.numeric(new), as.numeric(old[,rows]), tolerance = 1e-12)
  tab <- predict(fit, nd, kernel_inputs = f$inputs, row_id = "key")
  expect_identical(tab$point_id, as.character(nd$key))
  expect_false(any(tab$new_scale_level | tab$new_response_level))
  expect_equal(nrow(attr(tab, "new_groups")), 0)
  moved <- f$inputs; order <- rev(seq_along(moved$binned$point_ids))
  moved$binned$point_ids <- moved$binned$point_ids[order]
  moved$binned$vsum <- lapply(moved$binned$vsum, function(x) x[order,,drop=FALSE])
  moved$binned$csum <- lapply(moved$binned$csum, function(x) x[order,,drop=FALSE])
  expect_equal(as.numeric(predict(fit, nd, kernel_inputs = moved, row_id = "key", summary = FALSE)), as.numeric(new))
  units <- f$inputs; units$unit_conv <- units$unit_conv / 10
  units$binned$bin_rep <- units$binned$bin_rep * 10
  expect_equal(as.numeric(predict(fit, nd, kernel_inputs = units, row_id = "key", summary = FALSE)), as.numeric(new), tolerance=1e-12)
  twice <- nd; twice$effort <- nd$effort * 2
  expect_equal(as.numeric(predict(fit, twice, kernel_inputs=f$inputs, row_id="key", type="link", summary=FALSE)-new), rep(log(2),length(new)))
  changed <- nd; changed$count <- c(NA, -100, 1e6, 12)
  expect_identical(predict(fit, changed, kernel_inputs=f$inputs,row_id="key",summary=FALSE), new)
  expect_error(predict(fit, nd), "kernel_inputs")
  expect_error(predict(fit, kernel_inputs=f$inputs), "require newdata")
  bad <- f$inputs; bad$kernel <- "exp"
  expect_error(predict(fit,nd,kernel_inputs=bad,row_id="key"), "fitted kernel")
  bad <- f$inputs
  bad$scale_vars <- data.frame(covariate = "habitat", type = "surface")
  expect_error(predict(fit,nd,kernel_inputs=bad,row_id="key"), "kernel means")
  bad <- f$inputs; bad$binned$point_ids[2] <- bad$binned$point_ids[1]
  expect_error(predict(fit,nd,kernel_inputs=bad,row_id="key"), "unique")
  bad <- nd; bad$key[1] <- "not-a-map"
  expect_error(predict(fit,bad,kernel_inputs=f$inputs,row_id="key"), "do not match")
  expect_error(predict(fit,nd[FALSE,],kernel_inputs=f$inputs,row_id="key"), "nonempty")
})

test_that("saved prediction metadata and inherited IDs support an exponential kernel", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3, kernel = "exp")
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f, "gaussian"), 20)
  fit$draws$beta[, ncol(fit$draws$beta)] <- 1
  nd <- f$data[c(8, 1, 1), ]; nd$count <- NULL
  expected <- predict(fit, summary = FALSE)[, c(8, 1, 1)]
  saved <- unserialize(serialize(fit, NULL))
  actual <- predict(saved, nd, kernel_inputs = f$inputs, summary = FALSE)
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-12)
  expect_identical(predict(saved, nd, kernel_inputs = f$inputs)$point_id,
    as.character(nd$key))
  units <- f$inputs; units$unit_conv <- units$unit_conv / 5
  units$binned$bin_rep <- units$binned$bin_rep * 5
  expect_equal(as.numeric(predict(saved, nd, kernel_inputs = units, summary = FALSE)),
    as.numeric(expected), tolerance = 1e-12)
})

test_that("fitted simulation keeps each observation family's original RNG sequence", {
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    if (family == "binomial") f$data$count <- as.numeric(seq_len(nrow(f$data)) %% 2)
    fit <- fake_bayesian_draws(prepare_bayesian_fixture(f, family), 5)
    set.seed(614); selected <- sample(1:5, 3, replace = TRUE)
    reference <- sapply(selected, function(i) {
      p <- multiScaleR:::.msr_bayes_draw(fit, i)
      eta <- multiScaleR:::.msr_bayes_components(fit$specification, p)$eta
      switch(family,
        gaussian = rnorm(length(eta), eta, exp(p$log_dispersion)),
        poisson = rpois(length(eta), exp(eta)),
        binomial = rbinom(length(eta), fit$specification$trials, plogis(eta)),
        nbinom2 = rnbinom(length(eta), mu = exp(eta), size = exp(p$log_dispersion)))
    })
    expect_identical(as.numeric(as.matrix(simulate(fit, nsim = 3, seed = 614))),
      as.numeric(reference))
  }
})

test_that("seeded prediction restores the absence of a caller RNG state", {
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f, "gaussian"))
  state <- .Random.seed
  on.exit(assign(".Random.seed", state, envir = .GlobalEnv))
  rm(".Random.seed", envir = .GlobalEnv)
  nd <- data.frame(year = "new")
  first <- predict_scales(fit, nd, new_levels = "sample", seed = 319)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_identical(predict_scales(fit, nd, new_levels = "sample", seed = 319), first)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("zeroing a new count-model effect differs from population averaging", {
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f, "poisson"), 6000)
  fit$draws$beta[,] <- 0; fit$draws$beta[, 1] <- log(2)
  fit$draws$log_response_sd[,] <- log(.6)
  nd <- f$data[1, ]; nd$site <- "new_site"; nd$effort <- 1
  sampled <- predict(fit, nd, kernel_inputs = f$inputs, summary = FALSE,
    new_response_levels = "sample", seed = 721)
  zero <- predict(fit, nd, kernel_inputs = f$inputs, summary = FALSE,
    new_response_levels = "zero")
  expect_equal(as.numeric(zero), rep(2, nrow(zero)))
  expect_lt(abs(mean(sampled) - 2 * exp(.6^2 / 2)), .06)
  expect_gt(mean(sampled), mean(zero) + .25)
})

test_that("standard and matrix contrasts retain factor order in both formulas", {
  old <- options("contrasts"); on.exit(options(old))
  for (coding in c("contr.treatment", "contr.sum", "contr.helmert", "contr.poly", "contr.SAS", "matrix")) {
    f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
    f$data$season <- factor(rep(c("early", "middle", "late"), 4),
      levels = c("late", "early", "middle"), ordered = coding == "contr.poly")
    options(contrasts = rep(if (coding == "matrix") "contr.treatment" else coding, 2))
    if (coding == "matrix") contrasts(f$data$season) <-
      cbind(c(1, -1, 0), c(.5, .5, -1))
    fit <- multiScale_bayes(count ~ habitat + season, f$inputs, f$data,
      scales = list(habitat = ~season + (1 | year)), row_id = "key",
      family = "gaussian", method = "prepare")
    fit <- fake_bayesian_draws(fit, 3)
    fit$draws$beta[,] <- rep(seq_len(ncol(fit$draws$beta)) / 2, each = 3)
    fit$draws$scale_beta[,] <- rep(c(-2, .2, .3), each = 3)
    rows <- c(8, 1, 6)
    expected <- predict(fit, summary = FALSE)[, rows]
    sm <- fit$specification$scales$habitat
    log_distance <- drop(sm$X[rows, , drop = FALSE] %*% c(-2, .2, .3)) +
      log(fit$specification$unit_conv)
    options(contrasts = c("contr.treatment", "contr.treatment"))
    nd <- f$data[rows, ]; nd$season <- as.character(nd$season)
    expect_equal(as.numeric(predict(fit, nd, kernel_inputs = f$inputs, summary = FALSE)),
      as.numeric(expected), tolerance = 1e-12, info = coding)
    d <- predict_scales(fit, nd, summary = FALSE, scale = "log")$habitat
    expect_equal(as.numeric(d), rep(unname(log_distance), each = 3), tolerance = 1e-12, info = coding)
  }
})

test_that("seeded calls also restore state after a failed formula evaluation", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  # A random transform deliberately exercises cleanup after RNG consumption.
  # Actual prediction formulas should use deterministic, portable transforms.
  fit <- multiScale_bayes(count ~ habitat + offset(stats::rnorm(length(temp))) + (1 | site),
    f$inputs, f$data, scales = list(habitat = ~1 + (1 | year)),
    row_id = "key", family = "gaussian", method = "prepare")
  fit <- fake_bayesian_draws(fit)
  nd <- f$data[1:3, ]; nd$year <- NULL
  state <- .Random.seed
  expect_error(predict(fit, nd, kernel_inputs = f$inputs, seed = 219), "year")
  expect_identical(.Random.seed, state)
  expect_error(simulate(fit, newdata = nd, kernel_inputs = f$inputs, seed = 219), "year")
  expect_identical(.Random.seed, state)
  scale_fit <- multiScale_bayes(count ~ habitat, f$inputs, f$data,
    scales = list(habitat = ~I(stats::rnorm(length(temp))) + (1 | year)),
    row_id = "key", family = "gaussian", method = "prepare")
  scale_fit <- fake_bayesian_draws(scale_fit); state <- .Random.seed
  expect_error(predict_scales(scale_fit, nd, seed = 219), "year")
  expect_identical(.Random.seed, state)
})

test_that("compact scale summaries agree with full draws for different predictor designs", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  set.seed(291)
  C <- f$inputs$binned$csum$habitat
  f$inputs$binned$vsum$forest <- matrix(rnorm(length(C)), nrow(C)) * C
  f$inputs$binned$csum$forest <- f$inputs$binned$csum$habitat
  fit <- multiScale_bayes(count ~ habitat + forest, f$inputs, f$data,
    scales = list(habitat = ~temp + (1 | year), forest = ~1 + (1 | site)),
    row_id = "key", family = "gaussian", method = "prepare")
  fit <- fake_bayesian_draws(fit, 31)
  fit$draws$scale_beta[,] <- rep(c(-2, .3, -1), each = 31)
  nd <- data.frame(temp = c(1, 1, 2, 2, 1),
    year = c("future_a", "future_a", "future_b", "future_b", "future_a"),
    site = c("new_a", "new_b", "new_a", "new_b", "new_a"))
  x <- predict_scales(fit, nd, new_levels = "sample", summary = FALSE, seed = 931)
  tab <- predict_scales(fit, nd, new_levels = "sample", seed = 931)
  for (name in names(x)) {
    rows <- tab[tab$covariate == name, ]
    rownames(rows) <- NULL
    expected <- multiScaleR:::.msr_bayes_intervals(x[[name]][, rows$row, drop = FALSE])
    expect_equal(rows[c("mean", "median", "sd", "lower", "upper")], expected)
    expect_identical(rows$row, if (name == "habitat") c(1L, 3L) else c(1L, 2L))
  }
  expect_identical(attr(tab, "new_groups"), attr(x, "new_groups"))
})

test_that("collapsing repeated rows preserves explicitly batch-dependent formula terms", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  fit <- multiScale_bayes(count ~ habitat, f$inputs, f$data,
    scales = list(habitat = ~I(temp - mean(temp)) + (1 | year)),
    row_id = "key", family = "gaussian", method = "prepare")
  fit <- fake_bayesian_draws(fit, 11)
  fit$draws$scale_beta[,] <- rep(c(-2, .5), each = 11)
  nd <- data.frame(temp = c(0, 0, 0, 2, 2), year = "future")
  tab <- predict_scales(fit, nd, new_levels = "zero")
  expect_identical(tab$row, c(1L, 4L))
  expect_equal(tab$median, exp(-2 + .5 * (c(0, 2) - mean(nd$temp))) * 1000)
  expect_true(all(tab$distribution == "zero_deviation_reference"))
})

test_that("new-group decisions are explicit and separate for scales and responses", {
  f <- bayesian_scale_fixture();fit <- fake_bayesian_draws(prepare_bayesian_fixture(f,"gaussian"))
  nd <- f$data[1:3,]; nd$year <- rep("future",3)
  nd$site <- c("new","new",as.character(f$data$site[1]))
  call <- function(...) predict(fit,nd,kernel_inputs=f$inputs,row_id="key",...)
  expect_error(call(), "New scale group")
  expect_error(call(new_scale_levels="zero"), "New response group")
  sampled <- call(new_scale_levels="sample",new_response_levels="sample",seed=33)
  expect_true(all(sampled$new_scale_level))
  expect_identical(sampled$new_response_level,c(TRUE,TRUE,FALSE))
  groups <- attr(sampled,"new_groups")
  expect_equal(nrow(groups),2)
  expect_setequal(groups$component,c("scale","response"))
  expect_true(all(groups$action=="sample"))
  nd$site <- NULL
  no_response <- predict(fit,nd,kernel_inputs=f$inputs,row_id="key",re.form=NA,
    new_scale_levels="sample",seed=33)
  expect_false(any(no_response$new_response_level))
  nd$year <- NULL
  reference <- predict(fit,nd,kernel_inputs=f$inputs,row_id="key",re.form=NA,scale_re.form=NA)
  expect_false(any(reference$new_scale_level | reference$new_response_level))
  expect_error(predict(fit,nd,kernel_inputs=f$inputs,row_id="key",re.form=~1), "NULL or NA")
  expect_error(predict(fit,nd,kernel_inputs=f$inputs,row_id="key",scale_re.form=TRUE), "NULL or NA")
})

test_that("new annual scales share effects and retain pooling-SD uncertainty", {
  f <- bayesian_scale_fixture();n <- 6000L
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f,"gaussian"),n)
  set.seed(401); alpha <- rnorm(n,log(.25),.3)
  sigma <- rep(c(.1,.8),length.out=n)
  fit$draws$scale_beta[,1] <- alpha
  fit$draws$log_scale_sd[,1] <- log(sigma)
  nd <- data.frame(year=c("new_a","new_a","new_b"))
  expect_error(predict_scales(fit,nd), "New scale group")
  x <- predict_scales(fit,nd,new_levels="sample",summary=FALSE,scale="log",seed=491)$habitat
  expect_identical(x[,1],x[,2])
  expect_lt(abs(var(x[,1])-(var(alpha)+mean(sigma^2))),.025)
  expect_lt(abs(cov(x[,1],x[,3])-var(alpha)),.025)
  permutation <- c(3,1,2)
  permuted <- predict_scales(fit,nd[permutation,,drop=FALSE],new_levels="sample",summary=FALSE,scale="log",seed=491)$habitat
  expect_identical(as.numeric(permuted[,order(permutation)]),as.numeric(x))
  zero <- predict_scales(fit,nd,new_levels="zero",summary=FALSE,scale="log")$habitat
  expect_equal(zero[,1],alpha+log(fit$specification$unit_conv))
  expect_identical(zero[,1],zero[,3])
  tab <- predict_scales(fit,nd,new_levels="sample",draws=100,scale="log",seed=31)
  expect_equal(nrow(tab),2)
  expect_true(all(tab$distribution=="new_group_predictive"))
  expect_true(all(tab$new_level))
  expect_true(all(predict_scales(fit,nd,new_levels="zero",draws=100)$distribution=="zero_deviation_reference"))
  expect_true(all(predict_scales(fit,data.frame(dummy=1),re.form=NA,draws=100)$distribution=="zero_deviation_reference"))
  # Correlated location and SD draws must not be shuffled separately.
  fit$draws$scale_beta[,1] <- rep(log(c(.1,.4)),length.out=n)
  fit$draws$log_scale_sd[,1] <- rep(log(c(.1,.9)),length.out=n)
  distance <- predict_scales(fit,nd[1,,drop=FALSE],new_levels="sample",summary=FALSE,seed=92)$habitat
  expected <- mean(exp(fit$draws$scale_beta[,1]+exp(fit$draws$log_scale_sd[,1])^2/2))*fit$specification$unit_conv
  expect_lt(abs(mean(distance)-expected),20)
})

test_that("new response effects are shared within a replicate and observation noise is separate", {
  f <- bayesian_scale_fixture();fit <- fake_bayesian_draws(prepare_bayesian_fixture(f,"gaussian"),1)
  fit$draws$beta[,] <- 0;fit$draws$log_response_sd[,] <- log(.6)
  fit$draws$log_dispersion[,] <- log(.2)
  nd <- f$data[1:3,];nd$site <- c("new_a","new_a","new_b");nd$effort <- 1
  response <- simulate(fit,nsim=6000,seed=81,newdata=nd,kernel_inputs=f$inputs,row_id="key",new_response_levels="sample")
  expect_lt(abs(var(as.numeric(response[1,]))-.4),.025)
  expect_lt(abs(cov(as.numeric(response[1,]),as.numeric(response[2,]))-.36),.025)
  expect_lt(abs(cov(as.numeric(response[1,]),as.numeric(response[3,]))),.025)
  ref <- simulate(fit,nsim=3000,seed=81,newdata=nd,kernel_inputs=f$inputs,row_id="key",new_response_levels="zero")
  expect_lt(abs(var(as.numeric(ref[1,]))-.04),.006)
  expect_identical(attr(response,"point_ids"),as.character(nd$key))
})

test_that("response families and new binomial trial totals have explicit units", {
  for (family in c("gaussian","poisson","nbinom2","binomial")) {
    f <- bayesian_scale_fixture()
    if(family=="binomial") f$data$count <- rep(0:1,length.out=nrow(f$data))
    fit <- fake_bayesian_draws(prepare_bayesian_fixture(f,family),1)
    fit$draws$beta[,] <- 0
    target <- if(family=="gaussian")1.2 else if(family=="binomial").25 else 3
    fit$draws$beta[,1] <- if(family=="gaussian")target else if(family=="binomial")qlogis(target) else log(target)
    if(family=="nbinom2")fit$draws$log_dispersion[,] <- log(2)
    nd <- f$data[1:3,];nd$effort <- 1;nd$site <- "new"
    mu <- predict(fit,nd,kernel_inputs=f$inputs,row_id="key",new_response_levels="zero")
    expect_equal(mu$median,rep(target,3))
    y <- simulate(fit,nsim=2000,seed=71,newdata=nd,kernel_inputs=f$inputs,row_id="key",new_response_levels="zero")
    expect_lt(abs(mean(as.numeric(y[1,]))-target),.2)
    if(family!="gaussian")expect_true(all(as.matrix(y)>=0 & as.matrix(y)==floor(as.matrix(y))))
    if(family=="binomial")expect_true(all(as.matrix(y)<=1))
  }
  f <- bayesian_scale_fixture();f$data$count <- rep(1:3,length.out=nrow(f$data));f$data$failures <- 6-f$data$count
  fit <- multiScale_bayes(cbind(count,failures)~habitat+temp+(1|site),f$inputs,f$data,
    scales=list(habitat=~1+(1|year)),row_id="key",family="binomial",method="prepare")
  fit <- fake_bayesian_draws(fit,1);fit$draws$beta[,] <- 0
  nd <- f$data[1:3,];nd$count <- nd$failures <- NULL;nd$site <- "new"
  expect_error(simulate(fit,newdata=nd,kernel_inputs=f$inputs,row_id="key",new_response_levels="zero"),"trials explicitly")
  y <- simulate(fit,nsim=2000,seed=82,newdata=nd,kernel_inputs=f$inputs,row_id="key",new_response_levels="zero",trials=c(2,5,9))
  expect_true(all(as.matrix(y)>=0 & as.matrix(y)<=c(2,5,9)))
  expect_lt(max(abs(rowMeans(y)-c(1,2.5,4.5))),.15)
  expect_error(simulate(fit,newdata=nd,kernel_inputs=f$inputs,row_id="key",new_response_levels="zero",trials=c(2,4)),"positive integer totals")
  expect_error(simulate(fit,newdata=nd,kernel_inputs=f$inputs,row_id="key",new_response_levels="zero",trials=2.5),"positive integer totals")
})

test_that("seeds, errors and observed-response changes preserve prediction state", {
  f <- bayesian_scale_fixture();fit <- fake_bayesian_draws(prepare_bayesian_fixture(f,"gaussian"),10);before <- fit
  nd <- f$data[1:3,];nd$year <- "future";nd$site <- "new"
  set.seed(805);state <- .Random.seed
  p <- predict(fit,nd,kernel_inputs=f$inputs,row_id="key",new_scale_levels="sample",new_response_levels="sample",seed=41)
  expect_identical(.Random.seed,state)
  expect_identical(p,predict(fit,nd,kernel_inputs=f$inputs,row_id="key",new_scale_levels="sample",new_response_levels="sample",seed=41))
  changed <- nd;changed$count <- NA
  expect_identical(p,predict(fit,changed,kernel_inputs=f$inputs,row_id="key",new_scale_levels="sample",new_response_levels="sample",seed=41))
  predict_scales(fit,nd,new_levels="sample",seed=4)
  expect_identical(.Random.seed,state)
  bad <- fit;bad$draws$log_response_sd[,] <- 1000
  expect_error(predict(bad,nd,kernel_inputs=f$inputs,row_id="key",new_scale_levels="sample",new_response_levels="sample",seed=41),"Nonfinite")
  expect_identical(.Random.seed,state)
  predict(fit,nd,kernel_inputs=f$inputs,row_id="key",new_scale_levels="sample",new_response_levels="sample")
  expect_false(identical(.Random.seed,state))
  expect_identical(fit,before)
  expect_error(predict_scales(fit,nd,new_levels="invent"),"arg")
  expect_error(predict_scales(fit,nd,seed=0,new_levels="sample"),"seed")
})

test_that("fixed-factor levels and trained contrasts cannot be invented or relearned", {
  f <- bayesian_scale_fixture(n_site=6,n_year=4)
  f$data$region <- factor(rep(c("a","b","c"),length.out=nrow(f$data)))
  old <- options("contrasts")
  tryCatch({
    options(contrasts=c("contr.sum","contr.poly"))
    fit <- multiScale_bayes(count~habitat+region,f$inputs,f$data,
      scales=list(habitat=~1+(1|year)),row_id="key",family="gaussian",method="prepare")
    fit <- fake_bayesian_draws(fit,3);fit$draws$beta[,] <- 0
    fit$draws$beta[,2:3] <- rep(c(2,-3),each=3)
    options(contrasts=c("contr.treatment","contr.poly"))
    nd <- f$data[c(1,2,3),];nd$count <- NULL
    expected <- c(a=2,b=-3,c=1)[as.character(nd$region)]
    p <- predict(fit,nd,kernel_inputs=f$inputs,row_id="key",summary=FALSE)
    expect_equal(as.numeric(p[1,]),unname(expected))
    legacy <- fit;legacy$specification$response$prediction <- NULL
    expect_error(predict(legacy,nd,kernel_inputs=f$inputs,row_id="key"),"contrast coding")
    nd$region <- "unseen"
    expect_error(predict(fit,nd,kernel_inputs=f$inputs,row_id="key",new_response_levels="sample"),"new level")
  },finally=options(old))
  annual <- multiScale_bayes(count~habitat+year,f$inputs,f$data,
    scales=list(habitat=~1+(1|year)),row_id="key",family="gaussian",method="prepare")
  annual <- fake_bayesian_draws(annual)
  nd <- f$data[1,,drop=FALSE];nd$year <- "future"
  expect_error(predict(annual,nd,kernel_inputs=f$inputs,row_id="key",new_scale_levels="sample",new_response_levels="sample"),"new level")
  independent <- multiScale_bayes(count~habitat,f$inputs,f$data,
    scales=list(habitat=~0+year),row_id="key",family="gaussian",method="prepare")
  independent <- fake_bayesian_draws(independent)
  expect_error(predict_scales(independent,nd,new_levels="sample"),"new level")
})

test_that("numeric bases and scaling retain their fitting definitions for a single new row", {
  f <- bayesian_scale_fixture();nd <- f$data[1,,drop=FALSE];nd$temp <- 3
  for(transform in c("poly","scale")) {
    formula <- if(transform=="poly")count~habitat+poly(temp,2) else count~habitat+scale(temp)
    fit <- multiScale_bayes(formula,f$inputs,f$data,scales=list(habitat=~1),
      row_id="key",family="gaussian",method="prepare")
    fit <- fake_bayesian_draws(fit,3);fit$draws$beta[,] <- 0
    if(transform=="poly") {
      fit$draws$beta[,2:3] <- rep(c(2,-1),each=3)
      basis <- predict(poly(f$data$temp,2),newdata=nd$temp)
      expected <- drop(basis%*%c(2,-1))
    } else {
      fit$draws$beta[,2] <- 2
      expected <- 2*(nd$temp-mean(f$data$temp))/sd(f$data$temp)
    }
    predicted <- predict(fit,nd,kernel_inputs=f$inputs,row_id="key",summary=FALSE)
    expect_equal(as.numeric(predicted),rep(expected,3),tolerance=1e-12)
  }
})

test_that("different predictor and grouping terms receive distinct new effects", {
  f <- bayesian_scale_fixture()
  C <- f$inputs$binned$csum$habitat
  set.seed(341);forest <- matrix(rnorm(length(C)),nrow(C))
  f$inputs$binned$vsum$forest <- forest*C;f$inputs$binned$csum$forest <- C
  fit <- multiScale_bayes(count~habitat+forest,f$inputs,f$data,
    scales=list(habitat=~1+(1|year),forest=~1+(1|year)+(1|site)),
    row_id="key",family="gaussian",method="prepare")
  fit <- fake_bayesian_draws(fit,5000)
  fit$draws$log_scale_sd[,] <- rep(log(c(.5,.5,.4)),each=5000)
  nd <- data.frame(year=c("future","future","future"),site=c("new_a","new_a","new_b"))
  x <- predict_scales(fit,nd,new_levels="sample",summary=FALSE,scale="log",seed=291)
  expect_identical(x$habitat[,1],x$habitat[,3])
  expect_identical(x$forest[,1],x$forest[,2])
  expect_lt(abs(cov(x$habitat[,1],x$forest[,1])),.025)
  expect_lt(abs(var(x$forest[,1])-.41),.025)
  expect_lt(abs(cov(x$forest[,1],x$forest[,3])-.25),.025)
  meta <- attr(x,"new_groups")
  expect_equal(sum(meta$covariate=="habitat"),1)
  expect_equal(sum(meta$covariate=="forest"),3)
})
