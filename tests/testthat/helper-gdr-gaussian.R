gdr_gaussian_fixture <- function(block, observer, K, unused = TRUE) {
  skip_if_not_installed("unmarked")
  skip_if(as.character(packageVersion("unmarked")) != "1.5.2", "Native layout needs qualification")
  n <- c(0, 1, 2, 0, 1, 2, 0, 1)
  yd <- cbind(pmin(n, 1), pmax(n - 1, 0))
  yr <- cbind(pmin(n, 1), pmax(n - 1, 0), 0)
  sc <- data.frame(h = c(-.1, .1, .3, .5, .65, .8, 1, 1.2),
    cover = factor(rep(c("A", "B"), each = 4)),
    Observer = factor(rep(c("X", "Y"), 4)),
    block = factor(rep(c("a", "b"), 4), levels = c("a", "b", "unused")))
  contrasts(sc$cover) <- contr.sum(2)
  frame <- unmarked::unmarkedFrameGDR(yd, yr, siteCovs = sc,
    dist.breaks = c(0, 25, 75), unitsIn = "m", period.lengths = c(1, 2, 1))
  ids <- paste0("toy", seq_len(8))
  frame <- expand_gdr_intervals(frame, data.frame(survey = ids), "survey", "minutes",
    site_covariates = c(rem_hab = "h"))$data
  forms <- list(lambda = ~cover * splines::ns(h, knots = c(.2, .6), Boundary.knots = c(0, 1)),
    phi = ~1, dist = if (observer) ~h + (1 | Observer) else ~h,
    rem = ~rem_hab + I(interval == 1))
  if (block) forms$lambda <- update(forms$lambda, ~. + (1 | block))
  original <- get("fit_TMB", asNamespace("unmarked"))
  constructor <- unmarked::gdistremoval
  local <- new.env(parent = environment(constructor))
  local$fit_TMB <- function(model, data, params, random, starts, method, ...)
    stop(structure(list(message = "native numeric preparation only", call = NULL,
      model = model, data = data, parameters = params, random = random),
      class = c("msr_test_gaussian_GDR", "error", "condition")))
  environment(constructor) <- local
  captured <- tryCatch(constructor(lambdaformula = forms$lambda, phiformula = forms$phi,
    distanceformula = forms$dist, removalformula = forms$rem, data = frame,
    K = K, keyfun = "halfnorm", mixture = "P", output = "density", unitsOut = "kmsq",
    engine = "TMB", threads = 1), msr_test_gaussian_GDR = identity)
  stopifnot(inherits(captured, "msr_test_gaussian_GDR"),
    identical(get("fit_TMB", asNamespace("unmarked")), original))
  if (block && unused) {
    captured$data$Z_lambda <- cbind(captured$data$Z_lambda,
      Matrix::Matrix(0, 8L, 1L, sparse = TRUE))
    captured$data$n_grouplevels_lambda <- captured$data$n_grouplevels_lambda + 1L
    captured$parameters$b_lambda <- c(captured$parameters$b_lambda, 0)
  }
  gd <- get("getDesign", asNamespace("unmarked"))(frame, forms)
  map <- function(Z) {
    if (!ncol(Z)) return(integer())
    z <- Matrix::summary(Z)
    stopifnot(nrow(z) == 8L, all(z$x == 1), !anyDuplicated(z$i))
    as.integer(z$j[match(seq_len(8L), z$i)])
  }
  pars <- captured$parameters
  pars$beta_lambda <- seq_len(ncol(gd$X_lambda)) * .01
  pars$beta_lambda[1L] <- log(8/(pi * .075^2))
  pars$beta_dist <- c(log(40), .1)
  pars$beta_rem <- c(qlogis(.25), .1, .2)
  pars$b_lambda <- if (block) c(.1, -.2, if (unused) .3) else numeric()
  pars$b_dist <- if (observer) c(.08, -.05) else numeric()
  pars$lsigma_lambda <- if (block) log(.3) else numeric()
  pars$lsigma_dist <- if (observer) log(.4) else numeric()
  data <- list(Xlambda = gd$X_lambda, Xdist = gd$X_dist, Xrem = gd$X_rem,
    yd = unname(frame@yDistance), yr = unname(frame@yRemoval),
    block = map(gd$Z_lambda), observer = map(gd$Z_dist),
    periods = c(1, 2, 1), breaks = c(0, 25, 75), area = captured$data$A, K = K,
    annulus_area_m2 = captured$data$a, annulus_weights = captured$data$u,
    row_ids = list(survey = ids, removal_survey = rep(ids, each = 3L),
      removal_interval = rep(1:3, 8L)))
  list(data = data, pars = pars, captured = captured)
}

gdr_gaussian_require_RTMB <- function() {
  skip_if_not_installed("RTMB")
  skip_if(as.character(packageVersion("RTMB")) != "2.0", "Private derivative layout needs qualification")
}

gdr_gaussian_numeric <- function(data, p, uniform = FALSE) {
  eta <- as.vector(data$Xlambda %*% p$beta_lambda)
  if (length(data$block)) eta <- eta + p$b_lambda[data$block]
  eta_dist <- as.vector(data$Xdist %*% p$beta_dist)
  if (length(data$observer)) eta_dist <- eta_dist + p$b_dist[data$observer]
  pd <- .msr_gdr_exact_distance_bins(exp(eta_dist), data$breaks,
    if (uniform) "uniform" else "halfnorm")
  removal <- matrix(plogis(as.vector(data$Xrem %*% p$beta_rem)), 8L, 3L, byrow = TRUE)
  pr <- .msr_gdr_exact_removal_bins(removal, data$periods)
  ll <- .msr_gdr_finite_poisson_score_bins(exp(eta) * data$area,
    pd, pr, data$yd, data$yr, data$K)$logLik
  if (length(data$block)) ll <- ll + sum(dnorm(p$b_lambda, 0, exp(p$lsigma_lambda), log = TRUE))
  if (length(data$observer)) ll <- ll + sum(dnorm(p$b_dist, 0, exp(p$lsigma_dist), log = TRUE))
  -ll
}

