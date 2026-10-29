profile_glm_fixture <- function(n = 240L, seed = 723L, theta = 2, effect = 1.1,
                                kernel = "gaussian") {
  set.seed(seed)
  d <- seq(0, 1, length.out = 30)
  C <- matrix(sample(1:6, n * length(d), TRUE), n)
  radial <- outer(rnorm(n), sin(7 * d)) + outer(rnorm(n), cos(13 * d)) +
    outer(rnorm(n), d) + matrix(rnorm(n * length(d), sd = .3), n)
  b <- list(bin_rep = d, vsum = list(hab = radial * C), csum = list(hab = C),
            point_ids = as.character(seq_len(n)), n_points = n)
  inputs <- list(raw_cov = NULL, d_list = NULL, binned = b, min_D = 50,
                 max_D = 800, unit_conv = 1000, kernel = kernel)
  engine <- multiScaleR:::.msr_prepare_kernel_engine(b, kernel = kernel)
  kp <- c(log(.25), if (kernel == "expow") log(1.8))
  hab <- multiScaleR:::.msr_eval_kernel_engine(engine, kp)$value[, 1]
  dat <- data.frame(hab = hab, z = rnorm(n), group = factor(rep(c("a", "b"), length.out = n)),
                    exposure = runif(n, .8, 1.2))
  mu <- exp(.5 + effect * hab + .2 * dat$z + .3 * (dat$group == "b") + log(dat$exposure))
  dat$y <- if (is.infinite(theta)) stats::rpois(n, mu) else
    stats::rnbinom(n, mu = mu, size = theta)
  model <- stats::glm(y ~ hab + z + group + offset(log(exposure)), data = dat,
                     family = MASS::negative.binomial(2))
  list(inputs = inputs, data = dat, model = model, mu = mu,
       truth = c(log_scale = log(.25), beta = c(.5, .2, .3, effect), log_theta = log(theta)))
}

profile_glm_fit <- function(fixture, ...) {
  multiScale_optim(fixture$model, fixture$inputs, backend = "profile_glm",
                  n_cores = 1, verbose = FALSE, ...)
}
