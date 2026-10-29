bayesian_scale_fixture <- function(n_site = 8L, n_year = 4L, annual = TRUE,
                                   kernel = "gaussian", seed = 916L) {
  set.seed(seed)
  dat <- expand.grid(site = paste0("s", seq_len(n_site)),
                      year = paste0("y", seq_len(n_year)),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  dat$site <- factor(dat$site)
  dat$year <- factor(dat$year)
  dat$key <- if (annual) paste(dat$site, dat$year, sep = "_") else as.character(dat$site)
  ids <- unique(dat$key)
  d <- seq(0.02, 1, length.out = 15L)
  n <- length(ids)
  C <- matrix(sample(1:5, n * length(d), replace = TRUE), n)
  C[seq_len(min(3L, n)), 1:2] <- 0
  habitat <- outer(rnorm(n), sin(8 * d)) + outer(rnorm(n), cos(13 * d)) +
    matrix(rnorm(n * length(d), sd = 0.2), n)
  b <- list(bin_rep = d, vsum = list(habitat = habitat * C),
            csum = list(habitat = C), point_ids = ids, n_points = n)
  inputs <- list(binned = b, min_D = 40, max_D = 1000, unit_conv = 1000,
                  kernel = kernel)
  dat$temp <- rnorm(nrow(dat))
  dat$effort <- runif(nrow(dat), 0.7, 1.3)
  dat$count <- rpois(nrow(dat), 3)
  list(data = dat, inputs = inputs)
}

prepare_bayesian_fixture <- function(f, family = "nbinom2", scales = NULL, ...) {
  if (is.null(scales)) scales <- list(habitat = ~1 + (1 | year))
  multiScale_bayes(count ~ habitat + temp + (1 | site) + offset(log(effort)),
    f$inputs, f$data, scales = scales, row_id = "key", family = family,
    method = "prepare", ...)
}

fake_bayesian_draws <- function(fit, n = 8L) {
  fit$draws <- lapply(fit$parameters, function(z)
    matrix(rep(z, each = n), nrow = n))
  fit$inference <- "test fixture"
  fit
}
