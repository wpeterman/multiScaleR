# Independent calculations using only returned parameter draws and design data.
# These references do not call a package prediction, kernel, or likelihood routine.
log_mean_exp <- function(x) {
  m <- max(x)
  if (is.infinite(m)) m else m + log(mean(exp(x - m)))
}

log_add <- function(a, b) {
  m <- pmax(a, b)
  ok <- is.finite(m)
  m[ok] <- m[ok] + log1p(exp(-abs(a[ok] - b[ok])))
  m
}

bin_transform <- function(input, ids) {
  cuts <- seq(0, 1 + 1e-9, length.out = 13L)
  point <- unique(input$data[ids, "key"])
  train <- unlist(input$inputs$d_list[point], use.names = FALSE)
  code <- cut(train, cuts, include.lowest = TRUE, labels = FALSE)
  reps <- vapply(seq_len(12L), function(j) {
    x <- train[code == j]
    if (length(x)) mean(x) else mean(cuts[j + 0:1])
  }, numeric(1L))
  w <- exp(-.5 * (reps / .26)^2)
  x <- drop(input$values[input$data[ids, "key"], , drop = FALSE] %*% (w / sum(w)))
  c(list(distances = reps), list(center = mean(x), spread = sd(x)))
}

fit_predictor <- function(fit, input, fold, model) {
  draws <- fit$draws
  split <- input$folds[[fold]]
  train <- input$data[split$analysis, ]
  test <- input$data[split$assessment, ]
  ref <- bin_transform(input, split$analysis)
  spec <- fit$specification
  expected_beta <- c("(Intercept)", "habitat")
  stopifnot(identical(spec$beta_names, expected_beta),
    identical(spec$observation_ids, as.character(train$key)),
    identical(spec$y, as.numeric(train$y)),
    isTRUE(all.equal(unname(spec$center), ref$center, tolerance = 1e-12)),
    isTRUE(all.equal(unname(spec$spread), ref$spread, tolerance = 1e-12)),
    identical(spec$response$random[[1L]]$levels, sort(unique(as.character(train$site)))))
  n <- nrow(draws$beta)
  site <- match(test$site, spec$response$random[[1L]]$levels)
  site_sd <- exp(draws$log_response_sd[, 1L])
  test_sites <- unique(test$site[is.na(site)])
  stopifnot(length(test_sites) == as.integer(fold == "seen_year_new_site"))
  values <- input$values[as.character(test$key), , drop = FALSE]
  lscale <- matrix(draws$scale_beta[, 1L], n, nrow(test))
  bandwidth_new <- logical(nrow(test))
  scale_sd <- rep(0, n)
  if (model == "annual") {
    term <- spec$scales$habitat$random[[1L]]
    stopifnot(identical(term$levels, sort(unique(as.character(train$year)))))
    year <- match(test$year, term$levels)
    scale_sd <- exp(draws$log_scale_sd[, 1L])
    bandwidth_new <- is.na(year)
    for (j in which(!bandwidth_new)) lscale[, j] <- lscale[, j] + scale_sd *
      draws$scale_z[, year[j]]
  }
  eta_fixed <- outer(draws$beta[, 1L], test$offset, "+")
  for (j in which(!is.na(site))) eta_fixed[, j] <- eta_fixed[, j] +
    site_sd * draws$response_z[, site[j]]
  centered_square <- ref$distances^2 - min(ref$distances)^2
  kernel <- function(x, val) {
    w <- exp(-.5 * outer(exp(-2 * x), centered_square))
    drop(w %*% val) / rowSums(w)
  }
  cached <- eta_fixed
  for (j in which(!bandwidth_new)) cached[, j] <- eta_fixed[, j] +
    draws$beta[, 2L] * (kernel(lscale[, j], values[j, ]) - ref$center) / ref$spread
  # There is at most one unobserved latent effect in the selected assessment
  # blocks. The same node is applied to every visit sharing that effect.
  predictor <- function(z) {
    eta <- cached
    if (any(bandwidth_new)) {
      multiplier <- exp(-2 * (lscale[, which(bandwidth_new)[1L]] + z * scale_sd))
      w <- exp(-.5 * outer(multiplier, centered_square))
      x <- (w %*% t(values[bandwidth_new, , drop = FALSE])) / rowSums(w)
      eta[, bandwidth_new] <- eta_fixed[, bandwidth_new] +
        draws$beta[, 2L] * (x - ref$center) / ref$spread
    }
    if (anyNA(site)) eta[, is.na(site)] <- eta[, is.na(site)] + z * site_sd
    eta
  }
  list(predict = predictor, dispersion = exp(draws$log_dispersion[, 1L]),
    n = n, uses_latent = any(bandwidth_new) || anyNA(site),
    analytic_site_sd = if (input$family == "gaussian" && anyNA(site)) site_sd else NULL,
    new_scale = bandwidth_new, new_response = is.na(site), ref = ref,
    y = test$y, y_stress = if (input$family == "gaussian") test$y + 2.4 else test$y * 2 + 8L,
    family = input$family)
}

observation_ll <- function(eta, y, dispersion, family) {
  target <- matrix(y, nrow(eta), ncol(eta), byrow = TRUE)
  ll <- if (family == "gaussian") dnorm(target, eta, dispersion, log = TRUE) else
    dnbinom(target, mu = exp(eta), size = dispersion, log = TRUE)
  rowSums(ll)
}

# Conditional on the returned draws, these quadratures integrate exactly one
# new-group standard normal. The refinement difference audits the result.
reference_score <- function(f, uniforms, order) {
  if (length(f$analytic_site_sd)) {
    eta <- f$predict(0)
    n <- ncol(eta)
    v <- f$dispersion^2
    b <- f$analytic_site_sd^2
    # Integrate the rank-one covariance v * I + b * J without cancellation.
    ll <- function(y) {
      residual <- matrix(y, nrow(eta), n, byrow = TRUE) - eta
      average <- rowMeans(residual)
      -0.5 * (n * log(2 * pi) + (n - 1L) * log(v) + log(v + n * b) +
        rowSums((residual - average)^2) / v + n * average^2 / (v + n * b))
    }
    joint <- ll(f$y)
    tail <- ll(f$y_stress)
    y <- matrix(f$y, nrow(eta), n, byrow = TRUE)
    cdf <- pnorm(y, eta, sqrt(v + b))
  } else {
    g <- if (f$uses_latent) normal_quadrature(order) else
    list(nodes = 0, weights = 1)
    joint <- tail <- rep(-Inf, f$n)
    cdf <- matrix(0, f$n, length(f$y))
    y <- matrix(f$y, f$n, length(f$y), byrow = TRUE)
    for (k in seq_along(g$nodes)) {
      eta <- f$predict(g$nodes[k])
      w <- g$weights[k]
      joint <- log_add(joint, log(w) + observation_ll(eta, f$y, f$dispersion, f$family))
      tail <- log_add(tail, log(w) + observation_ll(eta, f$y_stress, f$dispersion, f$family))
      prob <- if (f$family == "gaussian") pnorm(y, eta, f$dispersion) else
        pnbinom(y - 1L, size = f$dispersion, mu = exp(eta)) +
        sweep(dnbinom(y, size = f$dispersion, mu = exp(eta)), 2L, uniforms, "*")
      cdf <- cdf + w * prob
    }
  }
  list(observed = -log_mean_exp(joint), stress = -log_mean_exp(tail), pit = colMeans(cdf),
    effective_density_draws = c(observed = exp(2 * log_mean_exp(joint) -
      log_mean_exp(2 * joint)) * f$n, stress = exp(2 * log_mean_exp(tail) -
      log_mean_exp(2 * tail)) * f$n))
}

# Composite Gauss-Legendre concentrates nodes where a narrow likelihood peaks.
# Using the normal density as a weight permits a dense, refined central grid.
normal_quadrature <- function(order) {
  g <- statmod::gauss.quad(order, "legendre")
  edges <- c(-12, -8, -6, seq(-5, 5, by = .5), 6, 8, 12)
  half_width <- diff(edges) / 2
  center <- head(edges, -1L) + half_width
  nodes <- c(outer(g$nodes, half_width) + rep(center, each = order))
  weights <- c(outer(g$weights, half_width)) * dnorm(nodes)
  stopifnot(abs(sum(weights) - 1) < 1e-10)
  list(nodes = nodes, weights = weights / sum(weights))
}

# Match the package seed/draw order for budget 1, then use a nested stream for
# the precision check. Repeat scoring at independent seeds without refitting.
augmented_scores <- function(f, seed, budgets = c(1L, 8L, 32L)) {
  set.seed(seed)
  z <- if (f$uses_latent) matrix(rnorm(f$n * max(budgets)), nrow = f$n) else
    matrix(0, f$n, 1L)
  if (!f$uses_latent) budgets <- 1L
  observed <- stress <- rep(-Inf, f$n)
  out <- matrix(NA_real_, length(budgets), 2L, dimnames = list(NULL, c("observed", "stress")))
  for (j in seq_len(max(budgets))) {
    eta <- f$predict(z[, j])
    observed <- log_add(observed, observation_ll(eta, f$y, f$dispersion, f$family))
    stress <- log_add(stress, observation_ll(eta, f$y_stress, f$dispersion, f$family))
    if (j %in% budgets) out[match(j, budgets), ] <-
      c(-log_mean_exp(observed), -log_mean_exp(stress)) + log(j)
  }
  data.frame(budget = budgets, seed = seed, out, row.names = NULL)
}

# The generative oracle knows the true global constants and realized training
# year bandwidths. It conditions on observed training visits for known sites,
# and integrates the held-out year and/or site. No held-out outcomes update it.
oracle_score <- function(input, fold, stress, uniforms) {
  test <- input$data[input$folds[[fold]]$assessment, ]
  train <- input$data[input$folds[[fold]]$analysis, ]
  p <- input$truth
  gaussian <- input$family == "gaussian"
  eta <- function(d, scale) {
    w <- exp(-.5 * outer(rep_len(exp(-2 * scale), nrow(d)), input$distances^2))
    x <- rowSums(w * input$values[d$key, , drop = FALSE]) / rowSums(w)
    p$beta0 + p$beta1 * (x - p$center) / p$spread + d$offset
  }
  train_eta <- eta(train, input$realized_year[as.character(train$year)])
  by_site <- split(seq_len(nrow(test)), test$site, drop = TRUE)
  conditional <- lapply(names(by_site), function(s) {
    ii <- which(train$site == s)
    if (gaussian) {
      sd <- 1 / sqrt(1 / p$site_sd^2 + length(ii) / p$dispersion^2)
      return(list(mean = sd^2 * sum(train$y[ii] - train_eta[ii]) / p$dispersion^2, sd = sd))
    }
    # The NB2 log conditional is strictly concave. Adapt the quadrature to its
    # mode and curvature, with a slightly widened Gaussian proposal.
    log_site <- function(b) vapply(b, function(x) sum(dnbinom(train$y[ii],
      mu = exp(train_eta[ii] + x), size = p$dispersion, log = TRUE)) +
      dnorm(x, 0, p$site_sd, log = TRUE), numeric(1L))
    if (!length(ii)) return(list(mean = 0, sd = p$site_sd, log_site = log_site))
    opt <- optimize(log_site, interval = c(-8, 8), maximum = TRUE, tol = 1e-11)
    mu <- exp(train_eta[ii] + opt$maximum)
    curvature <- sum((train$y[ii] + p$dispersion) * mu * p$dispersion /
      (mu + p$dispersion)^2) + 1 / p$site_sd^2
    list(mean = opt$maximum, sd = 1.25 / sqrt(curvature), log_site = log_site)
  })
  evaluate <- function(order) {
    node <- statmod::gauss.quad.prob(order, "normal")
    year <- if (fold == "whole_year" && p$annual_sd > 0) node else
      list(nodes = 0, weights = 1)
    log_w <- log(year$weights)
    outer_ll <- outer_stress <- rep(0, length(log_w))
    cdf <- numeric(nrow(test))
    for (j in seq_along(log_w)) {
      scale <- if (fold == "whole_year") p$log_scale + p$annual_sd * year$nodes[j] else
        input$realized_year[as.character(test$year)]
      test_eta <- eta(test, scale)
      for (k in seq_along(by_site)) {
        b <- conditional[[k]]
        r <- by_site[[k]]
        z <- b$mean + b$sd * node$nodes
        lp <- log(node$weights)
        if (!gaussian) lp <- lp + b$log_site(z) - dnorm(z, b$mean, b$sd, log = TRUE)
        lp <- lp - log_mean_exp(lp) - log(length(lp))
        e <- outer(z, test_eta[r], "+")
        ll <- function(y) observation_ll(e, y[r], p$dispersion, input$family)
        outer_ll[j] <- outer_ll[j] + log_mean_exp(lp + ll(test$y)) + log(length(lp))
        outer_stress[j] <- outer_stress[j] + log_mean_exp(lp + ll(stress)) + log(length(lp))
        y <- matrix(test$y[r], nrow(e), ncol(e), byrow = TRUE)
        prob <- if (gaussian) pnorm(y, e, p$dispersion) else
          pnbinom(y - 1L, size = p$dispersion, mu = exp(e)) +
          sweep(dnbinom(y, size = p$dispersion, mu = exp(e)), 2L, uniforms[r], "*")
        cdf[r] <- cdf[r] + year$weights[j] * colSums(exp(lp) * prob)
      }
    }
    list(observed = -log_mean_exp(log_w + outer_ll) - log(length(log_w)),
      stress = -log_mean_exp(log_w + outer_stress) - log(length(log_w)), pit = cdf)
  }
  coarse <- evaluate(65L)
  fine <- evaluate(129L)
  list(score = fine, refine = c(observed = abs(coarse$observed - fine$observed),
    stress = abs(coarse$stress - fine$stress), pit = max(abs(coarse$pit - fine$pit))))
}
