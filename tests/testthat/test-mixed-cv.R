mixed_cv_fixture <- function() {
  set.seed(10)
  group <- factor(rep(letters[1:18], each = 8))
  x <- rnorm(length(group))
  exposure <- runif(length(x), .8, 1.25)
  b <- rnorm(18, sd = .7)[group]
  slope <- rnorm(18, sd = .55)[group]
  eta <- .9 + .4 * x + b
  trials <- rep(c(4L, 5L), length.out = length(x))
  successes <- rbinom(length(x), trials, plogis(eta - 1))
  data <- data.frame(group, x, exposure, y = 1 + x + b + slope * x +
    log(exposure) + rnorm(length(x), sd = .5),
    poisson = rpois(length(x), exp(eta)),
    successes, failures = trials - successes,
    nb = rnbinom(length(x), size = 1.1, mu = exp(eta)),
    gamma = rgamma(length(x), shape = 2, scale = exp(eta) / 2),
    row.names = paste0("r", seq_along(x)))
  ids <- rownames(data)
  list(data = data, new = grouped_cv_folds(data$group, k = 3, ids = ids)[1],
    existing = list(one = list(analysis = ids[seq_along(ids) %% 4 != 0L],
      assessment = ids[seq_along(ids) %% 4 == 0L])))
}

test_that("Gaussian mixed CV agrees with independent dense conditioning", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  f <- mixed_cv_fixture()
  d <- f$data
  formula <- y ~ I(scale(x)) + offset(log(exposure)) + (1 + x | group)
  for (target in c("new_group", "existing_group")) {
    folds <- if (target == "new_group") f$new else f$existing
    cv <- cross_validate_mixed(formula, d, folds, prediction = target, keep_fits = TRUE)
    expect_true(cv$complete, info = paste(unlist(cv$errors), collapse = " "))
    fit <- cv$fits[[1]]
    split <- folds[[1]]
    tr <- d[split$analysis, ]
    te <- d[split$assessment, ]
    co <- glmmTMB::fixef(fit)$cond
    C <- glmmTMB::VarCorr(fit)$cond$group
    sigma2 <- stats::sigma(fit)^2
    mu <- function(v) co[1] + co[2] * (v$x - mean(tr$x)) / sd(tr$x) + log(v$exposure)
    score <- vapply(unique(as.character(te$group)), function(g) {
      a <- tr[tr$group == g, ]
      t <- te[te$group == g, ]
      zt <- cbind(1, t$x)
      S <- zt %*% C %*% t(zt) + diag(sigma2, nrow(t))
      expected <- mu(t)
      if (nrow(a)) {
        za <- cbind(1, a$x)
        V <- za %*% C %*% t(za) + diag(sigma2, nrow(a))
        cross <- zt %*% C %*% t(za)
        expected <- expected + drop(cross %*% solve(V, a$y - mu(a)))
        S <- S - cross %*% solve(V, t(cross))
      }
      residual <- t$y - expected
      expect_equal(cv$predictions$prediction[match(rownames(t), cv$predictions$row_id)],
        unname(expected), tolerance = 1e-9)
      as.numeric((nrow(t) * log(2 * pi) + determinant(S, logarithm = TRUE)$modulus +
        crossprod(residual, solve(S, residual))) / 2)
    }, numeric(1))
    actual <- cv$joint_scores[[1]]
    expect_equal(actual$neg_log_score, unname(score[match(actual$group, unique(as.character(te$group)))]),
      tolerance = 1e-9)
    expect_equal(cv$selection_records[[1]]$K, attr(stats::logLik(fit), "df"))
  }
  expect_error(cross_validate_mixed(formula, d, f$existing, prediction = "new_group"),
    "does not match")
  expect_error(cross_validate_mixed(formula, d, f$new, prediction = "existing_group"),
    "does not match")
})

test_that("non-Gaussian CV matches direct integration over the group effect", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  f <- mixed_cv_fixture()
  d <- f$data
  forms <- list(poisson = poisson ~ x + (1 | group),
    nbinom2 = nb ~ x + (1 | group),
    binomial = cbind(successes, failures) ~ x + (1 | group),
    Gamma = gamma ~ x + (1 | group))
  for (family in names(forms)) for (target in c("new_group", "existing_group")) {
    folds <- if (target == "new_group") f$new else f$existing
    cv <- cross_validate_mixed(forms[[family]], d, folds, family, target, keep_fits = TRUE)
    expect_true(cv$complete, info = paste(family, target, paste(unlist(cv$errors), collapse = " ")))
    model <- cv$fits[[1]]
    b <- glmmTMB::fixef(model)$cond
    tau <- sqrt(glmmTMB::VarCorr(model)$cond$group[1, 1])
    phi <- stats::sigma(model)
    group <- cv$joint_scores[[1]]$group[1]
    te <- cv$predictions[cv$predictions$group == group, ]
    tr <- cv$targets[match(folds[[1]]$analysis, cv$targets$row_id), ]
    tr <- tr[tr$group == group, ]
    link <- if (family == "binomial") stats::plogis else exp
    eta <- function(rows) b[1] + b[2] * d[rows$row_id, "x"]
    dens <- function(rows, z) {
      mu <- link(eta(rows) + z)
      y <- rows$observed
      prod(switch(family, poisson = dpois(y, mu), binomial = dbinom(y, rows$trials, mu),
        nbinom2 = dnbinom(y, size = phi, mu = mu),
        Gamma = dgamma(y, shape = 1 / phi^2, scale = mu * phi^2)))
    }
    integrate_group <- function(fn) stats::integrate(function(z) vapply(z,
      function(u) fn(u) * stats::dnorm(u, sd = tau), numeric(1)),
      -10 * tau, 10 * tau, rel.tol = 1e-9, abs.tol = 0)$value
    training <- integrate_group(function(z) dens(tr, z))
    joint <- integrate_group(function(z) dens(tr, z) * dens(te, z))
    expected <- vapply(eta(te), function(e) integrate_group(function(z)
      link(e + z) * dens(tr, z)) / training, numeric(1))
    expect_equal(cv$joint_scores[[1]]$neg_log_score[1], log(training) - log(joint), tolerance = 1e-7)
    expect_equal(te$prediction, unname(expected), tolerance = 1e-7)
  }
})

test_that("mixed CV comparison, selection and scoring use training outcomes only", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  f <- mixed_cv_fixture()
  d <- f$data
  fixed <- y ~ I(scale(x)) + offset(log(exposure)) + (1 | group)
  null <- y ~ (1 | group)
  full <- cross_validate_mixed(fixed, d, f$new, keep_fits = TRUE)
  base <- cross_validate_mixed(null, d, f$new)
  hold <- f$new[[1]]$assessment
  d[hold, "y"] <- d[hold, "y"] + 20
  perturbed <- cross_validate_mixed(fixed, d, f$new)
  expect_equal(full$training, perturbed$training)
  expect_equal(full$predictions$prediction, perturbed$predictions$prediction)
  expect_false(isTRUE(all.equal(full$joint_scores, perturbed$joint_scores)))
  expect_error(cv_tab(list(base = base, changed = perturbed), "neg_log_score"),
    "same observation measure")
  expect_equal(nrow(cv_tab(list(base = base, full = full), "neg_log_score", by_fold = TRUE)), 2L)
  selected <- select_cv(list(base = base, full = full), "AIC")
  expect_true(selected$complete)
  aic <- selected$selection$tables[[1]]
  expect_equal(aic$AIC, -2 * aic$logLik + 2 * aic$K)
  expect_identical(unname(selected$selection$chosen), aic$Modnames[which.min(aic$AIC)])
  changed <- full
  changed$joint_scores[[1]]$neg_log_score[1] <- 0
  expect_error(cv_tab(list(changed = changed), "neg_log_score"), "changed")
  d[hold, "x"] <- NA_real_
  failed <- cross_validate_mixed(fixed, d, f$new)
  expect_false(failed$complete)
  expect_true(all(is.na(failed$summary$estimate)))
  expect_true(all(!failed$predictions$eligible))
  expect_identical(failed$errors[[1]]$stage, "group_prediction")
})

test_that("mixed CV refits ecological scales within the training split", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  set.seed(60)
  n <- 72L
  ids <- paste0("r", seq_len(n))
  distance <- seq(.03, .98, length.out = 20)
  radial <- outer(rnorm(n), sin(8 * distance)) + outer(rnorm(n), distance) +
    matrix(rnorm(n * length(distance), sd = .2), n)
  buffers <- stats::setNames(lapply(seq_len(n), function(i)
    matrix(radial[i, ], ncol = 1L, dimnames = list(NULL, "cover"))), ids)
  distances <- stats::setNames(rep(list(distance), n), ids)
  spec <- msr_vars(habitat = kernel_var("cover"))
  inputs <- structure(list(kernel_dat = data.frame(habitat = rep(999, n), row.names = ids),
    d_list = distances, raw_cov = buffers, kernel = "gaussian", shape = NULL,
    min_D = 40, max_D = 800, unit_conv = 800, n_covs = 1L, sigma = .25,
    scale_vars = spec, resolution = 20, n_cols = 40L, binned = NULL,
    scl_params = list(mean = c(habitat = 999), sd = c(habitat = 100)), cell_data_stored = TRUE),
    class = "multiScaleR_data")
  group <- factor(rep(letters[1:12], each = 6))
  weight <- exp(-0.5 * (distance / .32)^2)
  true <- drop(radial %*% (weight / sum(weight)))
  d <- data.frame(group, y = 1 + 3 * true + rnorm(12, sd = .5)[group] + rnorm(n, sd = .2),
    row.names = ids)
  fold <- list(one = list(analysis = ids[seq_len(n) %% 3 != 0L],
    assessment = ids[seq_len(n) %% 3 == 0L]))
  misaligned <- inputs
  misaligned$raw_cov <- rev(misaligned$raw_cov)
  expect_error(cross_validate_mixed(y ~ habitat + (1 | group), d, fold,
    prediction = "existing_group", kernel_inputs = misaligned), "matching raw-cell buffers")
  cv <- cross_validate_mixed(y ~ habitat + (1 | group), d, fold, prediction = "existing_group",
    kernel_inputs = inputs, keep_fits = TRUE)
  expect_true(cv$complete, info = paste(unlist(cv$errors), collapse = " "))
  expect_equal(cv$selection_records$one$ecological_K, 1)
  expect_true(is.finite(cv$training$one$scale_est[1, "Mean"]))
  expect_false(isTRUE(all.equal(unname(cv$training$one$normalization$mean), 999)))
  d[fold$one$assessment, "y"] <- d[fold$one$assessment, "y"] + 20
  alternative <- cross_validate_mixed(y ~ habitat + (1 | group), d, fold,
    prediction = "existing_group", kernel_inputs = inputs)
  expect_equal(cv$training, alternative$training)
  expect_equal(cv$predictions$prediction, alternative$predictions$prediction)
})
