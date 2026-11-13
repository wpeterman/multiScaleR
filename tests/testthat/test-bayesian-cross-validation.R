bayes_cv_fixture <- function() {
  f <- bayesian_scale_fixture(n_site = 5L, n_year = 4L)
  b <- f$inputs$binned
  ids <- b$point_ids
  values <- lapply(seq_along(ids), function(i) {
    keep <- b$csum$habitat[i, ] > 0
    cbind(cover = rep(b$vsum$habitat[i, keep] / b$csum$habitat[i, keep], b$csum$habitat[i, keep]))
  })
  distances <- lapply(seq_along(ids), function(i)
    rep(b$bin_rep * (1 - (i %% 5L) * .015), b$csum$habitat[i, ]))
  names(values) <- names(distances) <- ids
  f$inputs$raw_cov <- values
  f$inputs$d_list <- distances
  f$inputs$kernel_dat <- data.frame(habitat = rep(99, length(ids)), row.names = ids)
  f$inputs$scale_vars <- msr_vars(habitat = kernel_var("cover"))
  f$inputs$binned$nbins <- 16L
  d <- rbind(f$data, f$data[c(2L, 20L), , drop = FALSE])
  d$count[21:22] <- c(4L, 1L)
  rownames(d) <- paste0("obs", seq_len(nrow(d)))
  f$data <- d
  f$folds <- list(holdout = list(analysis = rownames(d)[d$year != "y4" & d$site != "s5"],
    assessment = rownames(d)[d$year == "y4"]))
  f
}

test_that("Bayesian CV refits training rows and scores joint new-group predictions", {
  f <- bayes_cv_fixture()
  d <- f$data
  args <- list(formula = count ~ habitat + poly(temp, 2) + (1 | site) + offset(log(effort)),
    data = d, kernel_inputs = f$inputs, folds = f$folds,
    scales = list(habitat = ~1 + (1 | year)), family = "gaussian", row_id = "key",
    seed = 24L, fit_args = list(chains = 2L, iter = 1500L, warmup = 500L), keep_fits = TRUE)
  set.seed(19)
  before <- .Random.seed
  cv <- do.call(cross_validate_bayes, args)
  expect_identical(.Random.seed, before)
  expect_s3_class(cv, "multiScaleR_bayes_cv")
  expect_true(all(is.finite(cv$joint_scores$neg_log_score)), info = paste(unlist(cv$errors), collapse = " "))
  fit <- cv$fits$holdout
  a <- d[f$folds$holdout$analysis, , drop = FALSE]
  b <- d[f$folds$holdout$assessment, , drop = FALSE]
  expect_identical(cv$predictions$row_id, rownames(b))
  expect_identical(cv$predictions$point_id, b$key)
  expect_identical(fit$specification$y, as.numeric(a$count))
  expect_identical(fit$specification$observation_ids, a$key)
  expect_true(all(cv$predictions$new_scale_level))
  expect_identical(cv$predictions$new_response_level, as.character(b$site) == "s5")
  expect_equal(cv$joint_scores$n, nrow(b))
  expect_identical(cv$prediction_contracts$holdout$new_groups$level, c("y4", "s5"))

  # Independent bin construction averages raw distances only across unique
  # analysis points. Fixed breaks come from the declared extraction domain.
  cuts <- seq(0, 1 + 1e-9, length.out = 17)
  train_d <- unlist(f$inputs$d_list[unique(a$key)], use.names = FALSE)
  train_code <- .bincode(train_d, cuts, include.lowest = TRUE)
  rep_d <- vapply(seq_len(16L), function(j) {
    keep <- train_code == j
    if (any(keep)) mean(train_d[keep]) else mean(cuts[j + 0:1])
  }, numeric(1))
  expect_equal(cv$training$holdout$bin_distances, rep_d, tolerance = 1e-15)
  point <- unique(b$key)
  sums <- counts <- matrix(0, length(point), length(rep_d))
  for (i in seq_along(point)) {
    id <- point[i]
    code <- .bincode(f$inputs$d_list[[id]], cuts, include.lowest = TRUE)
    for (j in seq_along(rep_d)) {
      use <- code == j
      counts[i, j] <- sum(use)
      sums[i, j] <- sum(f$inputs$raw_cov[[id]][use, "cover"])
    }
  }
  held <- list(binned = list(vsum = list(habitat = sums), csum = list(habitat = counts),
    bin_rep = rep_d, point_ids = point), kernel = "gaussian", unit_conv = 1000,
    min_D = 40, max_D = 1000)
  eta <- predict(fit, b, kernel_inputs = held, row_id = "key", type = "link",
    summary = FALSE, new_response_levels = "sample", new_scale_levels = "sample",
    seed = cv$training$holdout$score_seed)
  ll <- vapply(seq_len(nrow(b)), function(j)
    dnorm(b$count[j], eta[, j], exp(fit$draws$log_dispersion[, 1]), log = TRUE), numeric(nrow(eta)))
  expect_equal(cv$predictions$prediction, colMeans(eta), tolerance = 1e-12)
  expect_equal(cv$predictions$marginal_neg_log_score, -log(colMeans(exp(ll))), tolerance = 1e-12)
  expect_equal(cv$joint_scores$neg_log_score, -log(mean(exp(rowSums(ll)))), tolerance = 1e-12)
  expect_gt(abs(cv$joint_scores$neg_log_score - sum(cv$predictions$marginal_neg_log_score)), 1e-4)

  # Held-out outcomes and raster values cannot enter analysis draws or bins.
  d[rownames(b), "count"] <- b$count + 3
  args$data <- d
  changed <- do.call(cross_validate_bayes, args)
  expect_identical(changed$training, cv$training)
  expect_identical(changed$fits$holdout$draws, fit$draws)
  expect_identical(changed$predictions$prediction, cv$predictions$prediction)
  expect_false(isTRUE(all.equal(changed$joint_scores$neg_log_score, cv$joint_scores$neg_log_score)))
  args$kernel_inputs$d_list[[point[1L]]] <- args$kernel_inputs$d_list[[point[1L]]] / 2
  args$new_scale_levels <- "error"
  rejected <- do.call(cross_validate_bayes, args)
  expect_identical(rejected$training, cv$training)
  expect_identical(rejected$fits$holdout$draws, fit$draws)
  expect_identical(rejected$fold_results$status, "failed")
  expect_identical(rejected$errors$holdout$stage, "prediction")
  expect_match(rejected$errors$holdout$message, "New scale group levels")
  expect_false(rejected$complete)
  expect_true(is.na(rejected$summary$estimate))
  expect_error(cv_tab(list(failed = rejected), "neg_log_score"), "complete Bayesian")
  if (cv$complete) {
    expect_equal(cv_tab(list(candidate = cv), "neg_log_score")$Score, cv$summary$estimate)
    expect_error(cv_tab(list(first = cv, second = changed), "neg_log_score"), "identical splits")
    altered <- cv
    altered$joint_scores$neg_log_score <- 0
    expect_error(cv_tab(list(altered = altered), "neg_log_score"), "changed")
  } else {
    expect_identical(cv$fold_results$status, "diagnostic_failed")
    expect_false(cv$diagnostics$holdout$passed)
  }
  expect_error(select_cv(list(candidate = cv)), "no training ML")
})

test_that("Bayesian grouped count scores use normalized importance weights directly", {
  f <- bayesian_scale_fixture(n_site = 5L, n_year = 3L)
  d <- f$data
  d$success <- (seq_len(nrow(d)) - 1L) %% 3L
  d$failure <- 3L - d$success
  m <- multiScale_bayes(cbind(success, failure) ~ habitat, f$inputs, d,
    family = "binomial", row_id = "key", method = "prepare")
  # Three explicit posterior vectors with zero mass on a deliberately invalid
  # candidate. The nonzero vectors differ only in their response intercept.
  v <- unlist(m$parameters, use.names = FALSE)
  samples <- matrix(rep(v, each = 3L), nrow = 3L)
  samples[, 1L] <- c(-1, 2, NA)
  m$importance <- list(candidates = samples, weights = c(.25, .75, 0))
  x <- d[1:4, , drop = FALSE]
  context <- .msr_bayes_prediction_prepare(m, x, f$inputs, "key", "sample", "sample", TRUE, TRUE)
  target <- data.frame(observed = x$success, trials = x$success + x$failure,
    assessment_group = c("a", "b", "a", "b"))
  score <- .msr_bayes_cv_scores(m, context, target, 5L, 12L)
  p <- plogis(c(-1, 2))
  density <- vapply(p, function(pj) dbinom(target$observed, target$trials, pj), numeric(4L))
  expect_equal(score$prediction, rep(sum(c(.25, .75) * p), 4L))
  expect_equal(score$marginal_neg_log_score, -log(drop(density %*% c(.25, .75))))
  expect_equal(score$joint, -log(vapply(list(c(1L, 3L), c(2L, 4L)), function(i)
    sum(apply(density[i, , drop = FALSE], 2L, prod) * c(.25, .75)), numeric(1))))
  expect_identical(score$posterior_vectors, 2L)
  expect_identical(score$new_group_draws, 1L)
})

test_that("Bayesian CV checks fold and prepared point identity before fitting", {
  f <- bayes_cv_fixture()
  args <- list(formula = count ~ habitat, data = f$data, kernel_inputs = f$inputs,
    folds = f$folds, row_id = "key", family = "poisson")
  args$kernel_inputs$raw_cov <- NULL
  args$kernel_inputs$d_list <- NULL
  expect_error(do.call(cross_validate_bayes, args), "raw cells")
  args$kernel_inputs <- f$inputs
  args$data$key[1L] <- "unknown"
  expect_error(do.call(cross_validate_bayes, args), "prepared point ID")
  args$data <- f$data
  args$folds$holdout$assessment <- c(args$folds$holdout$assessment, args$folds$holdout$analysis[1L])
  expect_error(do.call(cross_validate_bayes, args), "overlapping")
})

test_that("Bayesian CV rebuilds fixed metrics and retains weighted surface buffers", {
  skip_if_not_installed("posterior", minimum_version = "1.6.0")
  f <- bayes_cv_fixture()
  f$inputs$scale_vars <- msr_vars(habitat = kernel_var("cover"),
    weighted = surface_var("cover", "sq", weighted = TRUE),
    fixed = surface_var("cover", "sq", radius = 240))
  d <- f$data
  # A short run checks the actual held-out computations and failed chain
  # diagnostic records without claiming it supplies a usable posterior.
  cv <- cross_validate_bayes(count ~ habitat + weighted + fixed, d, f$inputs, f$folds,
    row_id = "key", family = "gaussian", seed = 15L, keep_fits = TRUE,
    fit_args = list(chains = 2L, iter = 24L, warmup = 10L))
  expect_identical(cv$fold_results$status, "diagnostic_failed",
    info = paste(unlist(cv$errors), collapse = " "))
  expect_true(is.finite(cv$joint_scores$neg_log_score))
  expect_identical(names(cv$fits$holdout$specification$surface$buffers), "weighted")
  keys <- d[f$folds$holdout$analysis, "key"]
  raw <- vapply(keys, function(id)
    sd(f$inputs$raw_cov[[id]][f$inputs$d_list[[id]] <= .24, "cover"]), numeric(1))
  expect_equal(unname(cv$training$holdout$fixed_center), mean(raw), tolerance = 1e-12)
  expect_equal(unname(cv$training$holdout$fixed_spread), sd(raw), tolerance = 1e-12)
  expect_false(cv$diagnostics$holdout$passed)
  expect_true(is.na(cv$summary$estimate))
  expect_true(all(!cv$predictions$eligible))
})

test_that("Bayesian CV compares Poisson and NB2 mass on common count blocks", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 6L, n_year = 4L, annual = FALSE)
  d <- f$data
  folds <- list(last = list(analysis = rownames(d)[d$year != "y4"],
    assessment = rownames(d)[d$year == "y4"]))
  candidates <- lapply(c(poisson = "poisson", nb = "nbinom2", gaussian = "gaussian"),
    function(family) cross_validate_bayes(count ~ habitat, d, f$inputs, folds,
      row_id = "key", family = family, bin_reference = "fixed", seed = 24L,
      method = "laplace", fit_args = list(laplace_draws = 256L)))
  expect_true(all(vapply(candidates, function(x) x$complete, logical(1))),
    info = paste(unlist(lapply(candidates, function(x) x$errors)), collapse = " "))
  compared <- cv_tab(candidates[c("poisson", "nb")], "neg_log_score")
  expect_setequal(compared$Modnames, c("poisson", "nb"))
  expect_equal(compared$Score[match(c("poisson", "nb"), compared$Modnames)],
    vapply(candidates[c("poisson", "nb")], function(x) x$summary$estimate, numeric(1)),
    ignore_attr = TRUE)
  expect_error(cv_tab(candidates, "neg_log_score"), "identical splits")
})

test_that("Bayesian CV averages replicated new effects in the weighted joint density", {
  f <- bayesian_scale_fixture(n_site = 4L, n_year = 3L)
  d <- f$data
  m <- multiScale_bayes(count ~ habitat + (1 | site), f$inputs, d,
    family = "gaussian", row_id = "key", method = "prepare")
  beta <- c(-.5, .3); tau <- c(.2, 1.25); sigma <- c(.6, .9); weight <- c(.2, .8)
  p <- m$parameters
  samples <- matrix(rep(unlist(p, use.names = FALSE), each = 2L), nrow = 2L)
  block <- rep(names(p), lengths(p))
  samples[, which(block == "beta")[1L]] <- beta
  samples[, block == "log_response_sd"] <- log(tau)
  samples[, block == "log_dispersion"] <- log(sigma)
  m$importance <- list(candidates = samples, weights = weight)
  x <- d[1:3, , drop = FALSE]
  x$site <- c("new_a", "new_a", "new_b")
  context <- .msr_bayes_prediction_prepare(m, x, f$inputs, "key", "sample", "sample", TRUE, TRUE)
  target <- data.frame(observed = c(.1, .4, -.7), trials = 1, assessment_group = "all")
  seed <- 914L; repeats <- 4L
  score <- .msr_bayes_cv_scores(m, context, target, repeats, seed)

  # Draw the two new effects explicitly for every mixture component and
  # replicate. The first two observations share the first normal draw.
  set.seed(seed)
  density <- rep(0, length(weight)); marginal <- matrix(0, 3L, length(weight))
  for (i in seq_along(weight)) for (j in seq_len(repeats)) {
    b <- rnorm(2L, sd = tau[i])
    p_row <- dnorm(target$observed, beta[i] + b[c(1L, 1L, 2L)], sigma[i])
    density[i] <- density[i] + prod(p_row) / repeats
    marginal[, i] <- marginal[, i] + p_row / repeats
  }
  expect_identical(score$new_group_draws, repeats)
  expect_equal(score$joint, -log(sum(weight * density)), tolerance = 1e-12)
  expect_equal(score$marginal_neg_log_score, -log(drop(marginal %*% weight)), tolerance = 1e-12)
})
