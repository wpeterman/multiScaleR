# One deterministic software fixture; at most24 planned scale fits per lifecycle.
# Native fixed-scale GLMs are score/design oracles, not a recovery study.
cv_stage1_cache <- new.env(parent = emptyenv())
cv_stage1_inputs <- function(kernel = "gaussian", sq = FALSE) {
  set.seed(9117)
  n <- 72L; ids <- paste0("cv", seq_len(n)); d <- seq(.025, .99, length.out = 24L)
  radial <- outer(rnorm(n), sin(7 * d)) + outer(rnorm(n), cos(13 * d)) +
    outer(rnorm(n), d) + matrix(rnorm(n * length(d), sd = .2), n)
  specs <- if (sq) msr_vars(hab = surface_var("terrain", metric = "sq", weighted = TRUE)) else
    msr_vars(hab = kernel_var("terrain"))
  buffers <- lapply(seq_len(n), function(i) matrix(radial[i, ], ncol = 1L,
    dimnames = list(NULL, "terrain")))
  distances <- lapply(seq_len(n), function(i) d * (1 - (i %% 5L) * .005))
  names(buffers) <- names(distances) <- ids
  inputs <- structure(list(kernel_dat = data.frame(hab = rep(0, n), row.names = ids),
    d_list = distances, raw_cov = buffers, kernel = kernel, shape = if (kernel == "expow") 1.8 else NULL,
    min_D = 35, max_D = 800, unit_conv = 800, n_covs = 1L, sigma = .3,
    scale_vars = specs, resolution = 20, n_cols = 40L,
    binned = if (sq) NULL else .msr_build_kernel_bins(distances, buffers, "hab", "terrain", 32L, ids),
    scl_params = list(mean = c(hab = 999), sd = c(hab = 100)), cell_data_stored = TRUE),
    class = "multiScaleR_data")
  raw <- vapply(seq_len(n), function(i) .msr_eval_scale_vars(distances[[i]], buffers[[i]],
    specs, .3, inputs$shape, kernel, 800, 20, 40L)[1L], 0)
  hab <- as.numeric(scale(raw))
  inputs$kernel_dat$hab <- hab
  inputs$scl_params <- list(mean = c(hab = mean(raw)), sd = c(hab = sd(raw)))
  dat <- data.frame(hab = rep(999, n), z = rnorm(n),
    group = factor(rep(c("a", "b"), n / 2L)), exposure = runif(n, .8, 1.2), row.names = ids)
  gaussian_eta <- 1 + 1.3 * hab + .2 * dat$z + .25 * (dat$group == "b") +
    .15 * hab * (dat$group == "b") + log(dat$exposure)
  dat$continuous <- gaussian_eta + rnorm(n, sd = .25)
  count_eta <- .7 + .7 * hab + .15 * dat$z^2 + .1 * (dat$group == "b") + log(dat$exposure)
  dat$count <- rpois(n, exp(count_eta))
  p <- plogis(-.3 + .8 * hab + .2 * dat$z + .1 * (dat$group == "b") + log(dat$exposure))
  dat$binary <- as.logical(rbinom(n, 1L, p))
  trials <- rep(c(4L, 6L, 5L), length.out = n)
  dat$yes <- rbinom(n, trials, p); dat$no <- trials - dat$yes
  folds <- stats::setNames(lapply(0:2, function(k) list(
    analysis = ids[seq_len(n) %% 3L != k], assessment = ids[seq_len(n) %% 3L == k])), c("a", "b", "c"))
  list(data = dat, inputs = inputs, folds = folds)
}

cv_stage1_fixture <- function() {
  if (exists("failure", cv_stage1_cache, inherits = FALSE)) stop(cv_stage1_cache$failure)
  if (!exists("value", cv_stage1_cache, inherits = FALSE)) {
    cv_stage1_cache$value <- tryCatch({
      f <- cv_stage1_inputs(); planned <- 0L
      run <- function(formula, data = f$data, inputs = f$inputs, folds = f$folds,
                      family = "gaussian", backend = "profile_glm", weights = NULL, keep_fits = TRUE) {
        planned <<- planned + length(folds); stopifnot(planned <= 24L)
        cross_validate(formula, data, inputs, folds, family = family, backend = backend,
          weights = weights, keep_fits = keep_fits)
      }
      formulas <- list(gaussian = continuous ~ hab * group + poly(z, 2) + offset(log(exposure)),
        poisson = count ~ hab + I(z^2) + group + offset(log(exposure)),
        binomial = cbind(yes, no) ~ hab + z * group + offset(log(exposure)))
      w <- rep(c(1, 2, 3), length.out = nrow(f$data))
      fits <- list(gaussian = run(formulas$gaussian),
        poisson = run(formulas$poisson, family = "poisson", weights = w),
        binomial = run(formulas$binomial, family = "binomial", weights = w),
        logical = run(I(binary) ~ hab + z + group, family = "binomial", folds = f$folds[1L]))
      changed <- f$data
      held <- f$folds$a$assessment
      changed[held, "continuous"] <- changed[held, "continuous"] + 100
      changed[held, "z"] <- changed[held, "z"] + 50
      changed$hab <- -1e5
      changed_inputs <- f$inputs
      changed_inputs$kernel_dat$hab <- 1e6
      changed_inputs$scl_params <- list(mean = c(hab = -1e9), sd = c(hab = .001))
      changed_inputs$binned$bin_rep[] <- 1e-5
      for (id in held) {
        changed_inputs$raw_cov[[id]] <- changed_inputs$raw_cov[[id]] * 10
        changed_inputs$d_list[[id]] <- changed_inputs$d_list[[id]] * .7
      }
      leak <- run(formulas$gaussian, changed, changed_inputs, folds = f$folds[1L], keep_fits = FALSE)
      refit <- run(formulas$gaussian, folds = f$folds[1L], backend = "refit")
      other <- list()
      for (kernel in c("exp", "expow")) {
        k <- cv_stage1_inputs(kernel)
        other[[kernel]] <- run(formulas$gaussian, k$data, k$inputs, k$folds[1L])
      }
      s <- cv_stage1_inputs(sq = TRUE)
      sq <- run(continuous ~ hab + z + group, s$data, s$inputs, s$folds[1L])
      rare <- f$data
      rare$group <- factor(c(as.character(rare$group[1:48]), rep("c", 24L)))
      ids <- rownames(rare)
      failed <- run(formulas$gaussian, rare, folds = list(
        good = list(analysis = ids[13:48], assessment = ids[1:12]),
        bad = list(analysis = ids[1:48], assessment = ids[49:72])), keep_fits = FALSE)
      m <- cv_stage1_inputs()
      for (i in seq_along(m$inputs$raw_cov)) {
        terrain <- m$inputs$raw_cov[[i]][, "terrain"]
        canopy <- .4 * terrain^2 + sin(m$inputs$d_list[[i]] * 15) * cos(i * .2)
        m$inputs$raw_cov[[i]] <- cbind(terrain = terrain, canopy = canopy)
      }
      m$inputs$scale_vars <- msr_vars(hab = kernel_var("terrain"), canopy = kernel_var("canopy"))
      m$inputs$sigma <- c(.3, .4); m$inputs$n_covs <- 2L
      m$inputs$kernel_dat$canopy <- 0
      raw_canopy <- vapply(seq_along(m$inputs$d_list), function(i)
        scale_type(d = m$inputs$d_list[[i]], kernel = "gaussian", sigma = .45,
          r_stack.df = m$inputs$raw_cov[[i]][, "canopy", drop = FALSE])[1L], 0)
      m$data$continuous <- m$data$continuous + .7 * as.numeric(scale(raw_canopy))
      multiple <- run(continuous ~ hab + canopy + z + group, m$data, m$inputs, m$folds[1L])
      # Six fixed-map metadata cases reuse the existing raw-cell fixtures.
      # No new draws, seeds or study panels; the original lifecycle cap is24.
      tag <- function(inputs) {
        inputs$raster_group <- stats::setNames(rep(c("mapA", "mapB"), length.out = nrow(inputs$kernel_dat)),
          rownames(inputs$kernel_dat))
        inputs
      }
      grouped <- list(
        gaussian = run(formulas$gaussian, inputs = tag(f$inputs), folds = f$folds[1L]),
        poisson = run(formulas$poisson, inputs = tag(f$inputs), folds = f$folds[1L], family = "poisson", weights = w),
        binomial = run(formulas$binomial, inputs = tag(f$inputs), folds = f$folds[1L], family = "binomial", weights = w),
        refit = run(formulas$gaussian, inputs = tag(f$inputs), folds = f$folds[1L], backend = "refit"),
        sq = run(continuous ~ hab + z + group, s$data, tag(s$inputs), s$folds[1L]),
        expow = run(formulas$gaussian, k$data, tag(k$inputs), k$folds[1L]))
      list(fixture = f, formulas = formulas, weights = w, fits = fits, leak = leak,
        refit = refit, other = other, sq = sq, failed = failed, multiple = multiple, grouped = grouped, planned_scale_fits = planned)
    }, error = function(e) { cv_stage1_cache$failure <- e; stop(e) })
  }
  cv_stage1_cache$value
}

test_that("explicit-fold GLM CV refits scales and preserves every planned row", {
  z <- cv_stage1_fixture(); f <- z$fixture
  expect_lte(z$planned_scale_fits, 24L)
  for (name in names(z$fits)) {
    x <- z$fits[[name]]
    expect_s3_class(x, "multiScaleR_cv")
    expect_true(x$complete, info = paste(name, paste(unlist(x$errors), collapse = "; ")))
    expect_true(all(x$fold_results$status == "scored"))
    expect_equal(nrow(x$predictions), sum(lengths(lapply(x$folds, `[[`, "assessment"))))
    expect_true(all(x$summary$complete))
    expect_true(all(is.finite(x$summary$estimate)))
    for (fold in names(x$folds)) {
      expect_identical(rownames(x$fits[[fold]]$opt_mod$model), x$folds[[fold]]$analysis)
      expect_identical(x$diagnostics[[fold]]$scale, x$fits[[fold]]$diagnostics)
      expect_identical(x$training[[fold]]$coefficients, stats::coef(x$fits[[fold]]$opt_mod))
      expect_identical(x$training[[fold]]$scale_est, x$fits[[fold]]$scale_est)
      expect_identical(x$training[[fold]]$start_results, x$fits[[fold]]$profile_fit$start_results)
      expect_gte(nrow(x$training[[fold]]$start_results), 6L)
      expect_equal(length(x$training[[fold]]$bin_reference$bin_rep), 32L)
    }
  }
  expect_identical(unname(z$fits$gaussian$assessment_counts), rep(1L, nrow(f$data)))
  expect_null(z$leak$fits)
  expect_match(paste(capture.output(print(z$fits$gaussian)), collapse = " "), "3 of 3 folds scored")
  expect_identical(z$fits$gaussian$settings$parameter_uncertainty, "not_integrated")
  expect_identical(z$fits$gaussian$settings$calibration, "not_established")
})

test_that("CV predictions and count scores agree with fixed-scale native GLM oracles", {
  z <- cv_stage1_fixture(); f <- z$fixture
  for (name in c("gaussian", "poisson", "binomial")) for (fold in names(f$folds)) {
    x <- z$fits[[name]]; fit <- x$fits[[fold]]
    train <- f$data[x$folds[[fold]]$analysis, , drop = FALSE]
    engine <- fit$profile_fit$contract$objective$engine
    kp <- fit$profile_fit$par[seq_along(engine$parameter_names)]
    kd <- .msr_eval_kernel_engine(engine, kp)
    train$hab <- kd$value[, "hab"]
    w <- if (name == "gaussian") rep(1, nrow(train)) else z$weights[match(rownames(train), rownames(f$data))]
    family <- switch(name, gaussian = gaussian(), poisson = poisson(), binomial = binomial())
    native <- do.call(stats::glm, list(formula = z$formulas[[name]], data = train,
      family = family, weights = w, control = glm.control(epsilon = 1e-11, maxit = 100L)))
    expect_equal(unname(coef(fit$opt_mod)), unname(coef(native)), tolerance = 1e-6)
    held <- f$data[x$folds[[fold]]$assessment, , drop = FALSE]
    inputs <- .msr_cv_inputs(f$inputs, rownames(held), "hab", x$training[[fold]]$bin_reference)
    en <- .msr_prepare_kernel_engine(inputs$binned, "hab", kernel = "gaussian", standardize = FALSE)
    raw <- .msr_eval_kernel_engine(en, kp)$raw[, 1L]
    held$hab <- (raw - kd$center) / kd$spread
    expect_equal(unname(predict(native, held, type = "response")),
      x$predictions$prediction[x$predictions$fold == fold], tolerance = 1e-6)
  }
  for (name in c("poisson", "binomial")) {
    x <- z$fits[[name]]; p <- x$predictions
    score <- if (name == "poisson") -dpois(p$observed, p$prediction, log = TRUE) else
      -dbinom(p$observed, p$trials, p$prediction, log = TRUE)
    expect_equal(p$neg_log_score, score, tolerance = 0)
    expect_equal(x$summary$estimate[x$summary$metric == "neg_log_score"],
      weighted.mean(score, p$weight), tolerance = 1e-12)
  }
  p <- z$fits$binomial$predictions
  manual <- sum(p$weight * (p$observed * (1 - p$prediction)^2 +
    (p$trials - p$observed) * p$prediction^2)) / sum(p$weight * p$trials)
  expect_equal(z$fits$binomial$summary$estimate[z$fits$binomial$summary$metric == "brier"], manual, tolerance = 1e-12)
  p <- z$fits$gaussian$predictions
  expect_equal(z$fits$gaussian$summary$estimate[1L], sqrt(mean((p$observed - p$prediction)^2)), tolerance = 1e-12)
})

test_that("held-out outcomes, predictors and full-preparation summaries cannot influence training", {
  z <- cv_stage1_fixture(); base <- z$fits$gaussian
  expect_true(z$leak$complete, info = paste(unlist(z$leak$errors), collapse = "; "))
  # Wall-clock duration is retained, but is not a fitted quantity or a timing benchmark.
  left <- z$leak$training$a; right <- base$training$a
  expect_identical(left[setdiff(names(left), "start_results")], right[setdiff(names(right), "start_results")])
  fields <- setdiff(names(left$start_results), "seconds")
  expect_identical(left$start_results[fields], right$start_results[fields])
  expect_identical(z$leak$diagnostics$a, base$diagnostics$a)
  expect_false(identical(z$leak$predictions$prediction,
    base$predictions$prediction[base$predictions$fold == "a"]))
  for (name in names(z$other)) {
    expect_true(z$other[[name]]$complete, info = paste(name, paste(unlist(z$other[[name]]$errors), collapse = "; ")))
    expect_equal(z$other[[name]]$settings$backend, "profile_glm")
  }
  expect_true(z$sq$complete, info = paste(unlist(z$sq$errors), collapse = "; "))
  expect_true(z$refit$complete, info = paste(unlist(z$refit$errors), collapse = "; "))
  expect_true(z$multiple$complete, info = paste(unlist(z$multiple$errors), collapse = "; "))
  expect_identical(rownames(z$multiple$training$a$scale_est), c("hab", "canopy"))
  expect_equal(nrow(z$multiple$training$a$scale_est), 2L)
  expect_identical(z$refit$diagnostics$a$scale, z$refit$fits$a$diagnostics)
})

test_that("CV retains failed folds and never silently ranks a partial result", {
  x <- cv_stage1_fixture()$failed
  expect_false(x$complete)
  expect_identical(x$fold_results$status, c("scored", "failed"))
  expect_equal(nrow(x$predictions), 36L)
  expect_equal(sum(x$predictions$fold == "bad"), 24L)
  expect_true(all(is.na(x$predictions$prediction[x$predictions$fold == "bad"])))
  expect_true(all(is.na(x$summary$estimate)))
  expect_true(all(!x$summary$complete))
  expect_true(all(is.finite(x$metrics$estimate[x$metrics$fold == "good"])))
  expect_identical(x$errors$bad$stage, "prediction")
  expect_match(x$errors$bad$message, "new level")
  expect_false(is.null(x$training$bad$scale_est))
  expect_false(is.null(x$diagnostics$bad$scale))
  expect_match(paste(capture.output(print(x)), collapse = " "), "Incomplete")
})

test_that("frozen training bins and exact raw initialization preserve row identity", {
  f <- cv_stage1_inputs(); before <- serialize(f$inputs, NULL)
  a <- .msr_cv_inputs(f$inputs, f$folds$a$analysis, "hab")
  b <- .msr_cv_inputs(f$inputs, f$folds$a$assessment, "hab", a$binned$cv_reference)
  expect_identical(b$binned$cv_reference, a$binned$cv_reference)
  expect_identical(b$binned$bin_rep, a$binned$bin_rep)
  expect_identical(a$binned$point_ids, f$folds$a$analysis)
  expect_identical(rownames(.msr_cv_initial_values(b)), f$folds$a$assessment)
  rev <- .msr_cv_inputs(f$inputs, rev(f$folds$a$analysis), "hab")
  expect_equal(rev$binned$bin_rep, a$binned$bin_rep, tolerance = 1e-14)
  expect_equal(rev$binned$vsum$hab, a$binned$vsum$hab[rev(seq_len(nrow(a$kernel_dat))), ], tolerance = 0)
  reference <- a$binned$cv_reference
  for (j in seq_len(32L)) {
    values <- unlist(a$d_list, use.names = FALSE)
    code <- .bincode(values, reference$breaks, include.lowest = TRUE)
    expected <- if (any(code == j)) mean(values[code == j]) else
      mean(reference$breaks[j + 0:1])
    expect_equal(a$binned$bin_rep[j], expected, tolerance = 1e-14)
  }
  expect_identical(serialize(f$inputs, NULL), before)
})

test_that("GLM scoring preserves likelihood constants, trials and infinite losses", {
  p <- data.frame(observed = c(0, 2, 4), trials = 4, prediction = c(.2, .5, .8), weight = c(1, 2, 0))
  x <- .msr_cv_score(p, "binomial", c("neg_log_score", "brier"))
  expect_equal(x$neg_log_score[1:2], -dbinom(p$observed[1:2], 4, p$prediction[1:2], log = TRUE), tolerance = 0)
  expect_true(is.na(x$neg_log_score[3L]))
  expect_identical(x$eligible, c(TRUE, TRUE, FALSE))
  expect_equal(.msr_cv_metric(x, "neg_log_score", "binomial", TRUE, 3L)$weight_sum, 3)
  expect_equal(.msr_cv_metric(x, "brier", "binomial", TRUE, 3L)$weight_sum, 12)
  p <- data.frame(observed = c(1, 0), trials = 1, prediction = c(0, 1), weight = 1)
  x <- .msr_cv_score(p, "binomial", c("neg_log_score", "brier"))
  expect_true(all(is.infinite(x$neg_log_score)))
  expect_equal(.msr_cv_metric(x, "neg_log_score", "binomial", TRUE, 2L)$estimate, Inf)
  expect_equal(x$brier, c(1, 1))
  p <- data.frame(observed = c(0, 2), trials = 1, prediction = 0, weight = 1)
  x <- .msr_cv_score(p, "poisson", c("neg_log_score", "deviance"))
  expect_equal(x$neg_log_score[1L], 0)
  expect_true(is.infinite(x$neg_log_score[2L]))
  expect_equal(x$deviance[1L], 0)
  expect_error(.msr_cv_score(transform(p, weight = 0), "poisson", "neg_log_score"), "no positive-weight")
})

test_that("invalid fold and model contracts fail before any scale fit", {
  f <- cv_stage1_inputs(); ids <- rownames(f$data)
  call <- function(folds = f$folds, inputs = f$inputs, ...) cross_validate(
    continuous ~ hab + z, f$data, inputs, folds, ...)
  expect_error(call(folds = list(a = list(analysis = ids[1:10], assessment = ids[10:12]))), "overlapping")
  expect_error(call(folds = list(a = list(analysis = ids[1:10], assessment = "unknown"))), "unknown IDs")
  expect_error(call(folds = list(a = list(analysis = ids[c(1, 1)], assessment = ids[2]))), "duplicate")
  expect_error(call(folds = unname(f$folds)), "uniquely named")
  expect_error(call(optim_args = list(n_cores = 2)), "worker overrides")
  expect_error(call(optim_args = list(refit_fn = identity)), "worker overrides")
  expect_error(call(optim_args = list(profile_control = list(worker_pool = list(1)))), "worker overrides")
  lean <- f$inputs; lean$d_list <- lean$raw_cov <- NULL
  expect_error(call(inputs = lean), "store_cell_data")
  expect_error(call(weights = rep(.5, 72L)), "integer frequencies")
  expect_error(call(family = "Gamma"), "arg")
  expect_error(call(metrics = "brier"), "supported for this family")
  expect_error(call(glm_control = list(unknown = 1)), "glm.control")
  expect_error(cross_validate(continuous ~ hab + absent, f$data, f$inputs, f$folds), "every formula variable")
  expect_error(.msr_cv_outcome(continuous ~ hab, transform(f$data, continuous = NA_real_), "gaussian"), "finite")
  expect_error(.msr_cv_outcome(I(scale(continuous)) ~ hab, f$data, "gaussian"), "predefine")
})

# These checks reuse the existing optimized fixture; they add no scale fits.
test_that("fold generators produce explicit splits accepted by scale CV", {
  z <- cv_stage1_fixture(); f <- z$fixture; ids <- rownames(f$data)
  generated <- grouped_cv_folds(as.character(seq_along(ids) %% 3L), 3, ids)
  spatial <- spatial_cv_folds(cbind((seq_along(ids) %% 3L) * 100, 0), 50, 3, ids, crs = 32617)
  for (folds in list(generated, spatial)) {
    actual <- lapply(folds, function(s) sort(s$assessment))
    original <- lapply(f$folds, function(s) sort(s$assessment))
    expect_equal(unname(actual), unname(original), ignore_attr = TRUE)
    expect_identical(.msr_cv_folds(folds, ids), folds)
  }
  expect_identical(z$fits$poisson$targets$row_id, sort(ids))
  expect_equal(z$fits$poisson$targets$observed, f$data[sort(ids), "count"])
})

test_that("ordinary GLM references, CV tables and training selection match native oracles", {
  z <- cv_stage1_fixture(); f <- z$fixture
  null <- cross_validate_glm(count ~ 1, f$data, f$folds, "poisson", weights = z$weights, keep_fits = TRUE)
  trend <- cross_validate_glm(count ~ poly(z, 2) + group + offset(log(exposure)), f$data,
    f$folds, "poisson", weights = z$weights, keep_fits = TRUE)
  expect_true(null$complete); expect_true(trend$complete)
  models <- list(spatial = z$fits$poisson, null = null, trend = trend)
  tab <- cv_tab(models, "neg_log_score")
  expect_setequal(tab$Modnames, names(models)); expect_true(all(diff(tab$Score) >= 0))
  for (nm in names(models)) {
    x <- models[[nm]]; p <- x$predictions
    expected <- weighted.mean(-dpois(p$observed, p$prediction, log = TRUE), p$weight)
    expect_equal(tab$Score[tab$Modnames == nm], expected)
    expect_equal(cv_tab(models, "rmse")$Score[cv_tab(models, "rmse")$Modnames == nm],
      sqrt(weighted.mean((p$observed - p$prediction)^2, p$weight)))
  }
  per <- cv_tab(models, "neg_log_score", by_fold = TRUE)
  expect_equal(nrow(per), 9L)
  expect_setequal(per$Fold, names(f$folds))
  for (criterion in c("AIC", "AICc", "BIC")) {
    chosen <- select_cv(models, criterion, ties = "first", keep_fits = TRUE)
    expect_true(chosen$complete, info = paste(unlist(chosen$errors), collapse = "; "))
    expect_identical(chosen$settings$selection_uncertainty, "not_integrated")
    for (fold in names(f$folds)) {
      fits <- lapply(models, function(x) x$fits[[fold]])
      native <- if (criterion == "BIC") bic_tab(fits, mod_names = names(models)) else
        aic_tab(fits, AICc = criterion == "AICc", mod_names = names(models))
      expect_equal(chosen$selection$tables[[fold]], native)
      winner <- native$Modnames[which.min(native[[criterion]])]
      expect_identical(chosen$selection$chosen[[fold]], winner)
      expect_identical(chosen$fits[[fold]], models[[winner]]$fits[[fold]])
      expect_equal(chosen$predictions[chosen$predictions$fold == fold, ],
        models[[winner]]$predictions[models[[winner]]$predictions$fold == fold, ], ignore_attr = TRUE)
      # Native GLM count plus optimized scale/shape counts remain provider-owned.
      expect_equal(native$K[native$Modnames == "spatial"],
        attr(logLik(z$fits$poisson$fits[[fold]]$opt_mod), "df") + 1)
    }
    expect_true(all(is.finite(cv_tab(list(selected = chosen, null = null), "neg_log_score")$Score)))
  }
  tied <- select_cv(list(first = null, second = null), "AIC")
  expect_false(tied$complete); expect_true(all(is.na(tied$summary$estimate)))
  expect_match(tied$errors$a$message, "tie")
  expect_true(all(is.na(tied$selection$chosen)))
  fixed <- select_cv(list(first = null, second = null), "AIC", ties = "first")
  expect_true(fixed$complete); expect_true(all(fixed$selection$chosen == "first")); expect_null(fixed$fits)
  # Valid arithmetic copies alter assessment predictions, never training fits.
  changed <- models
  changed$trend$predictions$prediction <- rep(100, nrow(changed$trend$predictions))
  expect_identical(select_cv(changed, "AIC", ties = "first")$selection$chosen,
    select_cv(models, "AIC", ties = "first")$selection$chosen)
  expect_error(select_cv(list(compact = cross_validate_glm(count ~ 1, f$data, f$folds, "poisson", weights = z$weights))), "retain every")
  expect_error(aic_tab(list(null)), "use cv_tab")
  expect_error(bic_tab(list(null)), "use cv_tab")
  replaced <- null
  replaced$fits$a <- do.call(glm, list(formula = count ~ 1, data = f$data, family = poisson(), weights = z$weights))
  expect_error(select_cv(list(original = null, replaced = replaced)), "analysis row IDs")
})

test_that("CV comparison binds split plans, training-only targets, trials and weights", {
  z <- cv_stage1_fixture(); x <- z$fits$poisson
  variants <- list(x, x, x, x)
  variants[[1]]$folds$a$analysis <- rev(variants[[1]]$folds$a$assessment)
  variants[[2]]$targets$observed[1] <- variants[[2]]$targets$observed[1] + 1
  variants[[3]]$targets$weight[1] <- variants[[3]]$targets$weight[1] + 1
  variants[[4]]$targets$trials[1] <- 2
  for (other in variants) expect_error(cv_tab(list(original = x, other = other), "rmse"))
  incomplete <- z$failed
  expect_error(cv_tab(list(bad = incomplete), "rmse"), "complete")
  legacy <- x; legacy$targets <- NULL
  expect_equal(cv_tab(list(current = x, legacy = legacy), "rmse")$Score, rep(x$summary$estimate[1], 2))
  partial <- legacy
  partial$folds <- partial$folds["a"]; partial$predictions <- partial$predictions[partial$predictions$fold == "a", ]
  partial$fold_results <- partial$fold_results[1, ]
  expect_error(cv_tab(list(partial = partial), "rmse"), "Legacy CV targets")
  badfit <- x; badfit$fits$a$opt_mod$model$count[1] <- 999
  expect_error(select_cv(list(good = x, bad = badfit)), "training responses")
  expect_error(cv_tab(list(x, x), "rmse"), "uniquely named")
  expect_error(cv_tab(list(good = x), "AIC"), "supported metric")
  expect_error(cv_tab(list(binomial = z$fits$binomial, poisson = x), "neg_log_score"), "identical splits|same response family")
})

test_that("ordinary binomial GLM scores preserve trial and frequency definitions", {
  z <- cv_stage1_fixture(); f <- z$fixture
  x <- cross_validate_glm(cbind(yes, no) ~ z * group + offset(log(exposure)), f$data, f$folds,
    "binomial", weights = z$weights, keep_fits = TRUE)
  expect_true(x$complete)
  p <- x$predictions
  expect_equal(cv_tab(list(native = x), "neg_log_score")$Score,
    weighted.mean(-dbinom(p$observed, p$trials, p$prediction, log = TRUE), p$weight))
  expect_equal(cv_tab(list(native = x), "brier")$Score,
    sum(p$weight * (p$observed * (1 - p$prediction)^2 + (p$trials - p$observed) * p$prediction^2)) /
      sum(p$weight * p$trials))
  expect_true(select_cv(list(binomial = x), "AIC")$complete)
  binary <- cross_validate_glm(I(binary) ~ z + group, f$data, f$folds, "binomial")
  expect_true(binary$complete); expect_null(binary$fits)
  # Zero-probability events remain infinite, rather than being clipped.
  zero <- z$fits$poisson; zero$predictions$prediction[] <- 0
  tab <- cv_tab(list(zero = zero, other = zero), "neg_log_score")
  expect_true(all(is.infinite(tab$Score))); expect_true(all(is.na(tab$Delta)))
  expect_equal(tab$Rank, c(1, 1))
})

test_that("grouped raster CV preserves row maps and native fixed-map predictions", {
  z <- cv_stage1_fixture(); f <- z$fixture
  expect_identical(z$planned_scale_fits, 24L)
  references <- list(gaussian = z$fits$gaussian, poisson = z$fits$poisson, binomial = z$fits$binomial,
    refit = z$refit, sq = z$sq, expow = z$other$expow)
  for (name in names(z$grouped)) {
    x <- z$grouped[[name]]; ref <- references[[name]]
    expect_true(x$complete, info = paste(name, paste(unlist(x$errors), collapse = "; ")))
    p <- ref$predictions[ref$predictions$fold == "a", , drop = FALSE]
    expect_equal(x$predictions, p, ignore_attr = TRUE)
    expect_equal(x$training$a$coefficients, ref$training$a$coefficients, tolerance = 0)
    expect_equal(x$training$a$scale_est, ref$training$a$scale_est, tolerance = 0)
    expect_equal(x$training$a$fitted_normalization, ref$training$a$fitted_normalization, tolerance = 0)
    for (part in c("analysis", "assessment")) {
      expect_identical(names(x$raster_groups$a[[part]]), f$folds$a[[part]])
      expected <- stats::setNames(rep(c("mapA", "mapB"), length.out = nrow(f$data)), rownames(f$data))
      expect_identical(x$raster_groups$a[[part]], expected[f$folds$a[[part]]])
    }
    expect_identical(x$fits$a$kernel_inputs$raster_group, x$raster_groups$a$analysis)
    expect_identical(x$settings$parameter_uncertainty, "not_integrated")
  }
  x <- z$grouped$poisson
  selected <- select_cv(list(grouped = x), "AIC")
  expect_true(selected$complete)
  expect_identical(selected$raster_groups, x$raster_groups)
  expect_identical(selected$selection$candidate_settings$grouped, x$settings)
  expect_identical(selected$settings$controls_from, "grouped")
  expect_equal(selected$selection$tables$a$K, attr(logLik(x$fits$a$opt_mod), "df") + 1)
})

test_that("malformed grouped raster identities stop before scale fitting", {
  f <- cv_stage1_inputs(); ids <- rownames(f$data)
  inputs <- f$inputs
  inputs$raster_group <- stats::setNames(rep(c("A", "B"), length.out = length(ids)), ids)
  call <- function(x) cross_validate(continuous ~ hab + z, f$data, x, f$folds)
  x <- inputs; names(x$raster_group) <- rev(ids)
  expect_error(call(x), "raster_group.*align")
  x <- inputs; x$raster_group[1] <- NA_character_
  expect_error(call(x), "raster_group.*nonmissing")
  x <- inputs; x$raster_group[1] <- ""
  expect_error(call(x), "raster_group.*nonmissing")
  x <- inputs; x$raster_group <- unname(x$raster_group)
  expect_error(call(x), "raster_group.*align")
})
