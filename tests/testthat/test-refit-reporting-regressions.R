test_that("large profile payloads have bounded keys and data-specific caches", {
  object <- structure(list(scale_est = data.frame(Mean = 2, SE = .2, row.names = "x"),
    min_D = .1, max_D = 10, kernel_inputs = list(kernel = "gaussian", unit_conv = 1),
    optim_results = list(par = 2, value = 4),
    call = as.call(list(as.name("fit"), rep("long_payload", 4000))),
    fitted_mod_original = list(y = 1:20)), class = "multiScaleR")
  key <- profile_scale_cache_key(object)
  expect_lte(nchar(key, type = "bytes"), 100)
  count <- 0
  profile <- function(...) {
    count <<- count + 1
    list(intervals = matrix(c(1, 3), 1, dimnames = list("x", c("2.5%", "97.5%"))),
      method = c(x = "profile"))
  }
  with_mocked_bindings({
    first <- scale_ci_table(object, df = 20, profile = TRUE)
    expect_equal(scale_ci_table(object, df = 20, profile = TRUE), first)
    expect_equal(count, 1)
    object$fitted_mod_original$y[1] <- 9
    scale_ci_table(object, df = 20, profile = TRUE)
    expect_equal(count, 2)
  }, profile_sigma_intervals = profile)
})

test_that("outer Hessian uncertainty requires finite positive curvature", {
  for (h in list(matrix(0, 2, 2), diag(c(-1, 1)), diag(c(NaN, 1)),
      diag(c(Inf, 1)), matrix(c(1, 0, .5, 1), 2))) {
    result <- .msr_refit_hessian_se(h, c(TRUE, TRUE), c(10, 1))
    expect_true(all(is.na(result$SE)))
    expect_true(result$diagnostic$triggered)
  }
  result <- .msr_refit_hessian_se(diag(c(4, 9)), c(TRUE, TRUE), c(10, 1))
  expect_equal(result$SE, c(5, 1/3))
  expect_false(result$diagnostic$triggered)
  result <- .msr_refit_hessian_se(diag(c(4, 0)), c(TRUE, FALSE), c(10, 1))
  expect_equal(result$SE, c(5, NA_real_))
})

test_that("quiet refits retain triggered diagnostics without printing blocks", {
  f <- make_core_fixture()
  fit <- function(verbose) suppressWarnings(multiScale_optim(f$fitted_mod,
    f$kernel_inputs, par = f$opt$optim_results$par, n_cores = 1,
    verbose = verbose))
  with_mocked_bindings({
    quiet <- capture.output(out <- fit(FALSE))
    expect_false(any(grepl("WARNING!!!|Optimization complete", quiet)))
    expect_true(out$diagnostics$max_distance$triggered)
    loud <- capture.output(fit(TRUE))
    expect_true(any(grepl("WARNING!!!", loud, fixed = TRUE)))
  }, .max_distance_diagnostic = function(...) list(code = "max_distance",
    triggered = TRUE, suggested_max_D = 1000))
})

test_that("selection labels distinguish kernel and response with native penalties", {
  dat <- data.frame(y = rep(0:1, 30), x = seq(-1, 1, length.out = 60))
  models <- list(glm(y ~ x, data = dat, family = binomial()),
    glm(y ~ x, data = dat, family = poisson()), lm(y ~ x, dat))
  for (model in models) {
    k <- attr(logLik(model), "df")
    expect_equal(.msr_parameter_count(model), as.integer(k))
    wrapper <- structure(list(opt_mod = model,
      kernel_inputs = list(kernel = "gaussian"),
      scale_est = data.frame(Mean = 2, SE = NA_real_, row.names = "x")),
      class = "multiScaleR")
    tab <- aic_tab(list(wrapper), AICc = FALSE)
    expect_equal(tab$K, k + 1)
    expect_equal(tab$AIC, -2 * as.numeric(logLik(model)) + 2 * (k + 1))
    expect_match(tab$Modnames, paste0("family=", .msr_selection_family_label(model)), fixed = TRUE)
    expect_match(tab$Modnames, "kernel=gaussian", fixed = TRUE)
    expect_equal(bic_tab(list(wrapper))$K, k + 1)
    wrapper$shape_est <- data.frame(Mean = NA_real_, SE = NA_real_, row.names = "x")
    rownames(wrapper$scale_est) <- rownames(wrapper$shape_est) <- "hab"
    wrapper$kernel_inputs$scale_vars <- msr_vars(hab = surface_var("hab", metric = "sq"))
    expect_equal(.msr_selection_parameter_count(wrapper), k + 1)
  }
  dat$duplicate <- dat$x
  alias <- glm(y ~ x + duplicate, data = dat, family = binomial())
  expect_equal(.msr_parameter_count(alias), alias$rank)
  skip_if_not_installed("MASS")
  dat$count <- rep(c(0, 1, 3, 8, 1), 12)
  nb <- suppressWarnings(MASS::glm.nb(count ~ x, dat))
  expect_equal(.msr_parameter_count(nb), nb$rank + 1L)
})

test_that("GAM penalties retain native effective degrees of freedom", {
  skip_if_not_installed("mgcv")
  set.seed(521)
  dat <- data.frame(hab = seq(-2, 2, length.out = 80))
  dat$response <- sin(dat$hab) + rnorm(80, sd = .3)
  fit <- mgcv::gam(response ~ s(hab, k = 5), data = dat, method = "ML")
  k <- attr(logLik(fit), "df")
  expect_equal(.msr_parameter_count(fit), k)
  aic <- aic_tab(list(fit), AICc = FALSE)
  expect_equal(aic$K, k)
  expect_equal(aic$AIC, as.numeric(AIC(fit)))
  expect_equal(bic_tab(list(fit))$BIC, as.numeric(BIC(fit)))
})
