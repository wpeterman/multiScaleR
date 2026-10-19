make_marginal_year_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(712)
    n_sites <- 30L
    site <- factor(rep(seq_len(n_sites), 3L))
    Year <- factor(rep(1:3, each = n_sites))
    habitat <- rep(stats::rnorm(n_sites), 3L)
    temperature <- stats::runif(length(site), 10, 20)
    effort <- stats::runif(length(site), 0.8, 1.2)
    N <- stats::rpois(length(site),
      exp(1.1 + 0.4 * habitat + 0.2 * (as.integer(Year) - 1)))
    y <- replicate(3L, stats::rbinom(length(N), N, 0.7))
    # Extra identifiers and a different column order exposed the original
    # positional selection bug. Year contrasts are not separate data columns.
    sc <- data.frame(site = site, effort = effort, temperature = temperature,
                     Year = Year, unused = seq_along(site), habitat = habitat)
    umf <- unmarked::unmarkedFramePCount(y = y, siteCovs = sc)
    mod <- unmarked::pcount(~Year ~ Year + habitat + temperature,
                            data = umf, K = 40)
    obj <- structure(list(opt_mod = mod,
      scl_params = list(mean = c(habitat = 0.4), sd = c(habitat = 0.2)),
      join_by = data.frame(site = levels(site))), class = "multiScaleR")
    cache <<- list(obj = obj, mod = mod, sc = sc, umf = umf)
    cache
  }
})

test_that("unmarked factor and unscaled effects use named site-year columns", {
  f <- make_marginal_year_fixture()
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  original <- f$obj$opt_mod@data@siteCovs
  plots <- expect_no_warning(plot_marginal_effects(f$obj, length.out = 7))
  expect_named(plots, c("Year", "habitat", "temperature"))
  expect_identical(f$obj$opt_mod@data@siteCovs, original)
  expect_identical(as.character(plots$Year$data$x), levels(f$sc$Year))

  # Compare both numeric and factor panels with independent native predictions.
  for (v in names(plots)) {
    nd <- f$sc[rep(1L, nrow(plots[[v]]$data)), , drop = FALSE]
    nd$habitat <- mean(f$sc$habitat)
    nd$temperature <- mean(f$sc$temperature)
    if (v == "habitat") {
      nd$habitat <- (plots[[v]]$data$x - 0.4) / 0.2
    } else nd[[v]] <- plots[[v]]$data$x
    native <- unmarked::predict(f$mod, type = "state", newdata = nd, re.form = NA)
    expect_equal(plots[[v]]$data$fit, native$Predicted)
    expect_equal(plots[[v]]$data$lwr, native$lower)
    expect_equal(plots[[v]]$data$upr, native$upper)
    expect_true(all(is.finite(plots[[v]]$data$fit)))
  }
  reordered <- f$obj
  reordered$opt_mod@data@siteCovs <- original[rev(names(original))]
  reordered_plots <- plot_marginal_effects(reordered, length.out = 7)
  expect_equal(reordered_plots$habitat$data, plots$habitat$data)
  expect_equal(reordered_plots$temperature$data, plots$temperature$data)
})

test_that("unmarked transformed terms, offsets and absent intercepts recover base variables", {
  f <- make_marginal_year_fixture()
  mod <- unmarked::pcount(~1 ~ 0 + Year + habitat + I(habitat^2) +
    Year:habitat + offset(log(effort)), data = f$umf, K = 40)
  obj <- f$obj
  obj$opt_mod <- mod
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- expect_no_warning(plot_marginal_effects(obj, length.out = 5))
  expect_named(plots, c("Year", "habitat"))
  nd <- f$sc[rep(1L, 5L), , drop = FALSE]
  nd$effort <- mean(f$sc$effort)
  nd$habitat <- (plots$habitat$data$x - 0.4) / 0.2
  native <- unmarked::predict(mod, type = "state", newdata = nd, re.form = NA)
  expect_equal(plots$habitat$data$fit, native$Predicted)
})

test_that("unmarked random-effect groups are excluded from panels and prediction", {
  f <- make_marginal_year_fixture()
  mod <- unmarked::pcount(~Year ~ Year + habitat + (1 | site),
                          data = f$umf, K = 40, engine = "TMB")
  obj <- f$obj
  obj$opt_mod <- mod
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- expect_no_warning(plot_marginal_effects(obj, length.out = 5))
  expect_named(plots, c("Year", "habitat"))
  nd <- f$sc[rep(1L, 5L), , drop = FALSE]
  nd$habitat <- (plots$habitat$data$x - 0.4) / 0.2
  native <- unmarked::predict(mod, type = "state", newdata = nd, re.form = NA)
  expect_equal(plots$habitat$data$fit, native$Predicted)
  expect_equal(plots$habitat$data$lwr, native$lower)
})

test_that("component variable recovery excludes only offset-only and grouping variables", {
  expr <- ~habitat + temperature + offset(log(effort)) + (habitat | site)
  expect_equal(.msr_marginal_variables(expr, offsets = TRUE),
               c("habitat", "temperature", "effort"))
  expect_equal(.msr_marginal_variables(expr, offsets = FALSE),
               c("habitat", "temperature"))
  expect_equal(.msr_marginal_variables(~effort + offset(log(effort)), FALSE),
               "effort")
})

test_that("intercept-only unmarked states return no covariate panels", {
  f <- make_marginal_year_fixture()
  obj <- f$obj
  obj$opt_mod <- unmarked::pcount(~Year ~ 1, data = f$umf, K = 40)
  expect_length(plot_marginal_effects(obj), 0L)
})

test_that("occupancy models use the state formula even without a detection intercept", {
  f <- make_marginal_year_fixture()
  set.seed(907)
  present <- stats::rbinom(nrow(f$sc), 1, 0.65)
  y <- replicate(3L, present * stats::rbinom(length(present), 1, 0.6))
  umf <- unmarked::unmarkedFrameOccu(y = y, siteCovs = f$sc)
  obj <- f$obj
  obj$opt_mod <- unmarked::occu(~0 + Year ~ Year + habitat, data = umf)
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- expect_no_warning(plot_marginal_effects(obj, length.out = 5))
  expect_named(plots, c("Year", "habitat"))
  nd <- f$sc[rep(1L, 5L), , drop = FALSE]
  nd$habitat <- (plots$habitat$data$x - 0.4) / 0.2
  native <- unmarked::predict(obj$opt_mod, type = "state", newdata = nd)
  expect_equal(plots$habitat$data$fit, native$Predicted)
})

test_that("generic numeric effects are independent of extra columns and response position", {
  f <- make_marginal_year_fixture()
  dat <- f$sc
  dat$id <- paste0("row", seq_len(nrow(dat)))
  dat$y <- rowSums(f$umf@y)
  dat <- dat[c("id", "temperature", "y", "habitat")]
  mod <- stats::glm(y ~ habitat + temperature, family = stats::poisson(), data = dat)
  obj <- f$obj
  obj$opt_mod <- mod
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- expect_no_warning(plot_marginal_effects(obj, length.out = 5))
  nd <- data.frame(habitat = mean(dat$habitat), temperature = plots$temperature$data$x)
  expect_equal(plots$temperature$data$fit,
    as.numeric(stats::predict(mod, newdata = nd, type = "response")))
})

test_that("constant nuisance factors retain explicit fitted contrasts", {
  f <- make_marginal_year_fixture()
  sc <- f$sc
  contrasts(sc$Year) <- stats::contr.sum(3L)
  umf <- unmarked::unmarkedFramePCount(y = f$umf@y, siteCovs = sc)
  mod <- unmarked::pcount(~Year ~ Year + habitat + temperature, data = umf, K = 40)
  obj <- f$obj
  obj$opt_mod <- mod
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- expect_no_warning(plot_marginal_effects(obj, length.out = 5))
  nd <- sc[rep(1L, 5L), , drop = FALSE]
  nd$habitat <- (plots$habitat$data$x - 0.4) / 0.2
  nd$temperature <- mean(sc$temperature)
  X <- stats::model.matrix(~Year + habitat + temperature, data = nd)
  expected <- exp(as.vector(X %*% unmarked::coef(mod, type = "state")))
  expect_equal(plots$habitat$data$fit, expected)
  se <- unname(sqrt(diag(X %*% unmarked::vcov(mod, type = "state") %*% t(X))))
  expect_equal(plots$habitat$data$lwr, exp(log(expected) - stats::qnorm(0.975) * se))
  expect_equal(plots$habitat$data$upr, exp(log(expected) + stats::qnorm(0.975) * se))
})

test_that("custom contrasts retain the occupancy link and zero random effects", {
  f <- make_marginal_year_fixture()
  sc <- f$sc
  contrasts(sc$Year) <- stats::contr.sum(3L)
  set.seed(907)
  present <- stats::rbinom(nrow(sc), 1, 0.65)
  y <- replicate(3L, present * stats::rbinom(length(present), 1, 0.6))
  occ <- unmarked::occu(~Year ~ Year + habitat,
    data = unmarked::unmarkedFrameOccu(y = y, siteCovs = sc))
  ran <- unmarked::pcount(~Year ~ Year + habitat + (1 | site),
    data = unmarked::unmarkedFramePCount(y = f$umf@y, siteCovs = sc),
    K = 40, engine = "TMB")
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  for (mod in list(occ, ran)) {
    obj <- f$obj
    obj$opt_mod <- mod
    plots <- expect_no_warning(plot_marginal_effects(obj, length.out = 5))
    nd <- sc[rep(1L, 5L), , drop = FALSE]
    nd$habitat <- (plots$habitat$data$x - 0.4) / 0.2
    X <- stats::model.matrix(~Year + habitat, data = nd)
    eta <- as.vector(X %*% unmarked::coef(mod, type = "state"))
    inverse <- if (methods::is(mod, "unmarkedFitOccu")) stats::plogis else exp
    expect_equal(plots$habitat$data$fit, inverse(eta))
  }
  # ZIP state prediction is derived from more than the state regression. The
  # custom-contrast fallback must not silently return conditional lambda alone.
  zip <- f$mod
  zip@mixture <- "ZIP"
  expect_error(.msr_unmarked_contrast_prediction(zip, ~habitat,
    data.frame(habitat = 0), "state"), "class-specific prediction")
})
