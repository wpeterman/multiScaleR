# Regression tests for the 2026-11-24 integration audit (refit backend).
# Every test starts from a real fit through kernel_prep() and
# multiScale_optim(), because the audited failures only appeared on real fits:
# hand-built objects carried data that real refitted models had lost.
skip_on_cran()

# Survival specials are found by name in the formula environment, as after
# library(survival) in a user script.
if (requireNamespace("survival", quietly = TRUE)) {
  strata <- survival::strata; cluster <- survival::cluster; coxph <- survival::coxph
  frailty <- survival::frailty; Surv <- survival::Surv
}

# One small landscape and Poisson data set shared by the glm-family tests.
audit_glm_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    rs <- terra::subset(sim_rast(dim = 60, resolution = 10, user_seed = 555), c(1, 4))
    terra::crs(rs) <- "EPSG:32617"
    sim <- sim_dat(alpha = 0.25, beta = c(0.75, -0.5), kernel = "gaussian",
                   sigma = c(60, 100), type = "count", n_points = 60,
                   raster_stack = rs, max_D = 250, user_seed = 555)
    pts <- sim$pts
    sf::st_crs(pts) <- 32617
    inputs <- kernel_prep(pts, rs, max_D = 250, verbose = FALSE)
    dat <- data.frame(y = sim$obs, inputs$kernel_dat)
    cache <<- list(rs = rs, inputs = inputs, dat = dat)
    cache
  }
})

draw_all <- function(plots) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (p in plots) print(p)
  invisible(plots)
}

test_that("a model fitted in a function with a local formula can be optimized", {
  f <- audit_glm_fixture()
  fit_local <- function() {
    d <- f$dat
    form <- y ~ bin1 + cont2
    multiScale_optim(glm(form, family = poisson, data = d), f$inputs, verbose = FALSE)
  }
  fit <- fit_local()
  expect_s3_class(fit, "multiScaleR")
  expect_equal(rownames(fit$scale_est), c("bin1", "cont2"))

  # The multi-species loop: the formula is built from the loop variable.
  fits <- lapply("abund", function(sp) {
    d <- f$dat
    d[[sp]] <- d$y
    multiScale_optim(glm(reformulate(c("bin1", "cont2"), sp), family = poisson, data = d),
                     f$inputs, verbose = FALSE)
  })
  expect_s3_class(fits[[1L]], "multiScaleR")
  expect_equal(fits[[1L]]$scale_est$Mean, fit$scale_est$Mean, tolerance = 1e-6)
})

test_that("a workspace object never replaces the model's own local formula", {
  f <- audit_glm_fixture()
  # A different `form` where an unqualified replay would look first.
  assign("form", y ~ bin1, envir = globalenv())
  on.exit(rm("form", envir = globalenv()), add = TRUE)
  fit_local <- function() {
    d <- f$dat
    form <- y ~ bin1 + cont2
    glm(form, family = poisson, data = d)
  }
  fit <- multiScale_optim(fit_local(), f$inputs, verbose = FALSE)
  expect_equal(attr(stats::terms(fit$opt_mod), "term.labels"), c("bin1", "cont2"))
  expect_equal(rownames(fit$scale_est), c("bin1", "cont2"))
})

test_that("an unreplayable fitting call stops before the search with a hint", {
  f <- audit_glm_fixture()
  make_model <- function() {
    fam <- poisson()
    glm(y ~ bin1 + cont2, family = fam, data = f$dat)
  }
  model <- make_model()
  # Remove the only binding of `fam`: no environment can replay this call.
  rm("fam", envir = environment(stats::formula(model)))
  expect_error(
    multiScale_optim(model, f$inputs, verbose = FALSE),
    "object 'fam' not found.*cannot find when it refits"
  )
})

test_that("a covariate shared by the count and zero parts of zeroinfl gets one scale", {
  skip_if_not_installed("pscl")
  f <- audit_glm_fixture()
  d <- f$dat
  d$y[seq(1, nrow(d), by = 4)] <- 0L
  fit <- multiScale_optim(pscl::zeroinfl(y ~ bin1 + cont2 | bin1, data = d),
                          f$inputs, verbose = FALSE)
  expect_s3_class(fit, "multiScaleR")
  expect_equal(rownames(fit$scale_est), c("bin1", "cont2"))
  expect_true(all(is.finite(fit$scale_est$Mean)))
})

test_that("refitted gls models keep data for effect curves, clamping and site layers", {
  skip_if_not_installed("nlme")
  f <- audit_glm_fixture()
  d <- f$dat
  d$yg <- as.numeric(scale(d$y))
  d$site <- seq_len(nrow(d)) / nrow(d) - 0.5
  fit <- multiScale_optim(nlme::gls(yg ~ bin1 + cont2 + site, data = d, method = "ML"),
                          f$inputs, verbose = FALSE)
  expect_true(is.data.frame(.model_data(fit$opt_mod)))
  expect_true(all(c("yg", "bin1", "cont2", "site") %in% names(.model_data(fit$opt_mod))))
  plots <- expect_no_error(plot_marginal_effects(fit))
  expect_named(plots, c("bin1", "cont2", "site"))
  draw_all(plots)
  clamped <- expect_no_warning(kernel_scale.raster(f$rs, multiScaleR = fit,
    scale_center = TRUE, clamp = TRUE, verbose = FALSE))
  expect_true("site" %in% names(clamped))
})

test_that("an offset given as an argument is held at its mean, not recycled", {
  f <- audit_glm_fixture()
  d <- f$dat
  d$eff <- log(seq(1, 3, length.out = nrow(d)))
  fit <- multiScale_optim(glm(y ~ bin1 + cont2, family = poisson, data = d, offset = eff),
                          f$inputs, verbose = FALSE)
  plots <- expect_no_warning(plot_marginal_effects(fit))
  curve <- plots$cont2$data
  nd <- data.frame(bin1 = mean(.model_data(fit$opt_mod)$bin1),
                   cont2 = (curve$x - fit$scl_params$mean[["cont2"]]) / fit$scl_params$sd[["cont2"]],
                   eff = mean(d$eff))
  expect_equal(curve$fit, unname(predict(fit$opt_mod, newdata = nd, type = "response")),
               tolerance = 1e-8)

  # Maps get an offset layer, so prediction never reuses the training offsets.
  r <- kernel_scale.raster(f$rs, multiScaleR = fit, scale_center = TRUE,
                           match_buffer = TRUE, verbose = FALSE)
  expect_true("eff" %in% names(r))
  expect_no_warning(predict(fit, newdata = r))
})

test_that("an offset argument that is not a data column gives a clear curve error", {
  f <- audit_glm_fixture()
  exposure <- log(seq(1, 3, length.out = nrow(f$dat)))
  model <- glm(y ~ bin1 + cont2, family = poisson, data = f$dat, offset = exposure)
  fit <- multiScale_optim(model, f$inputs, verbose = FALSE)
  expect_error(plot_marginal_effects(fit), "offset =.*exposure.*not a column")
})

test_that("clamped rasters work when a covariate enters only through a function", {
  f <- audit_glm_fixture()
  for (form in list(y ~ poly(cont2, 2) + bin1, y ~ splines::ns(cont2, df = 2) + bin1,
                    y ~ log(cont2 + 10) + bin1)) {
    fit <- multiScale_optim(glm(form, family = poisson, data = f$dat), f$inputs,
                            verbose = FALSE)
    r <- expect_no_error(kernel_scale.raster(f$rs, multiScaleR = fit,
      scale_center = TRUE, clamp = TRUE, verbose = FALSE))
    # Clamping limits are the scaled covariate's fitted range.
    observed <- range(.model_data(fit$opt_mod)$cont2)
    mapped <- terra::global(r[["cont2"]], range, na.rm = TRUE)
    expect_gte(min(unlist(mapped)), observed[1] - 1e-8)
    expect_lte(max(unlist(mapped)), observed[2] + 1e-8)
  }
})

test_that("a fixed-radius metric inside poly() supports include_fixed maps", {
  rs <- sim_rast(dim = 60, resolution = 10, user_seed = 555)
  terra::crs(rs) <- "EPSG:32617"
  lc <- terra::classify(rs[["cont1"]], c(-Inf, -1, 0, 1, Inf))
  names(lc) <- "lc"
  layers <- c(rs[["cont2"]], lc)
  sim <- sim_dat(alpha = 0.25, beta = 0.75, kernel = "gaussian", sigma = 80,
                 type = "count", n_points = 60, raster_stack = rs[["cont2"]],
                 max_D = 250, user_seed = 555)
  pts <- sim$pts
  sf::st_crs(pts) <- 32617
  specs <- msr_vars(cont2 = kernel_var("cont2"),
                    ed = landscape_var("lc", "ed", radius = 100))
  inputs <- kernel_prep(pts, layers, max_D = 250, scale_vars = specs, verbose = FALSE)
  d <- data.frame(y = sim$obs, inputs$kernel_dat)
  fit <- multiScale_optim(glm(y ~ cont2 + poly(ed, 2), family = poisson, data = d),
                          inputs, verbose = FALSE)
  r <- expect_no_error(kernel_scale.raster(layers, multiScaleR = fit, scale_center = TRUE,
    match_buffer = TRUE, include_fixed = TRUE, verbose = FALSE))
  expect_true(all(c("cont2", "ed") %in% names(r)))
  expect_s4_class(predict(fit, newdata = r), "SpatRaster")
})

# Step-selection data: 50 choice sets of 6 steps on a projected landscape.
audit_ssf_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    rs <- terra::subset(sim_rast(dim = 80, resolution = 10, user_seed = 555), c(1, 4))
    terra::crs(rs) <- "EPSG:32617"
    names(rs) <- c("land1", "land2")
    set.seed(11)
    n_sets <- 50L; n_alt <- 6L
    sx <- rep(stats::runif(n_sets, 250, 550), each = n_alt)
    sy <- rep(stats::runif(n_sets, 250, 550), each = n_alt)
    len <- stats::rgamma(n_sets * n_alt, 2, scale = 30)
    ang <- stats::runif(n_sets * n_alt, -pi, pi)
    d <- data.frame(step_id = rep(seq_len(n_sets), each = n_alt),
                    x = sx + len * cos(ang), y = sy + len * sin(ang),
                    sl_ = len, ta_ = ang)
    v <- terra::extract(rs[["land1"]], as.matrix(d[, c("x", "y")]))[, 1]
    score <- exp(1.2 * as.numeric(scale(v)) - 0.004 * d$sl_)
    d$case_ <- 0
    for (s in seq_len(n_sets)) {
      rows <- which(d$step_id == s)
      d$case_[sample(rows, 1, prob = score[rows])] <- 1
    }
    d$step_id_ <- d$step_id
    d$animal <- rep(c("a", "b"), each = n_sets * n_alt / 2)
    pts <- terra::vect(d, geom = c("x", "y"), crs = terra::crs(rs))
    inputs <- kernel_prep(pts, rs, max_D = 250, verbose = FALSE)
    cache <<- list(rs = rs, inputs = inputs, dat = data.frame(d, inputs$kernel_dat))
    cache
  }
})

test_that("cos(ta_) and log(sl_) terms support effect curves and maps on both backends", {
  skip_if_not_installed("survival")
  s <- audit_ssf_fixture()
  dat <- s$dat
  model <- survival::clogit(case_ ~ land1 + land2 + log(sl_) + cos(ta_) + strata(step_id),
                            data = dat, model = TRUE)
  for (backend in c("refit", "conditional_logit")) {
    fit <- multiScale_optim(model, s$inputs, verbose = FALSE, backend = backend)
    plots <- expect_no_error(plot_marginal_effects(fit))
    expect_true(all(c("land1", "land2", "sl_", "ta_") %in% names(plots)))
    draw_all(plots)
    # sl_ defaults to 0, where log(sl_) is undefined in every cell.
    expect_warning(
      r <- kernel_scale.raster(s$rs, multiScaleR = fit, scale_center = TRUE, verbose = FALSE),
      "where log\\(sl_\\) is not finite"
    )
    expect_true(all(c("sl_", "ta_") %in% names(r)))
    r2 <- kernel_scale.raster(s$rs, multiScaleR = fit, scale_center = TRUE,
                              match_buffer = TRUE, verbose = FALSE,
                              site_values = list(sl_ = mean(dat$sl_)))
    risk <- predict(fit, newdata = r2, type = "risk")
    expect_true(all(is.finite(range(terra::values(risk), na.rm = TRUE))))
  }
})

test_that("Cox fits without model = TRUE predict and plot after a refit", {
  skip_if_not_installed("survival")
  s <- audit_ssf_fixture()
  dat <- s$dat
  model <- survival::clogit(case_ ~ land1 + land2 + strata(step_id), data = dat)
  fit <- multiScale_optim(model, s$inputs, verbose = FALSE)
  expect_true(is.numeric(predict(fit)))
  draw_all(expect_no_error(plot_marginal_effects(fit)))
})

test_that("cluster() terms refit, and frailty() curves refuse clearly", {
  skip_if_not_installed("survival")
  s <- audit_ssf_fixture()
  dat <- s$dat
  model <- survival::coxph(Surv(rep(1, nrow(dat)), case_) ~ land1 + land2 +
                             strata(step_id) + cluster(animal), data = dat)
  expect_s3_class(multiScale_optim(model, s$inputs, verbose = FALSE), "multiScaleR")

  frail <- suppressWarnings(survival::coxph(Surv(rep(1, nrow(dat)), case_) ~
    land1 + land2 + frailty(animal), data = dat))
  fit <- suppressWarnings(multiScale_optim(frail, s$inputs, verbose = FALSE))
  expect_error(plot_marginal_effects(fit), "penalized Cox models")
})

test_that("amt fits name the missing data instead of failing after the search", {
  skip_if_not_installed("amt")
  skip_if_not_installed("survival")
  s <- audit_ssf_fixture()
  dat <- s$dat
  no_model <- amt::fit_clogit(dat, case_ ~ land1 + land2 + strata(step_id_))
  expect_error(multiScale_optim(no_model, s$inputs, verbose = FALSE),
               "model = TRUE")
  inline <- amt::fit_issf(dat, case_ ~ land1 + land2 + sl_ + log(sl_) + cos(ta_) +
                            strata(step_id_), model = TRUE)
  expect_error(multiScale_optim(inline, s$inputs, verbose = FALSE),
               "raw values of .*ta_.*cos\\(ta_\\) as data columns")
  # The documented alternative works: precomputed terms in the data.
  dat$cos_ta_ <- cos(dat$ta_)
  dat$log_sl_ <- log(dat$sl_)
  pre <- amt::fit_issf(dat, case_ ~ land1 + land2 + sl_ + log_sl_ + cos_ta_ +
                         strata(step_id_), model = TRUE)
  expect_s3_class(multiScale_optim(pre, s$inputs, verbose = FALSE), "multiScaleR")
})

test_that("model data recovery never returns a non-data-frame object", {
  skip_if_not_installed("survival")
  s <- audit_ssf_fixture()
  model <- local({
    data <- s$dat
    survival::clogit(case_ ~ land1 + strata(step_id), data = data)
  })
  # The call names `data`, which outside that frame is utils::data().
  out <- .model_data(model)
  expect_true(is.null(out) || is.data.frame(out))
})

test_that("summary() prints a profile_glm NB2 fit at the Poisson limit", {
  skip_if_not_installed("MASS")
  f <- profile_glm_fixture()
  # Underdispersed counts put theta on its Poisson limit, where MASS's own
  # summary printer fails on the NA theta standard error.
  f$data$y <- round(f$mu)
  f$model <- glm(y ~ hab + z + group + offset(log(exposure)), data = f$data,
                 family = MASS::negative.binomial(2))
  x <- suppressWarnings(profile_glm_fit(f))
  printed <- expect_no_error(utils::capture.output(print(summary(x))))
  expect_true(any(grepl("hab", printed)))
})

test_that("unmarked effect curves accept each class's state name and explain empty components", {
  skip_if_not_installed("unmarked")
  rs <- terra::subset(sim_rast(dim = 60, resolution = 10, user_seed = 555), c(1, 4))
  terra::crs(rs) <- "EPSG:32617"
  u <- sim_dat_unmarked(alpha = 0, beta = c(1, -0.5), kernel = "gaussian",
                        sigma = c(60, 100), type = "occ", n_points = 50, n_surv = 4,
                        det = 0.5, raster_stack = rs, max_D = 250, user_seed = 3)
  pts <- u$pts
  sf::st_crs(pts) <- 32617
  inputs <- kernel_prep(pts, rs, max_D = 250, verbose = FALSE)
  umf <- unmarked::unmarkedFrameOccu(y = u$y, siteCovs = inputs$kernel_dat)

  det_only <- multiScale_optim(unmarked::occu(~ cont2 ~ 1, umf), inputs, verbose = FALSE)
  expect_warning(out <- plot_marginal_effects(det_only),
                 "'state' component.*Components with covariates: det")
  expect_length(out, 0L)
  expect_length(plot_marginal_effects(det_only, type = "det"), 1L)

  multi <- unmarked::unmarkedMultFrame(y = cbind(u$y, u$y), siteCovs = inputs$kernel_dat,
                                       numPrimary = 2)
  col <- multiScale_optim(unmarked::colext(~ bin1 + cont2, ~1, ~1, ~1, multi), inputs,
                          verbose = FALSE)
  for (type in c("state", "psi", "lambda"))
    expect_named(draw_all(plot_marginal_effects(col, type = type)), c("bin1", "cont2"))
  expect_error(plot_marginal_effects(col, type = "col"), "site-level state component")

  rn <- multiScale_optim(unmarked::occuRN(~1 ~ bin1 + cont2, umf, K = 15), inputs,
                         verbose = FALSE)
  expect_named(plot_marginal_effects(rn), c("bin1", "cont2"))
  expect_named(plot_marginal_effects(rn, type = "lambda"), c("bin1", "cont2"))
})
