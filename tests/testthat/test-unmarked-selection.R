selection_unmarked_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(9109)
    raster <- terra::rast(nrows = 25, ncols = 25, xmin = 0, xmax = 500,
                          ymin = 0, ymax = 500, crs = "EPSG:26915")
    terra::values(raster) <- rnorm(625)
    names(raster) <- "hab"
    points <- terra::vect(cbind(runif(72, 80, 420), runif(72, 80, 420)),
                           crs = terra::crs(raster))
    prepared <- kernel_prep(points, raster, max_D = 70, sigma = 30,
                             verbose = FALSE)
    dat <- data.frame(Year = factor(rep(1:3, each = 24)),
                      hab = prepared$kernel_dat$hab,
                      site = factor(rep(seq_len(24), 3)))
    abundance <- rpois(72, exp(1 + 0.2 * dat$hab +
                         c(0, 0.1, -0.2)[dat$Year] +
                         rnorm(24, sd = 0.35)[dat$site]))
    y <- matrix(rbinom(216, rep(abundance, 3), 0.5), nrow = 72)
    umf <- unmarked::unmarkedFramePCount(y, siteCovs = dat)
    mod <- unmarked::pcount(~Year ~Year + hab + (1 | site), data = umf,
                             K = 25, engine = "TMB")
    plain <- unmarked::pcount(~1 ~hab, data = umf, K = 25, se = FALSE)
    fit <- structure(list(
      opt_mod = mod, kernel_inputs = prepared,
      scale_est = data.frame(Mean = 30, SE = NA_real_, row.names = "hab"),
      shape_est = NULL, min_D = prepared$min_D, max_D = prepared$max_D,
      opt_context = build_opt_context(mod, prepared$raw_cov,
        scale_vars = prepared$scale_vars, unit_conv = prepared$unit_conv,
        resolution = prepared$resolution, n_cols = prepared$n_cols,
        binned = prepared$binned), join_by = NULL
    ), class = "multiScaleR")
    cache <<- list(model = mod, plain = plain, fit = fit, umf = umf)
    cache
  }
})

test_that("TMB count includes all components and variance but excludes latent effects", {
  fix <- selection_unmarked_fixture()
  # Known designs: 4 abundance coefficients, 3 detection, 1 variance.
  expect_equal(.msr_parameter_count(fix$model), 8L)
  expect_equal(.msr_selection_parameter_count(fix$fit), 9L)
  expect_equal(.msr_parameter_count(structure(list(model = fix$model),
                                             class = "fit_clogit")), 8L)
  expect_equal(.msr_model_loglik(fix$model), as.numeric(unmarked::logLik(fix$model)))
  expect_equal((fix$model@AIC + 2 * unmarked::logLik(fix$model)) / 2, 8)
  expect_gt(length(fix$model@TMB$env$random), 8)
  prepared <- prepare_prediction(fix$fit, component = "state",
    native_args = list(re.form = NA))
  expect_equal(predict(prepared), unmarked::predict(fix$model, type = "state", re.form = NA))
  expect_equal(predict(unserialize(serialize(prepared, NULL))), predict(prepared))
  expow <- fix$fit
  expow$kernel_inputs$kernel <- "expow"
  expow$shape_est <- data.frame(Mean = 2, SE = NA_real_, row.names = "hab")
  expect_equal(.msr_selection_parameter_count(expow), 10L)
  # Fixed spatial descriptors do not add a row to scale_est or shape_est.
  expow$kernel_inputs$scale_vars <- msr_vars(
    hab = kernel_var("hab"), fixed_hab = surface_var("hab", "sq", radius = 25))
  expect_equal(.msr_selection_parameter_count(expow), 10L)
})

test_that("C and TMB counts handle factors, interactions and mixture parameters", {
  fix <- selection_unmarked_fixture()
  for (engine in c("C", "TMB")) {
    for (mixture in c("P", "NB", "ZIP")) {
      model <- suppressWarnings(unmarked::pcount(
        ~Year ~Year * hab, data = fix$umf, K = 25,
        engine = engine, mixture = mixture, se = FALSE))
      # Six abundance design columns plus three detection columns.
      expected <- 9L + as.integer(mixture != "P")
      expect_equal(.msr_parameter_count(model), expected)
      expect_equal((model@AIC + 2 * unmarked::logLik(model)) / 2,
                   as.numeric(expected), tolerance = 1e-7)
    }
  }
  other_K <- unmarked::pcount(~1 ~hab, data = fix$umf, K = 35, se = FALSE)
  expect_equal(.msr_parameter_count(other_K), 3L)
  expect_equal(.msr_parameter_count(fix$plain), 3L)
})

test_that("component formula classes count every free component and exclude fixed ones", {
  set.seed(9130)
  dat <- data.frame(Year = factor(rep(1:3, each = 10)), hab = rnorm(30))
  counts <- matrix(rpois(180, 2), nrow = 30)
  umf <- unmarked::unmarkedFramePCO(counts, siteCovs = dat, numPrimary = 3)
  open <- unmarked::pcountOpen(~Year, ~hab, ~hab, ~Year, data = umf,
                                K = 15, se = FALSE)
  expect_equal(.msr_parameter_count(open), 10L)
  fixed <- unmarked::pcountOpen(~Year, ~1, ~1, ~Year, data = umf,
                                 K = 15, fix = "gamma", se = FALSE)
  expect_equal(.msr_parameter_count(fixed), 7L)
  occ <- unmarked::unmarkedMultFrame(matrix(rbinom(180, 1, 0.5), nrow = 30),
                                     siteCovs = dat, numPrimary = 3)
  dynamic <- unmarked::colext(~Year, ~hab, ~hab, ~Year,
                               data = occ, se = FALSE)
  expect_equal(.msr_parameter_count(dynamic), 10L)
  expect_warning(aic_tab(list(open, fixed)), NA)
  simpler <- unmarked::colext(~1, ~1, ~1, ~1, data = occ, se = FALSE)
  expect_warning(bic_tab(list(dynamic, simpler)), NA)
})

test_that("selection scores, ranks and weights match independent unmarked oracles", {
  fix <- selection_unmarked_fixture()
  models <- list(random = fix$model, plain = fix$plain, spatial = fix$fit)
  k <- c(8, 3, 9)
  ll <- c(unmarked::logLik(fix$model), unmarked::logLik(fix$plain),
          unmarked::logLik(fix$model))
  n <- 72
  expected <- list(AIC = -2 * ll + 2 * k,
    AICc = -2 * ll + 2 * k + 2 * k * (k + 1) / (n - k - 1),
    BIC = -2 * ll + k * log(n))
  tabs <- list(AIC = aic_tab(models, AICc = FALSE, mod_names = names(models)),
               AICc = aic_tab(models, mod_names = names(models)),
               BIC = bic_tab(models, mod_names = names(models)))
  for (criterion in names(tabs)) {
    tab <- tabs[[criterion]]
    index <- order(expected[[criterion]])
    delta <- expected[[criterion]][index] - min(expected[[criterion]])
    weights <- exp(-delta / 2) / sum(exp(-delta / 2))
    expect_equal(tab$Modnames, names(models)[index])
    expect_equal(tab$K, k[index])
    expect_equal(tab$LL, ll[index])
    expect_equal(tab[[criterion]], expected[[criterion]][index])
    expect_equal(tab[[paste0("Delta_", criterion)]], delta)
    expect_equal(tab[[paste0(criterion, "Wt")]], weights)
    expect_equal(tab$Cum.Wt, cumsum(weights))
  }
  expect_warning(aic_tab(models), NA)
})

test_that("profile LL agrees with refits and AICc uses nine parameters", {
  fix <- selection_unmarked_fixture()
  sigmas <- c(25, 30, 40)
  prof <- profile_sigma(fix$fit, sigma_values = sigmas, metric = "LL",
                         verbose = FALSE)$profiles
  native_ll <- vapply(sigmas, function(s) {
    refit <- kernel_scale_fn(s / 70, fix$fit$kernel_inputs$d_list,
      fix$fit$kernel_inputs$raw_cov, "gaussian", fix$model,
      opt_context = fix$fit$opt_context, mod_return = TRUE)
    as.numeric(unmarked::logLik(refit$mod))
  }, numeric(1))
  expect_equal(prof$LL, native_ll, tolerance = 1e-7)
  expect_equal(prof$AICc, -2 * native_ll + 18 + 180 / (72 - 10),
               tolerance = 1e-7)
  core <- make_core_fixture()
  short <- with_mocked_bindings(
    profile_sigma(core$opt, n_pts = 3, verbose = FALSE),
    .msr_model_nobs = function(model) 1L, .package = "multiScaleR")
  expect_true(all(is.na(short$profiles$AICc)))
})

test_that("unreliable counts and likelihoods fail clearly without predictor guesses", {
  model <- selection_unmarked_fixture()$model
  malformed <- model
  malformed@AIC <- numeric()
  expect_error(.msr_parameter_count(malformed), "native AIC")
  malformed <- model
  malformed@AIC <- model@AIC + 0.4
  expect_error(.msr_parameter_count(malformed), "invalid count implied")
  malformed@AIC <- model@AIC + 2
  expect_error(.msr_parameter_count(malformed), "Inconsistent native AIC")
  malformed <- model
  malformed@opt$par <- "bad"
  expect_error(.msr_parameter_count(malformed), "malformed optimizer")
  malformed <- model
  malformed@TMB$par <- malformed@TMB$par[-1]
  expect_error(.msr_parameter_count(malformed), "TMB free parameter metadata")
  malformed <- model
  malformed@negLogLike <- Inf
  expect_error(.msr_parameter_count(malformed), "finite scalar log-likelihood")
  malformed <- model
  malformed@opt$par <- NULL
  expect_equal(.msr_parameter_count(malformed), 8L)
  boundary <- model
  boundary@opt$par[5] <- -100
  expect_equal(.msr_parameter_count(boundary), 8L)
  malformed <- model
  malformed@AIC <- NA_real_
  expect_error(.msr_parameter_count(malformed), "native AIC")
  malformed@AIC <- -2 * unmarked::logLik(model) - 2
  expect_error(.msr_parameter_count(malformed), "invalid count implied")
  expect_error(.msr_parameter_count(structure(list(), class = "unknownFit")),
                "unknownFit.*no reliable parameter metadata")
  expect_error(with_mocked_bindings(.msr_model_loglik(lm(mpg ~ wt, mtcars)),
                 logLik = function(...) c(-1, -2), .package = "multiScaleR"),
                "finite scalar")
  tiny <- lm(c(1, 3, 2, 5) ~ c(1, 2, 3, 4))
  expect_equal(.msr_parameter_count(tiny), 3L)
  expect_error(aic_tab(list(tiny, tiny)), "n > K")
  intercept <- lm(c(1, 3, 2, 5) ~ 1)
  expect_s3_class(aic_tab(list(tiny, intercept), AICc = FALSE), "aictab")
})

test_that("namespace and attached unmarked dispatch agree", {
  fix <- selection_unmarked_fixture()
  before <- .msr_parameter_count(fix$model)
  was_attached <- "package:unmarked" %in% search()
  if (!was_attached) {
    suppressPackageStartupMessages(library(unmarked))
    on.exit(detach("package:unmarked", unload = FALSE), add = TRUE)
  }
  expect_warning(after <- .msr_parameter_count(fix$model), NA)
  expect_identical(after, before)
  expect_identical(.msr_model_loglik(fix$model), as.numeric(unmarked::logLik(fix$model)))
})

test_that("Gaussian likelihood variance and spatial parameters are counted once", {
  dat <- data.frame(y = c(1, 3, 2, 5, 4, 7), hab = 1:6,
                    fixed_hab = c(0, 1, 0, 0, 1, 0))
  model <- lm(y ~ hab + fixed_hab, data = dat)
  # Three fitted coefficients plus residual variance; the fixed radius is not free.
  expect_equal(.msr_parameter_count(model), 4L)
  wrapper <- structure(list(opt_mod = model,
    kernel_inputs = list(scale_vars = msr_vars(hab = kernel_var("hab"),
      fixed_hab = surface_var("hab", "sq", radius = 25))),
    scale_est = data.frame(Mean = 30, SE = 1, row.names = "hab"),
    shape_est = NULL), class = "multiScaleR")
  expect_equal(.msr_selection_parameter_count(wrapper), 5L)
  wrapper$shape_est <- data.frame(Mean = 2, SE = 1, row.names = "hab")
  expect_equal(.msr_selection_parameter_count(wrapper), 6L)
})
