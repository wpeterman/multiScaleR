# multiScale_joint(): scales (and optionally coefficients) shared across
# independent data sets. Every fit here goes through the public API.

test_that("footprint_sd() converts sampling geometries to per-axis SDs", {
  expect_equal(footprint_sd("point"), 0)
  expect_equal(footprint_sd("disk", 100), 50)
  expect_equal(footprint_sd("square", sqrt(12)), 1)
  expect_equal(footprint_sd("rectangle", 6, 6), sqrt(72 / 24))
  expect_equal(footprint_sd("halfnormal", c(40, 60)), c(40, 60))
  expect_error(footprint_sd("disk", -1), "positive")
  expect_error(footprint_sd("rectangle", 100), "width")
})

test_that("multiScale_joint() validates its inputs before fitting", {
  expect_error(multiScale_joint(list()), "two or more")
  expect_error(multiScale_joint(list(1, 2)), "named list")
  expect_error(multiScale_joint(list(a = 1, b = 2)), "multiScale_optim")
})

skip_on_cran()

joint_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    rs <- sim_rast(dim = 60, resolution = 10, user_seed = 11)
    r <- rs[["cont2"]]
    names(r) <- "hab"
    terra::crs(r) <- "EPSG:32617"
    make <- function(type, n, seed, alpha) {
      sim <- sim_dat(alpha = alpha, beta = 0.8, kernel = "gaussian", sigma = 60,
                     type = type, n_points = n, raster_stack = r, max_D = 250,
                     user_seed = seed)
      pts <- sim$pts
      sf::st_crs(pts) <- 32617
      inputs <- kernel_prep(pts, r, max_D = 250, kernel = "gaussian", verbose = FALSE)
      inputs$min_D <- 10
      list(y = sim$obs, inputs = inputs)
    }
    a <- make("count", 50, 1, 0.3)
    b <- make("occ", 60, 2, 0)
    c <- make("count", 40, 3, 1)
    da <- data.frame(y = a$y, a$inputs$kernel_dat)
    db <- data.frame(y = b$y, b$inputs$kernel_dat)
    dc <- data.frame(y = c$y, c$inputs$kernel_dat)
    fa <- multiScale_optim(glm(y ~ hab, family = poisson, data = da), a$inputs, verbose = FALSE)
    fb <- multiScale_optim(glm(y ~ hab, family = binomial(link = "cloglog"), data = db),
                           b$inputs, verbose = FALSE)
    fc <- multiScale_optim(glm(y ~ hab, family = poisson, data = dc), c$inputs, verbose = FALSE)
    cache <<- list(raster = r, counts = fa, presence = fb, counts2 = fc,
                   presence_data = db, presence_inputs = b$inputs)
    cache
  }
})

test_that("a common scale maximizes the summed profile log-likelihood", {
  f <- joint_fixture()
  fits <- list(counts = f$counts, presence = f$presence)
  joint <- suppressWarnings(multiScale_joint(fits, scales = "common", verbose = FALSE))
  expect_s3_class(joint, "multiScaleR_joint")
  grid <- exp(seq(log(15), log(240), length.out = 60))
  ll <- rowSums(sapply(fits, function(x)
    profile_sigma(x, sigma_values = grid, metric = "LL", verbose = FALSE)$profiles$LL))
  expect_gte(joint$logLik, max(ll) - 1e-3)
  sigma <- unname(coef(joint)["sigma:hab"])
  expect_lt(abs(log(sigma) - log(grid[which.max(ll)])), 2 * diff(log(grid))[1])
  # Profile limits bracket the estimate, and both data sets report the same scale.
  p <- joint$parameters
  expect_true(p$profile_lower < p$estimate && p$estimate < p$profile_upper)
  expect_equal(joint$scales$sigma, rep(sigma, 2))
  expect_true(is.finite(joint$parameters$SE))
})

test_that("separate scales reproduce the separate fits", {
  f <- joint_fixture()
  fits <- list(counts = f$counts, presence = f$presence)
  joint <- suppressWarnings(multiScale_joint(fits, scales = "separate", profile = FALSE,
                                             verbose = FALSE))
  separate_ll <- sum(vapply(fits, function(x) as.numeric(logLik(x$opt_mod)), numeric(1)))
  expect_gte(joint$logLik, separate_ll - 1e-6)
  k <- sum(vapply(fits, function(x) attr(logLik(x$opt_mod), "df") + 1, numeric(1)))
  expect_equal(joint$df, k)
  expect_equal(joint$comparison$df_difference, 0)
})

test_that("zero footprints give the common scale and footprints add in variance", {
  f <- joint_fixture()
  fits <- list(counts = f$counts, presence = f$presence)
  common <- suppressWarnings(multiScale_joint(fits, profile = FALSE, verbose = FALSE))
  zero <- suppressWarnings(multiScale_joint(fits, scales = "footprint",
                                            footprint = c(counts = 0, presence = 0),
                                            profile = FALSE, verbose = FALSE))
  expect_equal(zero$logLik, common$logLik, tolerance = 1e-6)
  expect_equal(coef(zero)[["sigma:hab"]], coef(common)[["sigma:hab"]], tolerance = 1e-3)
  tau <- c(counts = 0, presence = footprint_sd("disk", 80))
  fp <- suppressWarnings(multiScale_joint(fits, scales = "footprint", footprint = tau,
                                          verbose = FALSE))
  eco <- coef(fp)[["sigma:hab"]]
  expect_equal(fp$scales$sigma, sqrt(eco^2 + tau[fp$scales$data_set]^2), ignore_attr = TRUE)
  # Profile limits of the apparent scale follow from the ecological limits.
  pres <- fp$scales$data_set == "presence"
  expect_equal(fp$scales$profile_lower[pres],
               sqrt(fp$parameters$profile_lower^2 + tau[["presence"]]^2))
  expect_error(multiScale_joint(fits, scales = "footprint", verbose = FALSE), "footprint")
  expect_error(multiScale_joint(fits, scales = "footprint_scaled",
                                footprint = c(counts = 0, presence = 0), verbose = FALSE),
               "two different footprint")
})

test_that("a shared coefficient equals one stacked model with data-set intercepts", {
  f <- joint_fixture()
  joint <- suppressWarnings(multiScale_joint(list(a = f$counts, c = f$counts2),
                                             scales = "common", share = "hab",
                                             profile = FALSE, verbose = FALSE))
  stacked <- rbind(data.frame(ds = "a", joint$fits$a$opt_mod$data[c("y", "hab")]),
                   data.frame(ds = "c", joint$fits$c$opt_mod$data[c("y", "hab")]))
  m <- glm(y ~ 0 + ds + hab, family = poisson, data = stacked)
  expect_equal(joint$logLik, as.numeric(logLik(m)), tolerance = 1e-6)
  expect_equal(coef(joint)[["beta:hab"]], coef(m)[["hab"]], tolerance = 1e-4)
  expect_equal(stats::sd(stacked$hab), 1)
  # Counts with a log link and presence with a cloglog link share log intensity.
  mixed <- suppressWarnings(multiScale_joint(list(counts = f$counts, presence = f$presence),
                                             share = "hab", profile = FALSE, verbose = FALSE))
  expect_true(is.finite(coef(mixed)[["beta:hab"]]))
})

test_that("shared coefficients refuse incomparable links and nonlinear terms", {
  f <- joint_fixture()
  db <- f$presence_data
  logit <- multiScale_optim(glm(y ~ hab, family = binomial, data = db),
                            f$presence_inputs, verbose = FALSE)
  expect_error(multiScale_joint(list(counts = f$counts, presence = logit), share = "hab",
                                verbose = FALSE), "comparable links")
  expect_error(multiScale_joint(list(counts = f$counts, presence = f$presence),
                                share = "road", verbose = FALSE), "no fit uses")
})

test_that("joint fits support logLik, AIC, coef, vcov, nobs, anova, print, summary", {
  f <- joint_fixture()
  fits <- list(counts = f$counts, presence = f$presence)
  common <- suppressWarnings(multiScale_joint(fits, profile = FALSE, verbose = FALSE))
  separate <- suppressWarnings(multiScale_joint(fits, scales = "separate", profile = FALSE,
                                                verbose = FALSE))
  ll <- logLik(common)
  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "df"), common$df)
  expect_equal(nobs(common), 110)
  expect_equal(AIC(common), -2 * common$logLik + 2 * common$df)
  expect_equal(nrow(AIC(common, separate)), 2L)
  expect_named(coef(common), "sigma:hab")
  expect_equal(dim(vcov(common)), c(1L, 1L))
  a <- anova(common, separate)
  expect_s3_class(a, "anova")
  expect_equal(a$df_difference[2], 1)
  expect_output(print(common), "Joint scale-of-effect fit")
  expect_output(print(summary(common)), "Native coefficients")
})

test_that("component fits work with the existing post-fit tools", {
  f <- joint_fixture()
  joint <- suppressWarnings(multiScale_joint(list(counts = f$counts, presence = f$presence),
                                             profile = FALSE, verbose = FALSE))
  comp <- joint$fits$presence
  expect_s3_class(comp, "multiScaleR")
  expect_equal(comp$scale_est$Mean, unname(coef(joint)["sigma:hab"]))
  expect_equal(comp$joint$data_set, "presence")
  expect_s3_class(summary(comp), "summary_multiScaleR")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(utils::capture.output(print(plot(comp))))
  expect_no_error(utils::capture.output(plot_marginal_effects(comp)))
  sm <- kernel_scale.raster(f$raster, multiScaleR = comp, scale_center = TRUE,
                            match_buffer = TRUE, verbose = FALSE)
  expect_s4_class(predict(comp, newdata = sm), "SpatRaster")
  # A component cannot be combined again as if it were a separate fit.
  expect_error(multiScale_joint(list(a = comp, b = f$counts)), "component of a joint fit")
})

test_that("cluster identifiers give residual correlations and jackknife SEs", {
  f <- joint_fixture()
  fits <- list(counts = f$counts, counts2 = f$counts2)
  ids <- list(counts = paste0("s", seq_len(50)), counts2 = paste0("s", seq_len(40)))
  joint <- suppressWarnings(multiScale_joint(fits, cluster = ids, uncertainty = "jackknife",
                                             groups = 5, seed = 1, profile = FALSE,
                                             verbose = FALSE))
  expect_equal(joint$residual_correlation$shared_clusters, 40)
  expect_true(is.finite(joint$parameters$SE))
  expect_true(joint$uncertainty$jackknife$replicates >= 2)
  expect_error(multiScale_joint(fits, uncertainty = "jackknife", verbose = FALSE), "cluster")
})
