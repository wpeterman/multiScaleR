# These tests check observed-data/criterion bookkeeping, not general MLE quality.
# A single failure-retaining fixture supplies six native fits per test lifecycle.
pcount_comparison_cache <- new.env(parent = emptyenv())
pcount_comparison_fixture <- function() {
  if (exists("failure", pcount_comparison_cache, inherits = FALSE))
    stop(pcount_comparison_cache$failure)
  if (!exists("value", pcount_comparison_cache, inherits = FALSE)) {
    pcount_comparison_cache$value <- tryCatch({
      set.seed(9109)
      raster <- terra::rast(nrows = 25, ncols = 25, xmin = 0, xmax = 500,
        ymin = 0, ymax = 500, crs = "EPSG:26915")
      terra::values(raster) <- rnorm(625)
      names(raster) <- "hab"
      points <- terra::vect(cbind(runif(72, 80, 420), runif(72, 80, 420)), crs = terra::crs(raster))
      prepared <- kernel_prep(points, raster, max_D = 70, sigma = 30, verbose = FALSE)
      dat <- data.frame(Year = factor(rep(1:3, each = 24)), hab = prepared$kernel_dat$hab,
        site = factor(rep(seq_len(24), 3)))
      abundance <- rpois(72, exp(1 + .2 * dat$hab + c(0, .1, -.2)[dat$Year] +
        rnorm(24, sd = .35)[dat$site]))
      y <- matrix(rbinom(216, rep(abundance, 3), .5), nrow = 72)
      a <- matrix(sin(seq_len(216) * .2), nrow = 72)
      b <- matrix(cos(seq_len(216) * .2), nrow = 72)
      a[2L, 1L] <- NA_real_; b[3L, 2L] <- NA_real_
      umf <- unmarked::unmarkedFramePCount(y, siteCovs = dat, obsCovs = list(a = a, b = b))
      cases <- list(
        NB_state = list(formula = ~Year ~Year + hab + (1 | site), mixture = "NB", K = 9L, sd = 1L),
        ZIP_state = list(formula = ~Year ~Year + hab + (1 | site), mixture = "ZIP", K = 9L, sd = 1L),
        P_det = list(formula = ~Year + (1 | site) ~Year + hab, mixture = "P", K = 8L, sd = 1L),
        P_both = list(formula = ~Year + (1 | site) ~Year + hab + (1 | site), mixture = "P", K = 9L, sd = 2L),
        missing_a = list(formula = ~Year + a ~Year + hab + (1 | site), mixture = "P", K = 9L, sd = 1L),
        missing_b = list(formula = ~Year + b ~Year + hab + (1 | site), mixture = "P", K = 9L, sd = 1L))
      fits <- list(); warnings <- list(); snapshots <- list(); calls <- 0L
      for (name in names(cases)) {
        case <- cases[[name]]; messages <- character()
        calls <- calls + 1L
        stopifnot(calls <= 6L)
        fit <- withCallingHandlers(unmarked::pcount(case$formula, data = umf, K = 25,
          mixture = case$mixture, engine = "TMB", se = TRUE, threads = 1), warning = function(w) {
            messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
          })
        fits[[name]] <- fit; warnings[[name]] <- messages
        snapshots[[name]] <- list(opt = fit@opt, estimates = fit@estimates,
          TMB_par = fit@TMB$par, response = fit@data@y, AIC = fit@AIC, negLogLike = fit@negLogLike)
      }
      list(fits = fits, cases = cases, warnings = warnings, snapshots = snapshots, calls = calls)
    }, error = function(e) {
      pcount_comparison_cache$failure <- e
      stop(e)
    })
  }
  pcount_comparison_cache$value
}

# No likelihoods, predictions or criteria are computed for these metadata prototypes.
pcount_metadata_fixture <- function(drop_site = FALSE, rows = seq_len(12L), transformed = FALSE,
                                    change_habitat = FALSE) {
  n <- 12L; y <- matrix(rep(0:2, n), nrow = n, byrow = TRUE)
  site <- data.frame(hab = seq(-1, 1, length.out = n), Year = factor(rep(1:3, each = 4L)),
    site = factor(rep(1:4, 3L)), row.names = paste0("s", seq_len(n)))
  a <- matrix(sin(seq_len(n * 3L) * .2), nrow = n)
  b <- matrix(cos(seq_len(n * 3L) * .2), nrow = n)
  a[2L, 1L] <- NA_real_; b[3L, 2L] <- NA_real_
  if (drop_site) a[4L, ] <- NA_real_
  if (change_habitat) site$hab <- site$hab * 1.1
  frame <- unmarked::unmarkedFramePCount(y[rows, , drop = FALSE],
    siteCovs = site[rows, , drop = FALSE],
    obsCovs = list(a = a[rows, , drop = FALSE], b = b[rows, , drop = FALSE]))
  formula <- if (transformed) ~Year * I(a^2) ~Year + hab + (1 | site) else
    ~Year + a ~Year + hab + (1 | site)
  split <- utils::getFromNamespace("split_formula", "unmarked")(formula)
  names(split) <- c("det", "state")
  design <- utils::getFromNamespace("getDesign", "unmarked")(frame, split, na.rm = FALSE)
  methods::new("unmarkedFitPCount", fitType = "pcount", formula = formula, formlist = split,
    data = frame, sitesRemoved = design$removed.sites, K = 25, mixture = "P",
    negLogLike = NA_real_, AIC = NA_real_)
}

test_that("mixed pcount penalties retain variance and mixture parameters", {
  f <- pcount_comparison_fixture()
  expect_identical(f$calls, 6L)
  for (name in names(f$fits)) {
    model <- f$fits[[name]]; k <- f$cases[[name]]$K
    ll <- as.numeric(unmarked::logLik(model))
    expect_identical(f$warnings[[name]], character())
    expect_true(is.finite(ll))
    expect_equal(.msr_parameter_count(model), k)
    expect_equal(length(model@opt$par), k)
    expect_equal(length(model@TMB$par), k)
    expect_equal((model@AIC + 2 * ll) / 2, as.numeric(k), tolerance = 1e-7)
    expect_equal(length(unmarked::coef(model)), k - f$cases[[name]]$sd)
    expect_gt(length(model@TMB$env$random), k)
    expect_equal(.msr_model_nobs(model), 72L)
    aic <- aic_tab(list(model), AICc = FALSE, mod_names = name)
    aicc <- aic_tab(list(model), mod_names = name)
    bic <- bic_tab(list(model), mod_names = name)
    expect_equal(aic$AIC, -2 * ll + 2 * k, tolerance = 1e-8)
    expect_equal(aicc$AICc, -2 * ll + 2 * k + 2 * k * (k + 1) / (72 - k - 1), tolerance = 1e-8)
    expect_equal(bic$BIC, -2 * ll + log(72) * k, tolerance = 1e-8)
    expect_equal(aic$K, as.numeric(k))
    expect_equal(bic$K, as.numeric(k))
  }
})

test_that("pcount comparisons check actual visits rather than raw frames or site totals", {
  fits <- pcount_comparison_fixture()$fits
  a <- fits$missing_a; b <- fits$missing_b
  expect_identical(a@data, b@data)
  expect_identical(a@data@y, b@data@y)
  expect_identical(a@sitesRemoved, b@sitesRemoved)
  expect_equal(.msr_model_nobs(a), .msr_model_nobs(b))
  native <- utils::getFromNamespace("fl_getY", "unmarked")
  expect_false(identical(native(a), native(b)))
  expect_error(unmarked::fitList(fits = list(a = a, b = b)), "missing covariate values")
  expect_error(aic_tab(list(a, b), AICc = FALSE), "response values")
  expect_error(aic_tab(list(a, b)), "response values")
  expect_error(bic_tab(list(a, b)), "response values")
  expect_error(aic_tab(list(fits$P_det, a), AICc = FALSE), "response values")
  expect_no_error(.msr_validate_comparable_models(fits[1:4]))
  expect_warning(aic <- aic_tab(list(a, a), AICc = FALSE), "redundant")
  expect_warning(bic <- bic_tab(list(a, a)), "redundant")
  expect_equal(aic$K, c(9, 9))
  expect_equal(bic$K, c(9, 9))
})

test_that("native response recovery preserves fitted objects and covariance metadata", {
  f <- pcount_comparison_fixture()
  native <- utils::getFromNamespace("fl_getY", "unmarked")
  for (name in names(f$fits)) {
    model <- f$fits[[name]]; before <- f$snapshots[[name]]
    rows <- setdiff(seq_len(nrow(model@data@y)), model@sitesRemoved)
    expect_equal(.msr_unmarked_pcount_response(model, rows), native(model), tolerance = 0)
    expect_no_error(.msr_selection_signature(model))
    expect_identical(model@opt, before$opt)
    expect_identical(model@estimates, before$estimates)
    expect_identical(model@TMB$par, before$TMB_par)
    expect_identical(model@data@y, before$response)
    expect_identical(model@AIC, before$AIC)
    expect_identical(model@negLogLike, before$negLogLike)
  }
})

test_that("pcount metadata retains row alignment, transforms and omitted-site masks", {
  original <- pcount_metadata_fixture()
  reversed <- pcount_metadata_fixture(rows = 12:1)
  transformed <- pcount_metadata_fixture(transformed = TRUE)
  changed <- pcount_metadata_fixture(change_habitat = TRUE)
  expect_equal(.msr_selection_signature(original), .msr_selection_signature(reversed), tolerance = 0)
  expect_equal(.msr_selection_signature(original), .msr_selection_signature(transformed), tolerance = 0)
  expect_equal(.msr_selection_signature(original), .msr_selection_signature(changed), tolerance = 0)
  expect_no_error(.msr_validate_comparable_models(list(original, reversed)))
  dropped <- pcount_metadata_fixture(drop_site = TRUE)
  subset <- pcount_metadata_fixture(drop_site = TRUE, rows = setdiff(seq_len(12L), 4L))
  expect_equal(.msr_model_nobs(dropped), 11L)
  expect_warning(signature <- .msr_selection_signature(dropped), NA)
  expect_equal(signature, .msr_selection_signature(subset), tolerance = 0)
  expect_no_error(.msr_validate_comparable_models(list(dropped, subset)))
  expect_warning(native <- utils::getFromNamespace("fl_getY", "unmarked")(dropped), "dropped due to missing values")
  expect_equal(.msr_unmarked_pcount_response(dropped, setdiff(seq_len(12L), 4L)), native, tolerance = 0)
})

test_that("inconsistent pcount metadata fails without raw-count fallback", {
  model <- pcount_metadata_fixture()
  corrupted <- model
  corrupted@sitesRemoved <- 4
  expect_error(.msr_selection_signature(corrupted), "retained-site metadata")
  corrupted <- model
  corrupted@formlist$det <- ~unknown_covariate
  expect_error(.msr_selection_signature(corrupted), "Could not recover used pcount responses")
  expect_length(model@sitesRemoved, 0L)
  expect_identical(model@data@y, pcount_metadata_fixture()@data@y)
})
