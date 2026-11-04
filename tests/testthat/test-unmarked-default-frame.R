# Two native fits are cached once; failures remain visible without refitting.
default_unmarked_cache <- new.env(parent = emptyenv())
default_unmarked_fixture <- function() {
  if (exists("failure", default_unmarked_cache, inherits = FALSE))
    stop(default_unmarked_cache$failure)
  if (!exists("value", default_unmarked_cache, inherits = FALSE)) {
    default_unmarked_cache$value <- tryCatch({
      set.seed(9111)
      observer <- outer(seq_len(72L) - 1L, c(0L, 7L, 14L), "+") %% 24L + 1L
      effect <- rnorm(24L, sd = .6)
      abundance <- rpois(72L, exp(1.4))
      y <- matrix(rbinom(216L, rep(abundance, 3L), plogis(-.4 + effect[observer])), nrow = 72L)
      rownames(y) <- paste0("s", seq_len(72L))
      frame <- unmarked::unmarkedFramePCount(y,
        obsCovs = data.frame(observer = factor(as.vector(observer))))
      stopifnot(is.null(frame@siteCovs))
      cases <- list(null = list(formula = ~1 ~1, engine = "C", K = 2L),
        detection = list(formula = ~1 + (1 | observer) ~1, engine = "TMB", K = 3L))
      fits <- list(); warnings <- list(); snapshots <- list(); calls <- 0L
      for (name in names(cases)) {
        case <- cases[[name]]; messages <- character(); calls <- calls + 1L
        stopifnot(calls <= 2L)
        fit <- withCallingHandlers(unmarked::pcount(case$formula, frame,
          K = 25, mixture = "P", engine = case$engine, se = TRUE, threads = 1),
          warning = function(w) {
            messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
          })
        fits[[name]] <- fit; warnings[[name]] <- messages
        snapshots[[name]] <- list(data = fit@data, opt = fit@opt, estimates = fit@estimates,
          AIC = fit@AIC, negLogLike = fit@negLogLike)
      }
      list(fits = fits, cases = cases, warnings = warnings, snapshots = snapshots, calls = calls)
    }, error = function(e) { default_unmarked_cache$failure <- e; stop(e) })
  }
  default_unmarked_cache$value
}

# These are metadata-only prototypes: no fit, criterion or prediction is made.
default_unmarked_metadata <- function(rows = seq_len(12L), explicit_ids = TRUE,
                                      missing = "a", drop_site = FALSE, occupancy = FALSE) {
  y <- matrix(rep(0:2, 12L), nrow = 12L, byrow = TRUE)
  if (occupancy) y <- 1L * (y > 0)
  if (explicit_ids) rownames(y) <- paste0("s", seq_len(12L))
  a <- matrix(sin(seq_len(36L) * .2), nrow = 12L)
  b <- matrix(cos(seq_len(36L) * .2), nrow = 12L)
  a[2L, 1L] <- NA_real_; b[3L, 2L] <- NA_real_
  if (drop_site) a[4L, ] <- NA_real_
  frame <- if (occupancy) unmarked::unmarkedFrameOccu(y[rows, , drop = FALSE],
    obsCovs = list(a = a[rows, , drop = FALSE], b = b[rows, , drop = FALSE])) else
    unmarked::unmarkedFramePCount(y[rows, , drop = FALSE],
      obsCovs = list(a = a[rows, , drop = FALSE], b = b[rows, , drop = FALSE]))
  formula <- if (missing == "a") ~a ~1 else ~b ~1
  split <- utils::getFromNamespace("split_formula", "unmarked")(formula)
  names(split) <- c("det", "state")
  design <- utils::getFromNamespace("getDesign", "unmarked")(frame, split, na.rm = FALSE)
  args <- list(Class = if (occupancy) "unmarkedFitOccu" else "unmarkedFitPCount",
    fitType = if (occupancy) "occu" else "pcount", formula = formula, formlist = split,
    data = frame, sitesRemoved = design$removed.sites, negLogLike = NA_real_, AIC = NA_real_)
  if (!occupancy) args <- c(args, list(K = 25, mixture = "P"))
  do.call(methods::new, args)
}

test_that("unmarked null and detection-only fits retain native sample sizes and penalties", {
  f <- default_unmarked_fixture()
  expect_identical(f$calls, 2L)
  for (name in names(f$fits)) {
    model <- f$fits[[name]]; k <- f$cases[[name]]$K
    ll <- as.numeric(unmarked::logLik(model))
    expect_null(model@data@siteCovs)
    expect_identical(f$warnings[[name]], character())
    expect_identical(.msr_model_nobs(model), 72L)
    expect_equal(.msr_model_nobs(model), unmarked::sampleSize(model))
    expect_identical(.msr_observation_ids(model), paste0("s", seq_len(72L)))
    expect_equal(.msr_parameter_count(model), k)
    aic <- aic_tab(list(model), AICc = FALSE, mod_names = name)
    aicc <- aic_tab(list(model), mod_names = name)
    bic <- bic_tab(list(model), mod_names = name)
    expect_equal(aic$AIC, -2 * ll + 2 * k, tolerance = 1e-8)
    expect_equal(aicc$AICc, -2 * ll + 2 * k + 2 * k * (k + 1) / (72 - k - 1), tolerance = 1e-8)
    expect_equal(bic$BIC, -2 * ll + log(72) * k, tolerance = 1e-8)
    expect_equal(aic$K, as.numeric(k))
    expect_equal(bic$K, as.numeric(k))
    expect_identical(model@data, f$snapshots[[name]]$data)
    expect_identical(model@opt, f$snapshots[[name]]$opt)
    expect_identical(model@estimates, f$snapshots[[name]]$estimates)
    expect_identical(model@AIC, f$snapshots[[name]]$AIC)
    expect_identical(model@negLogLike, f$snapshots[[name]]$negLogLike)
  }
  expect_no_error(.msr_validate_comparable_models(f$fits))
  expect_no_error(aic_tab(f$fits, AICc = FALSE))
  expect_no_error(bic_tab(f$fits))
})

test_that("default unmarked frames preserve original IDs through omissions and reordering", {
  for (occupancy in c(FALSE, TRUE)) {
    original <- default_unmarked_metadata(occupancy = occupancy)
    reversed <- default_unmarked_metadata(rows = 12:1, occupancy = occupancy)
    dropped <- default_unmarked_metadata(drop_site = TRUE, occupancy = occupancy)
    subset <- default_unmarked_metadata(rows = setdiff(seq_len(12L), 4L),
      drop_site = TRUE, occupancy = occupancy)
    expect_identical(.msr_model_nobs(original), 12L)
    expect_identical(.msr_observation_ids(original), paste0("s", seq_len(12L)))
    expect_identical(.msr_observation_ids(reversed), paste0("s", 12:1))
    expect_identical(.msr_observation_ids(dropped), paste0("s", setdiff(seq_len(12L), 4L)))
    expect_identical(.msr_model_nobs(dropped), 11L)
    expect_no_error(.msr_validate_comparable_models(list(original, reversed)))
    expect_no_error(.msr_validate_comparable_models(list(dropped, subset)))
    expect_equal(.msr_selection_signature(original), .msr_selection_signature(reversed), tolerance = 0)
    implicit <- default_unmarked_metadata(explicit_ids = FALSE, drop_site = TRUE, occupancy = occupancy)
    expect_identical(.msr_observation_ids(implicit), as.character(setdiff(seq_len(12L), 4L)))
    expect_identical(.msr_model_nobs(implicit), 11L)
  }
})

test_that("missing site covariates do not bypass used-visit checks", {
  a <- default_unmarked_metadata()
  b <- default_unmarked_metadata(missing = "b")
  expect_false(is.null(.msr_selection_signature(a)))
  expect_false(identical(.msr_selection_signature(a), .msr_selection_signature(b)))
  expect_error(.msr_validate_comparable_models(list(a, b)), "response values")
  for (omitted in list(c(4, 4), NA_real_, -1, 13, 1.5)) {
    corrupted <- a; corrupted@sitesRemoved <- omitted
    expect_error(.msr_model_nobs(corrupted), "omitted-site")
    expect_error(.msr_observation_ids(corrupted), "omitted-site")
    expect_error(.msr_selection_signature(corrupted), "omitted-site")
  }
  expect_null(a@data@siteCovs)
  expect_length(a@sitesRemoved, 0L)
})
