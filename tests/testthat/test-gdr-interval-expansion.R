gdr_expansion_fixture <- function(primary = 1L) {
  site <- data.frame(forest = c(.2, .8), Observer = factor(c("None", "B"),
    levels = c("None", "B", "unused")))
  contrasts(site$Observer) <- stats::contr.sum(3)
  site[["wind m/s"]] <- c(2, 4)
  site$matrix_field <- I(cbind(left = 1:2, right = 3:4))
  yearly <- data.frame(survey_wind = seq_len(2L * primary),
    season = ordered(rep(c("early", "late"), length.out = 2L * primary)))
  frame <- unmarked::unmarkedFrameGDR(matrix(0, 2, 2 * primary),
    matrix(0, 2, 3 * primary), numPrimary = primary, siteCovs = site,
    yearlySiteCovs = yearly, obsCovs = data.frame(recorded = seq_len(6 * primary)),
    dist.breaks = c(0, 25, 75), unitsIn = "m", period.lengths = c(1, 2, 1))
  keys <- data.frame(site = rep(c("a", "b"), each = primary),
    survey = paste0("survey", seq_len(2L * primary)))
  list(frame = frame, keys = keys)
}

test_that("explicit GDR expansion preserves every native field and source attribute", {
  for (primary in 1:2) {
    f <- gdr_expansion_fixture(primary)
    original <- f$frame
    out <- expand_gdr_intervals(f$frame, f$keys, c("site", "survey"), "minutes",
      site_covariates = c(rem_forest = "forest", rem_observer = "Observer",
        "wind m/s", rem_matrix = "matrix_field"),
      yearly_covariates = c(rem_wind = "survey_wind", rem_season = "season"))
    expect_s3_class(out, "msr_gdr_intervals")
    expect_false(out$fit_eligibility_established)
    expect_false(out$protocol$reordered)
    for (slot in setdiff(methods::slotNames(f$frame), "obsCovs"))
      expect_identical(methods::slot(out$data, slot), methods::slot(original, slot))
    actual <- unmarked::obsCovs(out$data)
    site_rows <- rep(1:2, each = primary * 3)
    event_rows <- rep(seq_len(2 * primary), each = 3)
    expect_identical(actual$recorded, unmarked::obsCovs(original)$recorded)
    expect_identical(actual$rem_forest, unmarked::siteCovs(original)$forest[site_rows])
    expect_identical(actual$rem_observer, unmarked::siteCovs(original)$Observer[site_rows])
    expect_identical(contrasts(actual$rem_observer), contrasts(unmarked::siteCovs(original)$Observer))
    expect_identical(actual[["wind m/s"]], unmarked::siteCovs(original)[["wind m/s"]][site_rows])
    expect_identical(actual$rem_matrix, unmarked::siteCovs(original)$matrix_field[site_rows, , drop = FALSE])
    expect_identical(actual$rem_wind, unmarked::yearlySiteCovs(original)$survey_wind[event_rows])
    expect_identical(actual$rem_season, unmarked::yearlySiteCovs(original)$season[event_rows])
    expect_identical(actual$interval, rep(1:3, 2L * primary))
    expect_identical(out$protocol$expanded_covariates$target,
      c("rem_forest", "rem_observer", "wind m/s", "rem_matrix", "rem_wind", "rem_season"))
    expect_identical(out$protocol$expanded_covariates$source_table, c(rep("site", 4), rep("yearly", 2)))
    expect_identical(out$protocol$shared_site_or_annual_names, "wind m/s")
    expect_identical(f$frame, original)
  }
})

test_that("native removal fixed and random designs receive all selected observation fields", {
  f <- gdr_expansion_fixture(2L)
  out <- expand_gdr_intervals(f$frame, f$keys, c("site", "survey"), "minutes",
    site_covariates = c(rem_forest = "forest", rem_observer = "Observer"),
    yearly_covariates = c(rem_wind = "survey_wind"))
  obs <- unmarked::obsCovs(out$data)
  formulas <- list(lambda = ~forest, phi = ~survey_wind, dist = ~forest,
    rem = ~ rem_forest * factor(interval) + I(rem_wind^2) + (1 | rem_observer))
  design <- get("getDesign", asNamespace("unmarked"))(out$data, formulas)
  fixed <- stats::model.matrix(~rem_forest * factor(interval) + I(rem_wind^2), obs)
  expect_equal(unname(design$X_rem), unname(fixed), tolerance = 0)
  expected_Z <- get("get_Z", asNamespace("unmarked"))(formulas$rem, obs)
  expect_equal(design$Z_rem, expected_Z, tolerance = 0)
  inputs <- get("get_ranef_inputs", asNamespace("unmarked"))(formulas,
    list(lambda = unmarked::siteCovs(out$data), phi = unmarked::yearlySiteCovs(out$data),
      dist = unmarked::siteCovs(out$data), rem = obs),
    design[c("X_lambda", "X_phi", "X_dist", "X_rem")],
    design[c("Z_lambda", "Z_phi", "Z_dist", "Z_rem")])
  expect_identical(inputs$rand_ef, "b_rem")
  expect_length(inputs$pars$beta_rem, ncol(fixed))
  expect_equal(inputs$data$Z_rem, expected_Z, tolerance = 0)
  expect_length(inputs$pars$lsigma_rem, 1L)
  # The actual provider evaluates the fixed fields when building random terms.
  only_group <- obs["rem_observer"]
  expect_error(get("get_reTrms", asNamespace("unmarked"))(formulas$rem, only_group), "not found")
  expect_false(out$fit_eligibility_established)
})

test_that("no-intercept interactions and aliases preserve distinct component roles", {
  f <- gdr_expansion_fixture()
  out <- expand_gdr_intervals(f$frame, f$keys, "survey", "minutes",
    site_covariates = c(rem_forest = "forest", rem_observer = "Observer", copy = "forest"))
  obs <- unmarked::obsCovs(out$data)
  formula <- ~0 + rem_forest:factor(interval) + log(1 + rem_forest) + rem_observer
  forms <- list(lambda = ~forest, phi = ~1, dist = ~forest, rem = formula)
  design <- get("getDesign", asNamespace("unmarked"))(out$data, forms)
  expect_equal(unname(design$X_rem), unname(stats::model.matrix(formula, obs)), tolerance = 0)
  expect_identical(unmarked::siteCovs(out$data)$forest, unmarked::siteCovs(f$frame)$forest)
  expect_identical(obs$copy, obs$rem_forest)
  expect_length(out$protocol$shared_site_or_annual_names, 0L)
  expect_false("rem_forest" %in% names(unmarked::siteCovs(out$data)))
})

test_that("empty selections and missing covariates retain original effort and missingness", {
  f <- gdr_expansion_fixture()
  empty <- expand_gdr_intervals(f$frame, f$keys, "survey", "minutes")
  expect_identical(names(unmarked::obsCovs(empty$data)), c("recorded", "interval"))
  expect_identical(names(empty$protocol$expanded_covariates), c("source_table", "source", "target"))
  expect_equal(nrow(empty$protocol$expanded_covariates), 0)
  f$frame@siteCovs$forest[1] <- NA_real_
  f$frame@yDistance[1, ] <- NA_real_
  f$frame@yRemoval[1, ] <- NA_real_
  f$frame@y[1, ] <- NA_real_
  out <- expand_gdr_intervals(f$frame, f$keys, "survey", "minutes",
    site_covariates = "forest", effort = c(0, 4))
  expect_true(all(is.na(unmarked::obsCovs(out$data)$forest[1:3])))
  expect_identical(out$data@yDistance, f$frame@yDistance)
  expect_identical(out$data@yRemoval, f$frame@yRemoval)
  expect_true(out$survey_audit$valid)
  expect_false(out$fit_eligibility_established)
})

test_that("selection, collision and event-order failures cannot overwrite or coerce inputs", {
  f <- gdr_expansion_fixture()
  original <- f$frame
  call <- function(...) expand_gdr_intervals(f$frame, f$keys, "survey", "minutes", ...)
  for (bad in list(NA_character_, "unknown", 1, list("forest"), ""))
    expect_error(call(site_covariates = bad), "existing covariate")
  expect_error(call(site_covariates = stats::setNames("forest", NA_character_)), "aliases")
  expect_error(call(site_covariates = c(interval = "forest")), "distinct names")
  expect_error(call(site_covariates = c(survey = "forest")), "distinct names")
  expect_error(call(site_covariates = c(recorded = "forest")), "cannot be overwritten")
  expect_error(call(site_covariates = c(same = "forest"), yearly_covariates = c(same = "survey_wind")), "distinct names")
  expect_error(call(site_covariates = c("forest", "forest")), "distinct names")
  expect_error(call(interval = "survey"), "distinct names")
  expect_error(expand_gdr_intervals(f$frame, f$keys[-1, ], "survey", "minutes"), "complete bounded")
  expect_error(expand_gdr_intervals(f$frame, f$keys[2:1, ], "survey", "minutes",
    expected_event_keys = f$keys), "survey audit has errors")
  expect_error(expand_gdr_intervals(f$frame, f$keys, "missing", "minutes"), "explicit event_id")
  expect_error(expand_gdr_intervals(list(), f$keys, "survey", "minutes"), "unmarkedFrameGDR")
  expect_identical(f$frame, original)
})
