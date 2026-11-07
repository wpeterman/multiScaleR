gdr_interval_fixture <- function(primary = 1L) {
  events <- 2L * primary
  n <- events * 3L
  frame <- unmarked::unmarkedFrameGDR(matrix(0, 2, 2 * primary),
    matrix(0, 2, 3 * primary), numPrimary = primary,
    siteCovs = data.frame(physical_site = factor(c("A", "B")), forest = c(.2, .8)),
    yearlySiteCovs = data.frame(season_cov = seq_len(events)),
    obsCovs = data.frame(existing = seq_len(n)),
    dist.breaks = c(0, 25, 75), unitsIn = "m", period.lengths = c(1, 2, 1))
  keys <- data.frame(physical_site = rep(c("A", "B"), each = primary),
    survey = paste0("survey", seq_len(events)))
  bands <- keys[rep(seq_len(events), each = 3), , drop = FALSE]
  rownames(bands) <- paste0("record", seq_len(n))
  bands$interval <- rep(1:3, events)
  bands$observer <- factor(rep(c("None", "A", "B"), events),
    levels = c("None", "A", "B", "unused"))
  contrasts(bands$observer) <- stats::contr.sum(4)
  bands$ordered_band <- ordered(rep(c("low", "mid", "high"), events))
  bands[["wind m/s"]] <- seq_len(n) / 10
  bands$mat <- I(cbind(left = seq_len(n), right = seq_len(n) + 7))
  list(frame = frame, keys = keys, bands = bands,
    event_id = c("physical_site", "survey"), effort = rep(4, events))
}

gdr_interval_prepare_fixture <- function(f, bands = f$bands, ...) {
  prepare_gdr_intervals(f$frame, f$keys, f$event_id, bands,
    interval_units = "minutes", site_id = "physical_site", effort = f$effort, ...)
}

test_that("keyed GDR intervals preserve all original slots and column attributes", {
  f <- gdr_interval_fixture()
  original <- f$frame
  rows <- c(6, 2, 4, 1, 5, 3)
  input <- f$bands[rows, , drop = FALSE]
  prepared <- gdr_interval_prepare_fixture(f, input)
  expect_s3_class(prepared, "msr_gdr_intervals")
  expect_true(prepared$survey_audit$valid)
  expect_false(prepared$fit_eligibility_established)
  expect_true(prepared$protocol$reordered)
  expect_equal(prepared$alignment$.input_row, match(seq_len(6), rows))
  expect_identical(prepared$alignment$.input_row_name, rownames(f$bands))
  expect_equal(prepared$alignment$.event, rep(1:2, each = 3))
  expect_equal(prepared$alignment$.site_row, rep(1:2, each = 3))
  expect_equal(prepared$alignment$.primary, rep(1L, 6))
  expect_equal(prepared$alignment$.interval, rep(1:3, 2))
  for (slot in setdiff(methods::slotNames(f$frame), "obsCovs"))
    expect_identical(methods::slot(prepared$data, slot), methods::slot(original, slot))
  actual <- unmarked::obsCovs(prepared$data)
  expect_identical(actual$existing, unmarked::obsCovs(original)$existing)
  for (name in setdiff(names(f$bands), f$event_id))
    expect_identical(actual[[name]], f$bands[[name]])
  expect_identical(levels(actual$observer), c("None", "A", "B", "unused"))
  expect_identical(contrasts(actual$observer), contrasts(f$bands$observer))
  expect_false(anyNA(actual$observer))
  expect_identical(prepared$protocol$interval_lengths, c(1, 2, 1))
  expect_identical(prepared$protocol$interval_units, "minutes")
  expect_identical(f$frame, original)
  printed <- capture.output(returned <- print(prepared))
  expect_lt(length(printed), 8L)
  expect_true(any(grepl("Input order changed", printed, fixed = TRUE)))
  expect_identical(returned, prepared)
  expect_false(gdr_interval_prepare_fixture(f)$protocol$reordered)
})

test_that("native GDR designs retain interval formulas and primary/site placement", {
  design <- get("getDesign", asNamespace("unmarked"))
  for (primary in 1:2) {
    f <- gdr_interval_fixture(primary)
    n <- nrow(f$bands)
    prepared <- gdr_interval_prepare_fixture(f, f$bands[rev(seq_len(n)), , drop = FALSE])
    canonical <- f$frame
    canonical_covs <- data.frame(unmarked::obsCovs(f$frame),
      f$bands[setdiff(names(f$bands), f$event_id)], check.names = FALSE)
    rownames(canonical_covs) <- rownames(unmarked::obsCovs(f$frame))
    unmarked::obsCovs(canonical) <- canonical_covs
    formulas <- list(lambda = ~forest, phi = ~season_cov, dist = ~season_cov,
      rem = ~I(interval == 1) + observer * `wind m/s` + ordered_band + mat)
    for (removal in list(formulas$rem, ~0 + factor(interval) + I(`wind m/s`^2),
        ~., ~observer / ordered_band, ~I(interval == 1) + (1 | observer))) {
      formulas$rem <- removal
      actual <- design(prepared$data, formulas)
      expected <- design(canonical, formulas)
      expect_identical(actual, expected)
      expect_equal(nrow(actual$X_rem), n)
      expect_equal(nrow(actual$X_lambda), 2L)
      expect_equal(nrow(actual$X_dist), 2L * primary)
      expect_identical(prepared$data@yDistance, f$frame@yDistance)
      expect_identical(prepared$data@yRemoval, f$frame@yRemoval)
    }
    expect_equal(prepared$alignment$.primary,
      rep(rep(seq_len(primary), each = 3), 2))
    expect_equal(prepared$alignment$.site_row, rep(1:2, each = 3 * primary))
  }
})

test_that("compound event identifiers do not collide at embedded separators or long labels", {
  f <- gdr_interval_fixture()
  f$keys$left <- c("a\rb", "a")
  f$keys$right <- c("c", "b\rc")
  f$event_id <- c("left", "right")
  bands <- f$keys[rep(1:2, each = 3), f$event_id, drop = FALSE]
  bands$interval <- rep(1:3, 2)
  bands$label <- factor(rep(c("None", "yes", "no"), 2))
  actual <- gdr_interval_prepare_fixture(f, bands[c(6, 2, 4, 1, 5, 3), ])
  expect_identical(unmarked::obsCovs(actual$data)$label, bands$label)
  expect_equal(actual$alignment$left, rep(f$keys$left, each = 3))
  f$keys$left <- c(paste0(rep("long", 3000), collapse = ""), "\u03b1")
  bands$left <- rep(f$keys$left, each = 3)
  actual <- gdr_interval_prepare_fixture(f, bands[6:1, ])
  expect_identical(unmarked::obsCovs(actual$data)$label, bands$label)
  expect_equal(actual$alignment$left, rep(f$keys$left, each = 3))
})

test_that("invalid interval keys cannot silently drop, duplicate or overwrite records", {
  f <- gdr_interval_fixture()
  original <- f$frame
  expect_error(gdr_interval_prepare_fixture(f, f$bands[-1, ]), "exactly one row")
  expect_error(gdr_interval_prepare_fixture(f, rbind(f$bands, f$bands[1, ])), "exactly one row")
  duplicate <- f$bands
  duplicate$interval[2] <- duplicate$interval[1]
  condition <- tryCatch(gdr_interval_prepare_fixture(f, duplicate), error = identity)
  expect_s3_class(condition, "msr_gdr_interval_error")
  expect_identical(condition$input_rows, c(1L, 2L))
  expect_true(condition$survey_audit$valid)
  unknown <- f$bands
  unknown$survey[1] <- "unknown"
  expect_error(gdr_interval_prepare_fixture(f, unknown), "unknown survey")
  unknown$survey[1] <- NA_character_
  expect_error(gdr_interval_prepare_fixture(f, unknown), "must be complete")
  for (bad in c(0, 4, 1.5, Inf, NA_real_)) {
    bands <- f$bands
    bands$interval[1] <- bad
    expect_error(gdr_interval_prepare_fixture(f, bands), "native removal range")
  }
  overwrite <- f$bands
  overwrite$existing <- 99
  expect_error(gdr_interval_prepare_fixture(f, overwrite), "cannot be overwritten")
  expect_error(prepare_gdr_intervals(f$frame, f$keys, f$event_id, f$bands,
    interval_units = "minutes", interval = "survey"), "different name")
  reordered <- f$keys[2:1, ]
  expect_error(prepare_gdr_intervals(f$frame, reordered, f$event_id, f$bands,
    interval_units = "minutes", expected_event_keys = f$keys), "survey audit has errors")
  expect_identical(f$frame, original)
})

test_that("missing histories, literal categories and shared native names stay explicit", {
  f <- gdr_interval_fixture()
  f$frame@yDistance[1, ] <- NA_real_
  f$frame@yRemoval[1, ] <- NA_real_
  f$frame@y[1, ] <- NA_real_
  f$effort[1] <- 0
  f$bands$unknown <- c(rep(NA_real_, 3), 1:3)
  prepared <- gdr_interval_prepare_fixture(f, f$bands[6:1, ])
  expect_true(prepared$survey_audit$valid)
  expect_true(all(is.na(prepared$data@yDistance[1, ])))
  expect_true(all(is.na(prepared$data@yRemoval[1, ])))
  expect_true(all(is.na(unmarked::obsCovs(prepared$data)$unknown[1:3])))
  expect_identical(unmarked::obsCovs(prepared$data)$observer, f$bands$observer)
  bad <- f
  bad$frame@yDistance[2, 1] <- 1
  condition <- tryCatch(gdr_interval_prepare_fixture(bad), error = identity)
  expect_s3_class(condition, "msr_gdr_interval_error")
  expect_false(condition$survey_audit$valid)
  expect_true("unequal_totals" %in% condition$survey_audit$findings$code)
  f <- gdr_interval_fixture()
  f$bands$forest <- seq_len(6)
  prepared <- gdr_interval_prepare_fixture(f)
  expect_identical(prepared$protocol$shared_site_or_annual_names, "forest")
  expect_identical(unmarked::siteCovs(prepared$data), unmarked::siteCovs(f$frame))
  expect_false(prepared$fit_eligibility_established)
})
