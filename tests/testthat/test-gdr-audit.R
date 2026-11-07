gdr_audit_fixture <- function() {
  frame <- unmarked::unmarkedFrameGDR(
    matrix(c(1, 2, 0, 0, 2, 1), 3, 2, byrow = TRUE),
    matrix(c(0, 1, 2, 0, 0, 0, 1, 1, 1), 3, 3, byrow = TRUE),
    siteCovs = data.frame(event = factor(c("a", "b", "c")), category = factor(c("None", "wind", "None"))),
    dist.breaks = c(0, 25, 75), unitsIn = "m", period.lengths = c(1, 2, 1))
  keys <- data.frame(event = c("a", "b", "c"), site = c("A", "A", "B"))
  list(frame = frame, keys = keys)
}

test_that("GDR audits preserve surveys, repeated physical sites, and literal None", {
  f <- gdr_audit_fixture(); before <- serialize(f$frame, NULL)
  a <- validate_unmarked_data(f$frame, f$keys, "event", "site", "minutes", c(4, 4, 4), f$keys)
  expect_s3_class(a, "msr_unmarked_audit")
  expect_true(a$valid)
  expect_equal(nrow(a$findings), 0)
  expect_identical(a$events, f$keys)
  expect_identical(a$protocol$interval_units, "minutes")
  expect_identical(serialize(f$frame, NULL), before)
  expect_identical(levels(f$frame@siteCovs$category), c("None", "wind"))
  expect_identical(as.character(f$frame@siteCovs$category), c("None", "wind", "None"))
  expect_equal(a$summary$events, 3)
})

test_that("equal paired totals cannot hide duplicate original event keys", {
  f <- gdr_audit_fixture(); f$keys$event[2] <- "a"
  f$frame@siteCovs$event <- f$keys$event
  expect_equal(rowSums(f$frame@yDistance), rowSums(f$frame@yRemoval))
  a <- validate_unmarked_data(f$frame, f$keys, "event", interval_units = "minutes", effort = rep(4, 3))
  expect_false(a$valid)
  expect_equal(a$findings$.event[a$findings$code == "duplicate_event"], c(1, 2))
  index <- data.frame(index = 1:3)
  a <- validate_unmarked_data(f$frame, index, "index", interval_units = "minutes")
  expect_true("index_like_event_key" %in% a$findings$code)
  expect_true("effort_unspecified" %in% a$findings$code)
})

test_that("GDR audit catches reordered, unmatched, missing, and misaligned keys", {
  f <- gdr_audit_fixture()
  audit <- function(keys = f$keys, frame = f$frame, reference = f$keys)
    validate_unmarked_data(frame, keys, "event", interval_units = "minutes", effort = rep(4, 3),
                           expected_event_keys = reference)
  a <- audit(f$keys[c(3, 2, 1), ])
  expect_equal(a$findings$.event[a$findings$code == "event_alignment"], c(1, 3))
  expect_true("covariate_key_alignment" %in% a$findings$code)
  keys <- f$keys; keys$event[3] <- "unknown"
  expect_true("event_alignment" %in% audit(keys)$findings$code)
  keys$event[2] <- NA
  expect_true("missing_event_key" %in% audit(keys)$findings$code)
  expect_true("event_layout" %in% audit(f$keys[1:2, ])$findings$code)
  reference <- f$keys; reference$event[2] <- "a"
  expect_true("duplicate_reference_event" %in% audit(reference = reference)$findings$code)
  bad <- f$frame; bad@siteCovs <- bad@siteCovs[1:2, ]
  expect_true("covariate_layout" %in% audit(frame = bad)$findings$code)
})

test_that("GDR audit checks counts, protocol, and inherited response alignment", {
  f <- gdr_audit_fixture()
  codes <- function(x) validate_unmarked_data(x, f$keys, "event", interval_units = "minutes")$findings$code
  for (value in c(-1, 0.5, Inf, NaN)) {
    x <- f$frame; x@yDistance[1, 1] <- value
    expect_true("invalid_counts" %in% codes(x))
  }
  x <- f$frame; x@yDistance[1, 1] <- 2
  expect_true("unequal_totals" %in% codes(x))
  x <- f$frame; x@yRemoval[1, 1] <- 1
  expect_true("response_slot_alignment" %in% codes(x))
  x <- f$frame; x@dist.breaks <- c(1, 25, 75)
  expect_true("distance_breaks" %in% codes(x))
  x <- f$frame; x@period.lengths[2] <- 1.5
  expect_true("interval_lengths" %in% codes(x))
  x <- f$frame; x@unitsIn <- "feet"
  expect_true("distance_units" %in% codes(x))
  x <- f$frame; x@survey <- "line"
  expect_true("survey_type" %in% codes(x))
  x <- f$frame; x@yDistance <- matrix("not_counts", 3, 2)
  expect_true("frame_layout" %in% codes(x))
  x <- f$frame; x@numPrimary <- 2
  expect_true("frame_layout" %in% codes(x))
})

test_that("effort distinguishes a surveyed zero from a missing history", {
  f <- gdr_audit_fixture()
  audit <- function(x, effort) validate_unmarked_data(x, f$keys, "event",
    interval_units = "minutes", effort = effort)
  expect_true(audit(f$frame, rep(4, 3))$valid)
  expect_true("unsurveyed_counts" %in% audit(f$frame, c(4, 0, 4))$findings$code)
  x <- f$frame; x@yDistance[2, ] <- NA; x@yRemoval[2, ] <- x@y[2, ] <- NA
  expect_true(audit(x, c(4, 0, 4))$valid)
  expect_equal(audit(x, c(4, 0, 4))$summary$missing_histories, 1)
  expect_true("surveyed_missing" %in% audit(x, rep(4, 3))$findings$code)
  expect_true("unknown_effort" %in% audit(x, c(4, NA, 4))$findings$code)
  x@yDistance[2, 1] <- 0
  expect_true("partial_history" %in% audit(x, c(4, 0, 4))$findings$code)
  expect_true("invalid_effort" %in% audit(f$frame, c(4, -1, 4))$findings$code)
})

test_that("GDR audit maps all primary and observation rows explicitly", {
  yd <- matrix(c(1, 0, 0, 0, 0, 1, 0, 0), 2, 4, byrow = TRUE)
  keys <- data.frame(site = factor(c("A", "A", "B", "B")), visit = c(1, 2, 1, 2))
  frame <- unmarked::unmarkedFrameGDR(yd, yd, numPrimary = 2,
    siteCovs = data.frame(site = factor(c("A", "B"))), yearlySiteCovs = keys,
    obsCovs = keys[rep(1:4, each = 2), ], dist.breaks = c(0, 25, 75), unitsIn = "m")
  audit <- function(x) validate_unmarked_data(x, keys, c("site", "visit"), "site", "minutes", rep(2, 4))
  expect_true(audit(frame)$valid)
  frame@obsCovs$visit[6] <- 2
  a <- audit(frame)
  expect_equal(a$findings$.event, 3)
  expect_equal(a$findings$.site_row, 2)
  expect_equal(a$findings$.primary, 1)
  frame@yearlySiteCovs$visit[1] <- 2
  expect_equal(sum(audit(frame)$findings$code == "covariate_key_alignment"), 2)
  changed_keys <- keys; changed_keys$site[2] <- "B"
  a <- validate_unmarked_data(frame, changed_keys, c("site", "visit"), "site", "minutes", rep(2, 4))
  expect_equal(a$findings$.event[a$findings$code == "physical_site_alignment"], c(1, 2))
})

test_that("GDR audit rejects ambiguous caller contracts", {
  f <- gdr_audit_fixture()
  expect_error(validate_unmarked_data(list(), f$keys, "event", interval_units = "minutes"), "unmarkedFrameGDR")
  expect_error(validate_unmarked_data(f$frame, f$keys, "absent", interval_units = "minutes"), "event_id")
  expect_error(validate_unmarked_data(f$frame, f$keys, "event", interval_units = ""), "interval_units")
  expect_error(validate_unmarked_data(f$frame, f$keys, "event", interval_units = "minutes", effort = 1), "effort")
  expect_error(validate_unmarked_data(f$frame, f$keys, "event", interval_units = "minutes",
    expected_event_keys = f$keys[1:2, ]), "expected_event_keys")
})

test_that("audit printing is compact without dropping the complete finding records", {
  frame <- unmarked::unmarkedFrameGDR(matrix(0, 12, 2), matrix(0, 12, 3),
    dist.breaks = c(0, 25, 75), unitsIn = "m")
  a <- validate_unmarked_data(frame, data.frame(event = rep("duplicate", 12)), "event",
    interval_units = "minutes", effort = rep(3, 12))
  text <- utils::capture.output(returned <- print(a))
  expect_match(paste(text, collapse = "\n"), "first 10 of 12")
  expect_equal(nrow(a$findings), 12)
  expect_identical(returned, a)
})
