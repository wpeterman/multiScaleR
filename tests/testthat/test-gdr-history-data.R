gdr_history_fixture <- function() {
  map <- expand.grid(interval = 1:5, distance_bin = 1:3)
  map$column <- paste0("source ", seq_len(15), " (mark)")
  h <- matrix(0, 10, 15, dimnames = list(paste0("birdrow", 1:10), map$column))
  at <- function(time, distance) (distance - 1L) * 5L + time
  h[1, c(at(1, 3), at(3, 1))] <- 1 # Outside first, near later.
  h[2, at(2, 1)] <- 1
  h[2, at(5, 2)] <- NA             # Valid first, malformed later.
  h[3, at(1, 2)] <- 1
  h[3, c(at(4, 1), at(4, 2))] <- 1 # Contradictory later distances.
  h[4, c(at(3, 2), at(5, 3))] <- 1 # Preserve later outside-radius mark.
  h[6, at(4, 1)] <- 1
  h[6, at(1, 3)] <- NA             # Earlier uncertainty quarantines event C.
  h[7, at(1, 1)] <- 1
  h[8, c(at(1, 1), at(1, 2))] <- 1 # Ambiguous first distance.
  h[9, ] <- NA                    # Unsurveyed, not a zero-count event.
  h[10, at(2, 1)] <- 1
  keys <- data.frame(survey = factor(c(rep("A", 4), "B", "C", "C", "D", "E", "F"),
      levels = c(LETTERS[1:7], "unused")),
    bird = c("outside", "robin_later_bad", "later_duplicate", "valid", "no_mark",
      "earlier_unknown", "valid", "first_duplicate", "unknown", "valid"),
    observer = factor(rep(c("None", "Ann"), 5), levels = c("None", "Ann", "unused")))
  contrasts(keys$observer) <- stats::contr.sum(3)
  events <- data.frame(survey = factor(LETTERS[1:7], levels = levels(keys$survey)))
  list(histories = h, history_keys = keys, event_keys = events, event_id = "survey",
    individual_id = "bird", cell_map = map, dist_breaks = c(0, 25, 75, Inf),
    radius = 75, interval_lengths = rep(75, 5), interval_units = "seconds",
    units_in = "m", surveyed = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE))
}

test_that("first encounters precede radius filtering and preserve malformed later records", {
  f <- gdr_history_fixture()
  original <- f
  out <- do.call(prepare_gdr_histories, f)
  expect_s3_class(out, "msr_gdr_histories")
  for (n in c("histories", "history_keys", "event_keys", "cell_map", "surveyed"))
    expect_identical(out[[n]], f[[n]])
  expect_identical(f, original)
  expect_identical(out$history_keys$observer, f$history_keys$observer)
  expect_equal(out$first$first_interval[c(1, 2, 3, 4, 10)], c(1, 2, 1, 3, 2))
  expect_equal(out$first$first_distance_bin[c(1, 2, 3, 4)], c(3, 1, 2, 2))
  expect_identical(out$first$status[1], "first_outside_radius")
  expect_false(out$first$primary_eligible[1])
  expect_true(out$first$primary_eligible[2])
  expect_false(out$first$complete_history_eligible[2])
  expect_true(out$first$primary_eligible[3])
  expect_false(out$first$history_complete[3])
  expect_true(out$first$complete_history_eligible[4])
  expect_equal(out$yDistance[1, ], c(1, 2))
  expect_equal(out$yRemoval[1, ], c(1, 1, 1, 0, 0))
  expect_equal(sum(out$first_cells[1, ]), 3)
  expect_equal(out$summary$primary_records, 5)
  expect_equal(out$summary$complete_histories, 3)
  expect_false(out$fit_eligibility_established)
  expect_false(out$density_identification_established)
  expect_match(out$fingerprint, "^[a-f0-9]{64}$")
  printed <- capture.output(returned <- print(out))
  expect_lte(length(printed), 5)
  expect_identical(returned, out)
})

test_that("unresolved first encounters quarantine whole events without inventing zeros", {
  out <- do.call(prepare_gdr_histories, gdr_history_fixture())
  expect_true(all(is.na(out$yDistance[c(3, 4, 5), ])))
  expect_true(all(is.na(out$yRemoval[c(3, 4, 5), ])))
  expect_true(all(is.na(out$first_cells[c(3, 4, 5), ])))
  expect_equal(out$yDistance[c(2, 7), ], matrix(0, 2, 2))
  expect_equal(out$yRemoval[c(2, 7), ], matrix(0, 2, 5))
  expect_equal(out$events$known_first_records, c(3, 0, 1, 0, 0, 1, 0))
  expect_equal(out$events$response_complete, c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_equal(out$events$unresolved_first_histories, c(0, 0, 1, 1, 1, 0, 0))
  expect_equal(out$summary$unresolved_events, 2)
  expect_true(is.na(out$first$first_interval[6]))
  expect_equal(out$first$candidate_interval[6], 4)
  expect_identical(out$first$status[8], "unresolved_first")
  expect_identical(out$first$status[9], "unsurveyed_history")
  expect_true(all(c("invalid_mark", "multiple_distance_bins", "first_outside_radius",
    "unresolved_first", "unsurveyed_history") %in% out$findings$code))
  for (j in 1:2)
    expect_equal(out$yDistance[, j], rowSums(out$first_cells[, (j - 1) * 5 + 1:5, drop = FALSE]))
  for (j in 1:5)
    expect_equal(out$yRemoval[, j], rowSums(out$first_cells[, c(j, j + 5), drop = FALSE]))
})

test_that("explicit maps preserve source permutations, identifiers and empty surveys", {
  f <- gdr_history_fixture()
  expected <- do.call(prepare_gdr_histories, f)
  row_order <- c(10, 8, 3, 6, 4, 2, 1, 9, 5, 7)
  col_order <- c(15:9, 1:8)
  f$histories <- f$histories[row_order, col_order]
  f$history_keys <- f$history_keys[row_order, , drop = FALSE]
  f$cell_map <- f$cell_map[c(9:15, 1:8), , drop = FALSE]
  out <- do.call(prepare_gdr_histories, f)
  expect_identical(out$histories, f$histories)
  expect_identical(out$cell_map, f$cell_map)
  expect_identical(out$first_cells, expected$first_cells)
  expect_identical(out$yDistance, expected$yDistance)
  expect_identical(out$yRemoval, expected$yRemoval)
  expect_equal(out$first$first_interval, expected$first$first_interval[row_order])
  f$histories <- f$histories[FALSE, , drop = FALSE]
  f$history_keys <- f$history_keys[FALSE, , drop = FALSE]
  empty <- do.call(prepare_gdr_histories, f)
  expect_equal(nrow(empty$first), 0)
  expect_equal(sum(empty$yDistance, na.rm = TRUE), 0)
  expect_true(all(is.na(empty$yDistance[5, ])))
  expect_identical(empty$history_keys, f$history_keys)
  # Compound keys are not concatenated with an ambiguous separator.
  f <- gdr_history_fixture()
  f$event_keys$left <- c("a\rb", "a", "C", "D", "E", "F", "G")
  f$event_keys$right <- c("c", "b\rc", rep("key", 5))
  event <- match(as.character(f$history_keys$survey), as.character(f$event_keys$survey))
  f$history_keys$left <- f$event_keys$left[event]
  f$history_keys$right <- f$event_keys$right[event]
  f$event_id <- c("left", "right")
  out <- do.call(prepare_gdr_histories, f)
  expect_identical(out$first$.event, expected$first$.event)
  f$event_keys$left[1] <- paste(rep("long", 3000), collapse = "")
  f$history_keys$left[event == 1] <- f$event_keys$left[1]
  expect_identical(do.call(prepare_gdr_histories, f)$yDistance, expected$yDistance)
})

test_that("malformed cells and structural mapping conflicts stay explicit", {
  f <- gdr_history_fixture()
  for (bad in c(2, -1, 1.5, Inf, NaN)) {
    changed <- f
    changed$histories[2, 15] <- bad # Later than the valid first.
    out <- do.call(prepare_gdr_histories, changed)
    expect_identical(out$histories[2, 15], bad)
    expect_true(out$first$primary_eligible[2])
    expect_false(out$first$history_complete[2])
    expect_true(any(out$findings$.history_row == 2 & out$findings$.cell == 15 & out$findings$code == "invalid_mark"))
  }
  bad <- f
  bad$surveyed[1] <- FALSE
  condition <- tryCatch(do.call(prepare_gdr_histories, bad), error = identity)
  expect_s3_class(condition, "msr_gdr_history_error")
  expect_true(all(1:4 %in% condition$history_rows))
  for (bad_status in list(NULL, rep(TRUE, 2), c(TRUE, rep(NA, 6)), rep(1, 7))) {
    bad <- f
    bad["surveyed"] <- list(bad_status)
    expect_error(do.call(prepare_gdr_histories, bad), "explicit complete logical")
  }
  for (radius in c(0, 40, Inf)) {
    bad <- f; bad$radius <- radius
    expect_error(do.call(prepare_gdr_histories, bad), "radius must match")
  }
  bad <- f; bad$cell_map$interval[1] <- 2
  expect_error(do.call(prepare_gdr_histories, bad), "combination exactly once")
  bad <- f; bad$cell_map$column[1] <- "unknown"
  expect_error(do.call(prepare_gdr_histories, bad), "every named history column")
  bad <- f
  bad$history_keys$survey <- as.character(f$history_keys$survey)
  bad$history_keys$survey[1] <- "unknown"
  expect_error(do.call(prepare_gdr_histories, bad), "unknown event")
  bad <- f; bad$history_keys[2, c("survey", "bird")] <- bad$history_keys[1, c("survey", "bird")]
  expect_error(do.call(prepare_gdr_histories, bad), "identifiers must be unique")
  bad <- f; bad$expected_event_keys <- f$event_keys[7:1, , drop = FALSE]
  expect_error(do.call(prepare_gdr_histories, bad), "independent event mapping")
  bad <- f; storage.mode(bad$histories) <- "character"
  expect_error(do.call(prepare_gdr_histories, bad), "numeric or logical matrix")
})

test_that("native optional validation checks exact margins and finite physical geometry", {
  f <- gdr_history_fixture()
  retained <- c(1:5, 10)
  f$histories <- f$histories[retained, , drop = FALSE]
  f$history_keys <- f$history_keys[retained, , drop = FALSE]
  f$event_keys <- f$event_keys[c(1, 2, 6, 7), , drop = FALSE]
  f$surveyed <- rep(TRUE, 4)
  f$interval_lengths <- c(1, 2, 1, 3, 1)
  expected <- do.call(prepare_gdr_histories, f)
  frame <- unmarked::unmarkedFrameGDR(expected$yDistance, expected$yRemoval,
    dist.breaks = c(0, .025, .075), unitsIn = "km", period.lengths = f$interval_lengths)
  f$native_frame <- frame
  f$expected_event_keys <- f$event_keys
  out <- do.call(prepare_gdr_histories, f)
  expect_true(out$survey_audit$valid)
  expect_identical(out$protocol$interval_lengths, f$interval_lengths)
  expect_equal(rowSums(out$yDistance), rowSums(out$yRemoval))
  expect_identical(out$protocol$response_conventions[["first_cells"]], "paired_first_distance_time")
  expect_identical(out$protocol$response_conventions[["native"]], "native_two_margins")
  f$native_frame@yDistance[1, ] <- rev(frame@yDistance[1, ])
  expect_error(do.call(prepare_gdr_histories, f), "exactly match")
  f$native_frame <- frame
  f$native_frame@period.lengths[2] <- 3
  expect_error(do.call(prepare_gdr_histories, f), "exactly match")
  f$native_frame <- frame
  f$native_frame@dist.breaks[2] <- .03
  expect_error(do.call(prepare_gdr_histories, f), "exactly match")
  f$native_frame <- frame
  expect_identical(frame@yDistance, expected$yDistance)
  expect_false(out$fit_eligibility_established)
})

test_that("GDR comparison canonicalizes equivalent geometry without changing native inputs", {
  model <- gdr_guard_fixture()$P
  original <- gdr_guard_state(model)
  meters <- .msr_selection_signature(model)
  km <- model
  km@data@unitsIn <- "km"
  km@data@dist.breaks <- km@data@dist.breaks / 1000
  legacy <- .msr_selection_signature(km)
  original_legacy <- legacy
  expect_identical(legacy$gdr_protocol$distance_units, "km")
  expect_null(.msr_assert_selection_signatures(list(meters, legacy)))
  expect_null(.msr_assert_selection_signatures(list(legacy, meters)))
  expect_identical(legacy, original_legacy)
  for (bad in list(NULL, c("0", "25", "75"), c(0, NA, 75), matrix(c(0, 25, 75)))) {
    malformed <- legacy
    malformed$gdr_protocol$distance_breaks <- bad
    expect_error(.msr_assert_selection_signatures(list(meters, malformed)), "distance bins are not qualified")
  }
  expect_identical(gdr_guard_state(model), original)
  expect_identical(km@data@unitsIn, "km")
  expect_equal(km@data@dist.breaks, c(0, .025, .075))
  expect_equal(.msr_parameter_count(km), .msr_parameter_count(model))
  km@data@dist.breaks[2] <- .03
  expect_error(.msr_assert_selection_signatures(list(meters, .msr_selection_signature(km))), "different GDR observation")
})

test_that("invisible-class and paired-margin identities do not identify corrected density", {
  y <- matrix(c(1, 0, 0, 2, 0, 1, 2, 0, 0, 0, 0, 0), 3, 4, byrow = TRUE)
  yd <- cbind(rowSums(y[, 1:2]), rowSums(y[, 3:4]))
  yr <- cbind(y[, 1] + y[, 3], y[, 2] + y[, 4])
  pd <- rbind(c(.2, .3), c(.1, .4), c(.3, .2))
  pr <- rbind(c(.3, .4), c(.2, .3), c(.4, .2))
  abundance <- c(8, 12, 4)
  n <- rowSums(y)
  constant <- sum(-lgamma(n + 1) - rowSums(lgamma(y + 1)) +
    rowSums(lgamma(yd + 1)) + rowSums(lgamma(yr + 1)))
  for (size in c(Inf, 2)) {
    reference <- .msr_gdr_exact_score_bins(abundance, pd, pr, yd, yr, size, "explicit_fixture")
    mean <- abundance * rowSums(pd) * rowSums(pr)
    total <- if (is.infinite(size)) dpois(n, mean, log = TRUE) else dnbinom(n, mu = mean, size = size, log = TRUE)
    paired <- sum(total + vapply(1:3, function(i) {
      probability <- as.vector(t(outer(pd[i, ], pr[i, ])))
      dmultinom(y[i, ], prob = probability / sum(probability), log = TRUE)
    }, numeric(1)))
    expect_equal(paired - reference$logLik, constant, tolerance = 1e-12)
    for (w in c(.3, .9, .999)) {
      invisible <- .msr_gdr_exact_score_bins(abundance / (1 - w), pd * (1 - w), pr,
        yd, yr, size, "explicit_invisible_class_fixture")
      expect_equal(invisible$logLik, reference$logLik, tolerance = 1e-12)
      expect_equal(invisible$per_event$expected_observed, reference$per_event$expected_observed, tolerance = 1e-12)
      expect_equal(invisible$conditional_distance, reference$conditional_distance, tolerance = 1e-12)
      expect_gt(invisible$per_event$expected_abundance[1], reference$per_event$expected_abundance[1])
    }
  }
})
