#' Convert Linked Bird Records into Paired Distance-Removal Counts
#'
#' Turn a matrix of bird records, one row per individual at a survey, into
#' the distance and removal count tables of a distance-removal model. The
#' function finds each bird's first detection across all distance bands,
#' including bands beyond the survey radius, before it applies the radius. It
#' keeps every source record, reports what happened to each bird, and never
#' deletes or repairs a record.
#'
#' @param histories Numeric or logical matrix with one row per individual
#'   bird at a survey and one named column per source cell (one distance band
#'   in one time interval). Enter 1 where the bird was recorded and 0
#'   elsewhere. Any other value, including `NA`, is kept as it is and reported
#'   in `findings`; it is never treated as 0.
#' @param history_keys Data frame with one row per row of `histories`, in the
#'   same order, holding the original survey and individual identifiers and
#'   any other notes you want to keep.
#' @param event_keys Data frame with one row per survey event, in the order
#'   the count tables should have. Include completed surveys that recorded no
#'   birds. Event keys must be unique.
#' @param event_id Character vector naming the column(s) that identify a
#'   survey in both `event_keys` and `history_keys`.
#' @param individual_id Character vector naming the column(s) of
#'   `history_keys` that identify a bird within a survey. Each survey and
#'   individual combination must be unique.
#' @param cell_map Data frame with columns `column` (a column name of
#'   `histories`), `interval` and `distance_bin` (positive whole numbers). It
#'   must map every column of `histories` to exactly one combination of
#'   interval and distance band, and every combination must appear once.
#'   Columns of `histories` can be in any order.
#' @param dist_breaks Increasing distance boundaries of the source bands,
#'   starting at zero, in `units_in`. Include bands beyond the survey radius
#'   when the records have them. Only the last boundary may be `Inf`.
#' @param radius Survey radius, in `units_in`. It must equal one of
#'   `dist_breaks`, because a band cannot be split without finer distances.
#' @param interval_lengths Positive whole numbers giving the duration of each
#'   time interval, in `interval_units`. Intervals may differ in length.
#' @param interval_units Character string naming the time unit, such as
#'   `"seconds"`. It is recorded, not converted.
#' @param units_in Distance units of `dist_breaks` and `radius`, `"m"` or
#'   `"km"`.
#' @param surveyed Logical vector with one value per event, in event order,
#'   with no `NA`: `TRUE` for a completed survey and `FALSE` for an event that
#'   was not surveyed. A recorded bird at an unsurveyed event is an error.
#' @param expected_event_keys Optional data frame, built independently, with
#'   the `event_id` columns in the intended event order. The function stops if
#'   `event_keys` differs from it.
#' @param native_frame Optional single-visit `unmarkedFrameGDR` built earlier
#'   from these records. It is audited with [validate_unmarked_data()], and
#'   the function stops unless its two count tables, its distance boundaries
#'   up to the radius and its interval lengths match the prepared results
#'   exactly. Distances in meters and kilometers are compared in meters,
#'   without changing either input.
#' @param x A prepared `msr_gdr_histories` object.
#' @param ... Not used.
#' @return A list of class `msr_gdr_histories`:
#'   * `histories`, `history_keys`, `event_keys` and `cell_map`: the inputs,
#'     unchanged.
#'   * `first`: one row per bird, with its row in `histories`
#'     (`.history_row`), its event (`.event`), the first interval with a
#'     recorded mark (`candidate_interval`), the resolved first detection
#'     (`first_interval`, `first_distance_bin`; `NA` when unresolved),
#'     whether the bird is counted (`primary_eligible`), whether every cell
#'     holds a valid 0 or 1 with at most one band per interval
#'     (`history_complete`), whether both hold (`complete_history_eligible`)
#'     and its `status` (see "Read the status of each bird").
#'   * `events`: one row per event, with `surveyed`, `response_complete`
#'     (whether its counts are observed), `known_first_records` (birds
#'     counted) and `unresolved_first_histories`.
#'   * `findings`: one row per problem, with the history row, the cell
#'     (`.cell`, `.column`, `interval`, `distance_bin`), a `code` and the
#'     original `value`. Findings about a whole history have missing cell
#'     fields.
#'   * `first_cells`: counts of first detections by event (rows) and cell
#'     inside the radius (columns), ordered by distance band, then interval.
#'   * `native_cells`: the rows of `cell_map` for those columns, in the same
#'     order.
#'   * `yDistance` and `yRemoval`: the two count tables, events by distance
#'     band and events by interval, ready for `unmarked::unmarkedFrameGDR()`.
#'   * `protocol`: the source and in-radius distance breaks, the radius,
#'     units, interval lengths and the column order of `first_cells`.
#'   * `surveyed`, `survey_audit` (the audit of `native_frame`, or `NULL`),
#'     `fingerprint` (a SHA-256 of the inputs) and `summary` (counts of
#'     histories, events, counted birds, complete histories, invalid cells,
#'     unresolved events and birds first detected outside the radius).
#'   * `fit_eligibility_established` and
#'     `density_identification_established`: always `FALSE`. Preparing counts
#'     does not show that a model fitted to them is valid.
#'
#'   Printing reports the summary counts and returns the object invisibly.
#'   Errors in the inputs stop with a condition of class
#'   `msr_gdr_history_error` whose `history_rows` element lists the rows
#'   involved, when they are known.
#' @details
#' # Find each bird's first detection
#'
#' The first detection is found across all source bands before the radius is
#' applied. A bird first recorded beyond the radius stays out of the counts,
#' even if a later mark places it near the point. If a bird has marks in more
#' than one distance band in the same interval, its distance in that interval
#' is not unique. Such marks are kept for review and never collapsed into one
#' band. They may call for a model of movement or of the recording process.
#'
#' A bird is counted when its first detection is known and inside the radius,
#' even if a later mark is malformed. That bird is then left out of the
#' complete-history subset. When a malformed or multi-band mark falls at or
#' before a bird's first detection, the first detection is unresolved, as it
#' is for a history with a malformed mark and no detection. The whole event's
#' counts then become missing, so that an incomplete count
#' cannot pass as an observed zero. The known partial count stays in
#' `events$known_first_records`, and the source rows stay in `histories`.
#'
#' # Read the status of each bird
#'
#' | `status` | Meaning |
#' |:---|:---|
#' | `retained_first` | Counted. The history is complete |
#' | `retained_first_incomplete_history` | Counted, but a later mark is malformed or spans two bands |
#' | `first_outside_radius` | Not counted. The bird was first detected beyond the radius |
#' | `unresolved_first` | The first detection cannot be determined. The event's counts are missing |
#' | `no_record` | No detections were recorded, and every cell is a valid 0 |
#' | `unsurveyed_history` | The history belongs to an event that was not surveyed |
#'
#' `findings` uses the codes `invalid_mark` (a cell that is not 0 or 1, with
#' its `value`), `multiple_distance_bins` (marks in more than one band in an
#' interval), and `unresolved_first`, `first_outside_radius` and
#' `unsurveyed_history` as in the table.
#'
#' # Read the count tables
#'
#' Summing the columns of `first_cells` over intervals gives `yDistance`, and
#' summing over distance bands gives `yRemoval`, so the two tables always
#' have equal row totals. A completed survey with no counted birds has a row
#' of zeros. An unsurveyed event, or an event with an unresolved bird, has a
#' row of `NA`.
#'
#' # Keep the other responses separate
#'
#' The paired first detections are one response. Complete linked histories of
#' the counted birds, re-detections, and marks beyond the radius are others.
#' Their likelihoods have different constants and conditioning, so they
#' cannot be compared in one AIC table. The time intervals within one count
#' are also not independent visits. Keep the source records for any separate
#' analysis. Better records of the detected birds do not by themselves
#' identify the birds that were never detected.
#'
#' See the "Prepare your own records" section of
#' `vignette("v12-distance-removal", package = "multiScaleR")` for a worked
#' example.
#' @seealso [validate_unmarked_data()], [prepare_gdr_intervals()],
#'   [expand_gdr_intervals()]
#' @examples
#' # Source cells: 2 time intervals by 3 distance bands (0-25 m, 25-75 m and
#' # beyond the 75 m radius)
#' map <- expand.grid(interval = 1:2, distance_bin = 1:3)
#' map$column <- paste0("cell", 1:6)
#' h <- matrix(0, 2, 6, dimnames = list(NULL, map$column))
#' h[1, c(5, 2)] <- 1 # Beyond 75 m first, then near: not counted
#' h[2, 1] <- 1       # Interval 1, 0-25 m: counted
#' h[2, 2] <- 2       # Malformed later mark: the valid first still counts
#' prepared <- prepare_gdr_histories(h,
#'   history_keys = data.frame(survey = c("a", "a"), bird = c("one", "two")),
#'   event_keys = data.frame(survey = c("a", "b", "c")), event_id = "survey",
#'   individual_id = "bird", cell_map = map, dist_breaks = c(0, 25, 75, Inf),
#'   radius = 75, interval_lengths = c(75, 75), interval_units = "seconds",
#'   units_in = "m", surveyed = c(TRUE, TRUE, FALSE))
#' # What happened to each bird
#' prepared$first[, c("first_interval", "first_distance_bin", "status")]
#' # Survey b is an observed zero. Survey c was not surveyed, so it is NA
#' prepared$yDistance
#' prepared$yRemoval
#' # Both tables are totals of the same first-detection cells
#' prepared$first_cells
#' @export
prepare_gdr_histories <- function(histories, history_keys, event_keys,
    event_id, individual_id, cell_map, dist_breaks, radius, interval_lengths,
    interval_units, units_in, surveyed, expected_event_keys = NULL,
    native_frame = NULL) {
  fail <- function(message, rows = integer()) stop(structure(list(message = message,
    call = NULL, history_rows = rows), class = c("msr_gdr_history_error", "error", "condition")))
  if (!is.matrix(histories) || !(is.numeric(histories) || is.logical(histories)) ||
      !ncol(histories) || is.null(colnames(histories)) || anyNA(colnames(histories)) ||
      anyDuplicated(colnames(histories)) || any(!nzchar(colnames(histories))))
    fail("histories must be a numeric or logical matrix with unique named source cells.")
  table_ok <- function(x) is.data.frame(x) && ncol(x) > 0L &&
    !anyNA(names(x)) && !anyDuplicated(names(x)) && all(nzchar(names(x))) &&
    all(vapply(x, function(y) is.atomic(y) && is.null(dim(y)), logical(1)))
  if (!table_ok(event_keys) || !nrow(event_keys) || !table_ok(history_keys) ||
      nrow(history_keys) != nrow(histories)) fail("Supply complete event and history key tables in their original row orders.")
  columns_ok <- function(x, data) is.character(x) && length(x) > 0L &&
    !anyNA(x) && !anyDuplicated(x) && all(x %in% names(data))
  if (!columns_ok(event_id, event_keys) || !columns_ok(event_id, history_keys) ||
      !columns_ok(individual_id, history_keys) || any(individual_id %in% event_id))
    fail("event_id and individual_id must name distinct original identifier columns.")
  encode <- function(data, columns) {
    # paste0 recycles a separator even for zero-length fields. Empty linked
    # tables must remain empty, so completed zero-count surveys can be retained.
    if (!nrow(data)) return(character())
    fields <- lapply(data[columns], function(x) {
      x <- enc2utf8(as.character(x))
      if (anyNA(x) || any(!nzchar(trimws(x)))) fail("Event and individual identifiers must be complete.")
      paste0(nchar(x, type = "bytes"), ":", x)
    })
    do.call(paste0, fields)
  }
  event_codes <- encode(event_keys, event_id)
  if (anyDuplicated(event_codes)) fail("Original event identifiers must be unique.")
  history_codes <- encode(history_keys, c(event_id, individual_id))
  if (anyDuplicated(history_codes)) fail("Event/individual history identifiers must be unique.",
    which(duplicated(history_codes) | duplicated(history_codes, fromLast = TRUE)))
  event <- match(encode(history_keys, event_id), event_codes)
  if (anyNA(event)) fail("History rows contain unknown event identifiers.", which(is.na(event)))
  if (!is.null(expected_event_keys) && (!table_ok(expected_event_keys) ||
      !columns_ok(event_id, expected_event_keys) ||
      !identical(encode(expected_event_keys, event_id), event_codes)))
    fail("Event keys differ from the independent event mapping.")
  m <- nrow(event_keys)
  if (!is.logical(surveyed) || !is.null(dim(surveyed)) || length(surveyed) != m || anyNA(surveyed))
    fail("Supply explicit complete logical surveyed status for every event.")
  if (!is.character(units_in) || length(units_in) != 1L || is.na(units_in) || !units_in %in% c("m", "km"))
    fail("units_in must be m or km.")
  validate_character_scalar(interval_units, "interval_units")
  if (!nzchar(trimws(interval_units))) fail("Supply a nonempty interval_units label.")
  if (!is.numeric(dist_breaks) || !is.null(dim(dist_breaks)) || length(dist_breaks) < 2L ||
      anyNA(dist_breaks) || dist_breaks[1L] != 0 || any(diff(dist_breaks) <= 0) ||
      any(!is.finite(utils::head(dist_breaks, -1L))) || utils::tail(dist_breaks, 1L) == -Inf)
    fail("Source distance breaks must increase from zero; only the final boundary may be infinite.")
  if (!is.numeric(radius) || !is.null(dim(radius)) || length(radius) != 1L || !is.finite(radius) || radius <= 0 ||
      !radius %in% dist_breaks) fail("The finite positive radius must match a source distance boundary.")
  if (!is.numeric(interval_lengths) || !is.null(dim(interval_lengths)) || !length(interval_lengths) ||
      any(!is.finite(interval_lengths)) || any(interval_lengths <= 0) ||
      any(interval_lengths != floor(interval_lengths))) fail("Supply positive integer native interval_lengths.")
  jd <- length(dist_breaks) - 1L
  tt <- length(interval_lengths)
  jj <- match(radius, dist_breaks) - 1L
  cells <- as.double(jd) * tt
  if (cells > .Machine$integer.max || as.double(m) * jj * tt > .Machine$integer.max ||
      ncol(histories) != cells) fail("The full source layout or requested native cell matrix is not representable.")
  if (!is.data.frame(cell_map) || anyNA(names(cell_map)) || any(!nzchar(names(cell_map))) || anyDuplicated(names(cell_map)) ||
      !all(c("column", "interval", "distance_bin") %in% names(cell_map)) ||
      nrow(cell_map) != cells || !is.character(cell_map$column) || !is.null(dim(cell_map$column)) || anyNA(cell_map$column) ||
      anyDuplicated(cell_map$column) || !setequal(cell_map$column, colnames(histories)))
    fail("cell_map must map every named history column exactly once.")
  index_ok <- function(x, upper) is.numeric(x) && is.null(dim(x)) &&
    all(is.finite(x)) && all(x == floor(x)) && all(x >= 1 & x <= upper)
  if (!index_ok(cell_map$interval, tt) || !index_ok(cell_map$distance_bin, jd) ||
      anyDuplicated((cell_map$distance_bin - 1) * tt + cell_map$interval))
    fail("cell_map needs every source distance/interval combination exactly once.")
  map <- cell_map[match(colnames(histories), cell_map$column), , drop = FALSE]
  valid <- is.finite(histories) & (histories == 0 | histories == 1)
  positive <- is.finite(histories) & histories == 1
  if (any(!surveyed[event] & rowSums(positive) > 0))
    fail("Positive records contradict unsurveyed event status.", which(!surveyed[event] & rowSums(positive) > 0))
  findings <- list()
  add <- function(rows, code, columns = rep(NA_integer_, length(rows)), band = rep(NA_integer_, length(rows))) {
    if (!length(rows)) return(invisible(NULL))
    known <- !is.na(columns)
    value <- rep(NA_real_, length(rows))
    value[known] <- histories[cbind(rows[known], columns[known])]
    findings[[length(findings) + 1L]] <<- data.frame(.history_row = as.integer(rows),
      .cell = as.integer(columns), .column = colnames(histories)[columns],
      interval = ifelse(known, map$interval[columns], band),
      distance_bin = map$distance_bin[columns], code = code, value = value)
    invisible(NULL)
  }
  bad_cells <- which(!valid, arr.ind = TRUE)
  add(bad_cells[, 1L], "invalid_mark", bad_cells[, 2L])
  nh <- nrow(histories)
  candidate <- first_interval <- first_bin <- rep(NA_integer_, nh)
  prior_bad <- rep(FALSE, nh)
  complete <- rep(TRUE, nh)
  for (band in seq_len(tt)) {
    selected <- which(map$interval == band)
    hits <- rowSums(positive[, selected, drop = FALSE])
    bad <- rowSums(!valid[, selected, drop = FALSE]) > 0L
    multiple <- hits > 1L
    complete <- complete & !bad & !multiple
    add(which(multiple), "multiple_distance_bins", band = rep(band, sum(multiple)))
    new <- is.na(candidate) & hits > 0L
    candidate[new] <- band
    resolved <- which(new & !prior_bad & !bad & !multiple)
    if (length(resolved)) {
      first_interval[resolved] <- band
      position <- max.col(positive[resolved, selected, drop = FALSE] * 1, ties.method = "first")
      first_bin[resolved] <- map$distance_bin[selected[position]]
    }
    prior_bad <- prior_bad | bad | multiple
  }
  unresolved <- is.na(first_interval) & (!is.na(candidate) | !complete)
  primary <- !is.na(first_interval) & first_bin <= jj & surveyed[event]
  outside <- !is.na(first_interval) & first_bin > jj
  status <- rep("no_record", nh)
  status[unresolved] <- "unresolved_first"
  status[outside] <- "first_outside_radius"
  status[primary & complete] <- "retained_first"
  status[primary & !complete] <- "retained_first_incomplete_history"
  status[!surveyed[event]] <- "unsurveyed_history"
  add(which(unresolved), "unresolved_first")
  add(which(outside), "first_outside_radius")
  add(which(!surveyed[event]), "unsurveyed_history")
  first <- data.frame(.history_row = seq_len(nh), .event = event,
    candidate_interval = candidate, first_interval = first_interval,
    first_distance_bin = first_bin, primary_eligible = primary,
    history_complete = complete, complete_history_eligible = primary & complete,
    status = status)
  canonical <- order(cell_map$distance_bin, cell_map$interval)
  native_cells <- cell_map[canonical[cell_map$distance_bin[canonical] <= jj], , drop = FALSE]
  position <- (event[primary] - 1) * (jj * tt) +
    (first_bin[primary] - 1) * tt + first_interval[primary]
  paired <- matrix(tabulate(position, nbins = m * jj * tt), m, jj * tt, byrow = TRUE,
    dimnames = list(NULL, native_cells$column))
  unknown_events <- tabulate(event[unresolved], nbins = m)
  response_complete <- surveyed & unknown_events == 0L
  paired[!response_complete, ] <- NA_integer_
  yd <- matrix(0, m, jj)
  yr <- matrix(0, m, tt)
  for (j in seq_len(jj)) yd[, j] <- rowSums(paired[, (j - 1L) * tt + seq_len(tt), drop = FALSE])
  for (band in seq_len(tt)) yr[, band] <- rowSums(paired[, seq.int(band, jj * tt, by = tt), drop = FALSE])
  events <- data.frame(.event = seq_len(m), surveyed = surveyed, response_complete = response_complete,
    known_first_records = tabulate(event[primary], nbins = m), unresolved_first_histories = unknown_events)
  audit <- NULL
  if (!is.null(native_frame)) {
    if (!methods::is(native_frame, "unmarkedFrameGDR") || native_frame@numPrimary != 1L)
      fail("native_frame must be a single-primary unmarkedFrameGDR.")
    audit <- validate_unmarked_data(native_frame, event_keys, event_id,
      interval_units = interval_units, effort = ifelse(surveyed, 1, 0), expected_event_keys = expected_event_keys)
    if (!isTRUE(audit$valid)) fail("The original native survey audit has errors; retain and inspect the source response.")
    equal <- function(a, b) isTRUE(all.equal(unname(a), unname(b), tolerance = 0, check.attributes = FALSE))
    native_meters <- native_frame@dist.breaks * if (native_frame@unitsIn == "m") 1 else 1000
    source_meters <- dist_breaks[seq_len(jj + 1L)] * if (units_in == "m") 1 else 1000
    if (!equal(native_frame@yDistance, yd) || !equal(native_frame@yRemoval, yr) ||
        !equal(native_meters, source_meters) || !equal(native_frame@period.lengths, interval_lengths))
      fail("Prepared paired first cells do not exactly match the native margins and survey protocol.")
  }
  empty <- data.frame(.history_row = integer(), .cell = integer(), .column = character(),
    interval = integer(), distance_bin = integer(), code = character(), value = numeric())
  findings <- if (length(findings)) do.call(rbind, findings) else empty
  rownames(findings) <- NULL
  protocol <- list(source_distance_breaks = dist_breaks, distance_breaks = dist_breaks[seq_len(jj + 1L)],
    radius = radius, distance_units = units_in, interval_lengths = interval_lengths, interval_units = interval_units,
    numPrimary = 1L, first_cell_order = "distance_major_then_interval",
    response_conventions = c(first_cells = "paired_first_distance_time", native = "native_two_margins",
      complete_history = "linked_record_histories_of_primary_cohort"),
    qualification = "source_response_preparation_only")
  structure(list(histories = histories, history_keys = history_keys, event_keys = event_keys,
    cell_map = cell_map, native_cells = native_cells, first = first, events = events, findings = findings,
    first_cells = paired, yDistance = yd, yRemoval = yr, protocol = protocol,
    surveyed = surveyed, survey_audit = audit, fingerprint = .msr_gdr_hash(list(histories = histories,
      history_keys = history_keys, event_keys = event_keys, cell_map = cell_map,
      surveyed = surveyed, protocol = protocol)),
    summary = list(histories = nh, events = m, primary_records = sum(primary),
      complete_histories = sum(primary & complete), invalid_cells = nrow(bad_cells),
      unresolved_events = sum(surveyed & unknown_events > 0), outside_first = sum(outside)),
    fit_eligibility_established = FALSE, density_identification_established = FALSE), class = "msr_gdr_histories")
}

#' @rdname prepare_gdr_histories
#' @export
print.msr_gdr_histories <- function(x, ...) {
  z <- x$summary
  cat("Linked GDR histories:", z$histories, "histories in", z$events, "events.\n")
  cat("Valid in-radius first records:", z$primary_records, "; complete primary histories:", z$complete_histories, ".\n")
  cat("Outside-first records:", z$outside_first, "; unresolved event responses:", z$unresolved_events, ".\n")
  cat("Original records and exclusions remain available for review.\n")
  invisible(x)
}
