#' Prepare Linked First-Encounter and Complete-History Responses
#'
#' Retain original bird histories while deriving paired first-distance/time
#' cells and native GDR margins from the first encounter across all distances.
#' @param histories Numeric or logical matrix with one row per linked individual
#'   within an event and named source-cell columns. Valid marks are zero or one.
#'   Missing, nonfinite and other numeric marks remain in the returned matrix
#'   and are reported for review; they are never converted to zero.
#' @param history_keys Data frame in history-row order, retaining the original
#'   event and individual identifiers and any additional source metadata.
#' @param event_keys Data frame in the intended native event order, including
#'   surveyed events with no recorded individuals. Keys must be unique.
#' @param event_id Character vector naming columns that identify an event in
#'   both key tables. Compound identifiers are supported.
#' @param individual_id Character vector naming individual identifier columns
#'   in `history_keys`. The event/individual combination must be unique.
#' @param cell_map Data frame with `column`, `interval` and `distance_bin`.
#'   It must map every named history column to exactly one source distance/time
#'   cell. Interval and distance indices are positive integers. Input columns
#'   may have any order; the complete rectangular source layout is required.
#' @param dist_breaks Increasing numeric source distance boundaries, beginning
#'   at zero, in `units_in`. Only the final boundary may be infinite. Retain
#'   outside-radius bins here when they occur in the original records.
#' @param radius Positive finite survey radius in `units_in`, equal to a source
#'   distance boundary. Bins cannot be split without finer source measurements.
#' @param interval_lengths Positive integer native removal durations, one per
#'   interval, measured in `interval_units`. Unequal durations are supported.
#' @param interval_units Nonempty time-unit label, such as `"seconds"`.
#'   Values are recorded without conversion.
#' @param units_in Source distance units, `"m"` or `"km"`.
#' @param surveyed Complete logical vector in event order. Supply `TRUE` for
#'   completed surveys and `FALSE` for unsurveyed events. There is no default
#'   assumption about survey effort. A positive mark at an unsurveyed event
#'   raises an error.
#' @param expected_event_keys Optional independent event mapping in intended
#'   order. It must contain the `event_id` columns. This detects mapping drift;
#'   a supplied key table cannot prove its own alignment.
#' @param native_frame Optional existing single-primary `unmarkedFrameGDR`.
#'   Its two margins, inherited response, finite distance geometry, units and
#'   durations must match the prepared response. The original survey audit is
#'   returned. This argument validates data and never evaluates saved callbacks.
#' @return A list of class `msr_gdr_histories` with:
#'   \itemize{
#'     \item `histories`, `history_keys`, `event_keys` and `cell_map`: unchanged
#'       source values, order, identifiers and attributes.
#'     \item `native_cells`: the in-radius cell map in distance-major then
#'       interval order. Rows correspond to columns of `first_cells`.
#'     \item `first`: original `.history_row`, matched `.event`, earliest known
#'       positive `candidate_interval`, unambiguous `first_interval` and
#'       `first_distance_bin`, `primary_eligible`, `history_complete`,
#'       `complete_history_eligible`, and the explicit preparation `status`.
#'       `history_complete` requires complete binary cells and a unique distance
#'       state per interval; it is not a judgment that other source marks are
#'       biological errors.
#'       Eligibility here concerns data layout, not a fitted model.
#'     \item `events`: `.event`, `surveyed`, `response_complete`,
#'       `known_first_records` and `unresolved_first_histories` in event order.
#'       Counts of known records remain available even in quarantined events.
#'     \item `findings`: `.history_row`, `.cell`, `.column`, `interval`,
#'       `distance_bin`, `code` and original numeric `value`. History-level
#'       findings have missing cell/value fields. No source row is deleted.
#'     \item `first_cells`, `yDistance`, `yRemoval`: paired first-cell counts
#'       and their exact two margins. Unsurveyed events and events with an
#'       unresolved first encounter have fully missing response rows; completed
#'       surveys with no eligible records have observed zero rows.
#'     \item `protocol`: original and finite-circle geometry, distance/time
#'       units, durations, one-primary layout, response conventions and
#'       preparation scope. No biological target population is inferred.
#'     \item `surveyed`: the original effort-status vector; `survey_audit`:
#'       optional native audit; `fingerprint`: SHA256 of source inputs and
#'       protocol; `summary`: history/event, retained and unresolved counts.
#'     \item `fit_eligibility_established` and
#'       `density_identification_established`: always `FALSE`.
#'   }
#'   Printing reports these preparation counts and returns the object invisibly.
#'   Structural errors carry the implicated `history_rows` when recoverable.
#' @details
#' Determine first encounters before truncating distance. A bird first recorded
#'   outside the radius remains outside the primary response even if a later
#'   mark is near the point. Multiple positive distance cells in one interval
#'   do not determine a unique state for this binned layout. They are retained
#'   for adjudication and never collapsed into a single distance. Such marks
#'   may require a different movement or recording model. If earlier or
#'   first-band marks make the first encounter ambiguous, the event's aggregate
#'   response is missing pending adjudication, rather than a partial zero count.
#'
#' A uniquely known in-radius first encounter is retained when only later marks
#'   are malformed. That history is excluded from the complete-history subset.
#'   Complete histories of the primary cohort and native first-time margins are
#'   different responses. Their likelihood constants and conditioning differ.
#'   Keep re-detections and original outside-radius cells for separate analyses;
#'   within-count bands are not independent primary visits.
#'
#' No source correction, probability model, density correction, parameter,
#'   prior, model-selection admission or uncertainty procedure is added. Empty
#'   linked histories do not establish the size of an unobserved population.
#'   Native geometry comparisons convert meters/kilometers to meters only in
#'   the comparison copy; inputs and native objects remain unchanged.
#' @seealso [validate_unmarked_data()], [prepare_gdr_intervals()],
#'   [audit_gdr_integration()]
#' @examples
#' map <- expand.grid(interval = 1:2, distance_bin = 1:3)
#' map$column <- paste0("cell", 1:6)
#' h <- matrix(0, 2, 6, dimnames = list(NULL, map$column))
#' h[1, c(5, 2)] <- 1 # Outside first, then near: excluded from the circle.
#' h[2, 1] <- 1
#' h[2, 2] <- 2       # Malformed later mark: the valid first still counts.
#' prepared <- prepare_gdr_histories(h,
#'   history_keys = data.frame(survey = c("a", "a"), bird = c("one", "two")),
#'   event_keys = data.frame(survey = c("a", "b")), event_id = "survey",
#'   individual_id = "bird", cell_map = map, dist_breaks = c(0, 25, 75, Inf),
#'   radius = 75, interval_lengths = c(75, 75), interval_units = "seconds",
#'   units_in = "m", surveyed = c(TRUE, FALSE))
#' prepared$first
#' prepared$yRemoval
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
#' @param x Prepared `msr_gdr_histories` data to print.
#' @param ... Unused.
#' @export
print.msr_gdr_histories <- function(x, ...) {
  z <- x$summary
  cat("Linked GDR histories:", z$histories, "histories in", z$events, "events.\n")
  cat("Valid in-radius first records:", z$primary_records, "; complete primary histories:", z$complete_histories, ".\n")
  cat("Outside-first records:", z$outside_first, "; unresolved event responses:", z$unresolved_events, ".\n")
  cat("Original records and exclusions remain available for review.\n")
  invisible(x)
}
