#' Audit Survey Identities and Paired Distance-Removal Counts
#'
#' Check an `unmarkedFrameGDR` against a table of original survey identifiers
#' before you fit a model. The audit looks for split or duplicated surveys,
#' covariates out of line with the surveys, distance and removal counts that
#' disagree, and counts that conflict with survey effort. It only reports. It
#' does not combine surveys, reorder rows, remove records or repair counts.
#'
#' @param object An `unmarkedFrameGDR` holding paired distance and removal
#'   counts for point surveys. Other frame classes are not supported.
#' @param event_keys A data frame with one row per survey event, in the
#'   frame's order. An event is one primary period at one frame row, so with
#'   one visit per point there is one row per frame row. With several primary
#'   periods, list all periods of the first frame row, then all periods of the
#'   second, and so on. Use the original survey identifiers, such as a point
#'   name and date, rather than a row number.
#' @param event_id Character vector naming the columns of `event_keys` that
#'   together identify one survey event. Their combinations should be unique;
#'   repeated combinations are reported.
#' @param site_id Optional character vector naming the columns of
#'   `event_keys` that identify the physical site. A site may appear in
#'   several events (repeat visits). The default `NULL` skips the site checks.
#' @param interval_units Character string naming the time unit of
#'   `object@period.lengths`, for example `"minutes"`. Removal probabilities
#'   are per unit of this time. The label is recorded, not converted.
#' @param effort Optional numeric vector with one value per event, in event
#'   order, such as minutes surveyed. A positive value marks a completed
#'   survey, zero marks an event that was not surveyed (its counts must be
#'   `NA`), and `NA` marks unknown effort. With the default `NULL`, the audit
#'   adds a warning, because counts alone cannot separate a completed survey
#'   with no detections from one that never happened.
#' @param expected_event_keys Optional data frame, built independently of
#'   `event_keys` (for example from the field schedule), with the `event_id`
#'   columns in the order the frame should have. Rows whose keys differ are
#'   reported. A key table cannot show by itself that it is in the right
#'   order.
#' @param x An object returned by `validate_unmarked_data()`.
#' @param ... Further arguments passed to the data-frame print method.
#' @return A list of class `msr_unmarked_audit`:
#'   * `valid`: `TRUE` when no errors were found. Warnings can remain, so read
#'     them too.
#'   * `findings`: a data frame with one row per finding. Columns are
#'     `severity` (`"error"` or `"warning"`), `code` (see "Finding codes"),
#'     `message`, the event's position in frame order (`.event`), its frame
#'     row (`.site_row`) and primary period (`.primary`), followed by the
#'     `event_keys` columns for that event. Findings about the whole frame
#'     have missing positions.
#'   * `events`: `event_keys` as supplied.
#'   * `summary`: the number of frame rows (`sites`), `primary_periods`,
#'     `events`, `errors` and `warnings`, and `missing_histories`, the number
#'     of events with both count tables fully missing.
#'   * `protocol`: `distance_units`, `interval_units`, `distance_breaks` and
#'     `interval_lengths` from the frame.
#'   * `effort`: the supplied effort vector.
#'
#'   Printing shows the error and warning counts and at most ten findings,
#'   and returns the object invisibly. The full table stays in `$findings`.
#' @details
#' # Identify survey events and physical sites
#'
#' A physical site can be visited more than once. Give each visit its own
#' event identifier, and keep the site in a separate column named by
#' `site_id`. A row number from a spreadsheet cannot show that a detection
#' record and a zero-detection record belong to the same survey, so resolve
#' survey identities when you build the frame from the original records, for
#' example with [prepare_gdr_histories()].
#'
#' The distance and removal tables count the same birds, so their totals
#' must match for each event. Matching totals are necessary but cannot reveal
#' a survey split across two frame rows. Supplying the original event keys
#' lets the audit report such a split as `duplicate_event`. Key names that
#' look like row indices (such as `index`, `row_id` or `X`) give a warning.
#' When an `event_id` or `site_id` column also appears in the site, yearly or
#' observation covariates, its values are compared with `event_keys` row by
#' row. Without `expected_event_keys` or such shared columns, unique keys
#' alone cannot prove that the rows are aligned.
#'
#' # Record survey effort
#'
#' A survey that took place and detected nothing is an observed zero. A
#' survey that never took place should have fully missing counts in both
#' tables. Supply `effort` so the audit can tell these apart. Counts that are
#' missing in only one table, or in only some columns, are reported because
#' their totals cannot be compared. Character and factor values, including the
#' literal string `"None"`, are kept as values and not treated as missing.
#'
#' # Finding codes
#'
#' | Code | Severity | Meaning |
#' |:---|:---|:---|
#' | `frame_layout` | error | The count tables do not have matching rows and complete primary-period blocks. Other checks are skipped |
#' | `survey_type` | error | The frame is not a point survey |
#' | `event_layout` | error | `event_keys` does not have one row per event. Other checks are skipped |
#' | `missing_event_key` | error | An event identifier is missing or blank |
#' | `duplicate_event` | error | One event key appears in more than one event, as when a survey is split |
#' | `index_like_event_key` | warning | An `event_id` column name looks like a row index |
#' | `missing_site_key` | error | A site identifier is missing or blank |
#' | `physical_site_alignment` | error | Primary periods in one frame row have different site keys |
#' | `event_alignment` | error | The key differs from `expected_event_keys` in that row |
#' | `duplicate_reference_event` | error | `expected_event_keys` repeats an event |
#' | `covariate_layout` | error | A covariate table has the wrong number of rows or repeated names |
#' | `covariate_key_alignment` | error | A key column in a covariate table differs from `event_keys` |
#' | `response_slot_alignment` | error | The frame's response and removal tables disagree |
#' | `distance_breaks` | error | Distance breaks do not start at zero, increase, and match the distance columns |
#' | `distance_units` | error | Distance units are not `"m"` or `"km"` |
#' | `interval_lengths` | error | Interval lengths are not positive whole numbers, one per removal column |
#' | `invalid_counts` | error | A count is negative, fractional, infinite or `NaN` |
#' | `partial_history` | error | An event's counts are only partly missing |
#' | `unequal_totals` | error | Distance and removal totals differ for the event |
#' | `effort_unspecified` | warning | No `effort` was supplied |
#' | `invalid_effort` | error | Effort is negative or infinite |
#' | `unsurveyed_counts` | error | An event with zero effort has counts, including zeros |
#' | `surveyed_missing` | error | An event with positive effort has fully missing counts |
#' | `unknown_effort` | warning | Effort is `NA` for the event |
#'
#' The audit does not decide how to treat a bird first detected beyond the
#' survey radius and later inside it, or how to turn raw bird records into
#' counts. [prepare_gdr_histories()] handles that step. See the "Check that
#' the records line up with the surveys" section of
#' `vignette("v12-distance-removal", package = "multiScaleR")` for the audit
#' in a full workflow.
#' @seealso [prepare_gdr_histories()], [prepare_gdr_intervals()],
#'   [expand_gdr_intervals()], [sim_dat_gdr()]
#' @examples
#' # Two visits to site A, one frame row each, with matching totals
#' frame <- unmarked::unmarkedFrameGDR(
#'   yDistance = matrix(c(1, 0, 0, 0), 2, byrow = TRUE),
#'   yRemoval = matrix(c(0, 1, 0, 0), 2, byrow = TRUE),
#'   dist.breaks = c(0, 25, 75), unitsIn = "m",
#'   period.lengths = c(1, 2))
#' keys <- data.frame(site = c("A", "A"), visit = c("spring", "summer"))
#' audit <- validate_unmarked_data(frame, keys, c("site", "visit"),
#'   site_id = "site", interval_units = "minutes", effort = c(3, 3))
#' audit$valid      # TRUE: no errors and, with effort supplied, no warnings
#' audit
#'
#' # Give both rows the same key, as if one survey had been split in two.
#' # Their totals still match, but the event mapping is now invalid.
#' split_keys <- keys
#' split_keys$visit[2] <- "spring"
#' split_audit <- validate_unmarked_data(frame, split_keys, c("site", "visit"),
#'   site_id = "site", interval_units = "minutes", effort = c(3, 3))
#' split_audit$valid
#' split_audit$findings[, c("severity", "code", ".event", "visit")]
#' @export
validate_unmarked_data <- function(object, event_keys, event_id, site_id = NULL,
                                   interval_units, effort = NULL,
                                   expected_event_keys = NULL) {
  if (!methods::is(object, "unmarkedFrameGDR"))
    stop("Only unmarkedFrameGDR survey audits are currently supported.", call. = FALSE)
  if (!is.data.frame(event_keys) || !ncol(event_keys) ||
      anyDuplicated(names(event_keys)) || any(!nzchar(names(event_keys))) ||
      any(vapply(event_keys, function(x) !is.atomic(x) || !is.null(dim(x)), logical(1))))
    stop("event_keys must be a data frame with uniquely named atomic columns.", call. = FALSE)
  check_columns <- function(x, label, optional = FALSE) {
    if (optional && is.null(x)) return(invisible(NULL))
    if (!is.character(x) || !length(x) || anyNA(x) || anyDuplicated(x) ||
        any(!x %in% names(event_keys)))
      stop(label, " must name distinct columns of event_keys.", call. = FALSE)
  }
  check_columns(event_id, "event_id")
  check_columns(site_id, "site_id", TRUE)
  if (!is.character(interval_units) || length(interval_units) != 1L ||
      is.na(interval_units) || !nzchar(trimws(interval_units)))
    stop("Supply an explicit, nonempty interval_units label.", call. = FALSE)
  reserved <- c("severity", "code", "message", ".event", ".site_row", ".primary")
  if (any(names(event_keys) %in% reserved))
    stop("event_keys uses reserved audit column names.", call. = FALSE)
  findings <- list()
  yd <- object@yDistance; yr <- object@yRemoval
  m <- nrow(yd); tt <- object@numPrimary
  good_layout <- is.numeric(yd) && is.numeric(yr) &&
    length(tt) == 1L && is.finite(tt) && tt >= 1 && tt == round(tt) &&
    m > 0L && nrow(yr) == m && ncol(yd) > 0L && ncol(yr) > 0L &&
    ncol(yd) %% tt == 0L && ncol(yr) %% tt == 0L
  add <- function(code, message, rows = NA_integer_, severity = "error") {
    if (!length(rows)) return(invisible(NULL))
    rows <- as.integer(rows)
    keys <- event_keys[rep(NA_integer_, length(rows)), , drop = FALSE]
    available <- !is.na(rows) & rows >= 1L & rows <= nrow(event_keys)
    if (any(available)) keys[available, ] <- event_keys[rows[available], , drop = FALSE]
    findings[[length(findings) + 1L]] <<- cbind(data.frame(
      severity = severity, code = code, message = message, .event = rows,
      .site_row = if (good_layout) (rows - 1L) %/% tt + 1L else NA_integer_,
      .primary = if (good_layout) (rows - 1L) %% tt + 1L else NA_integer_,
      stringsAsFactors = FALSE), keys)
    invisible(NULL)
  }
  finish <- function(missing_histories = NA_integer_) {
    empty <- data.frame(severity = character(), code = character(), message = character(),
      .event = integer(), .site_row = integer(), .primary = integer())
    tab <- if (length(findings)) do.call(rbind, findings) else
      cbind(empty, event_keys[FALSE, , drop = FALSE])
    rownames(tab) <- NULL
    errors <- sum(tab$severity == "error"); warnings <- sum(tab$severity == "warning")
    structure(list(valid = errors == 0L, findings = tab, events = event_keys,
      summary = list(sites = m, primary_periods = tt,
        events = if (good_layout) m * tt else NA_integer_, errors = errors,
        warnings = warnings, missing_histories = missing_histories),
      protocol = list(distance_units = object@unitsIn, interval_units = interval_units,
        distance_breaks = object@dist.breaks, interval_lengths = object@period.lengths),
      effort = effort), class = "msr_unmarked_audit")
  }
  if (!good_layout) {
    add("frame_layout", "Distance/removal matrices need matching sites and complete primary-period blocks.")
    return(finish())
  }
  ne <- m * tt; jd <- ncol(yd) / tt; jr <- ncol(yr) / tt
  if (!identical(object@survey, "point"))
    add("survey_type", "Native GDR validation currently requires point surveys.")
  if (nrow(event_keys) != ne) {
    add("event_layout", "event_keys must have one row per site and primary period.")
    return(finish())
  }
  key_text <- function(x, columns) {
    as.data.frame(lapply(x[columns], as.character), stringsAsFactors = FALSE)
  }
  keys <- key_text(event_keys, event_id)
  bad <- Reduce(`|`, lapply(keys, function(x) is.na(x) | !nzchar(trimws(x))))
  add("missing_event_key", "Survey event identifiers must be complete.", which(bad))
  dup <- duplicated(keys) | duplicated(keys, fromLast = TRUE)
  add("duplicate_event", "Original event keys occur in multiple frame events; inspect split surveys.", which(dup))
  if (any(grepl("(^|[._])index$|^index$|^row([._]?id|num|number)?$|^X$|^Unnamed", event_id,
                ignore.case = TRUE)))
    add("index_like_event_key", "An index-like event key cannot establish original survey identity.",
        severity = "warning")
  if (!is.null(site_id)) {
    sites <- key_text(event_keys, site_id)
    bad_site <- Reduce(`|`, lapply(sites, function(x) is.na(x) | !nzchar(trimws(x))))
    add("missing_site_key", "Physical-site identifiers must be complete.", which(bad_site))
    # Primary periods within one native row share that row's physical site.
    if (tt > 1L) {
      first <- rep(seq.int(1L, ne, by = tt), each = tt)
      changed <- Reduce(`|`, lapply(sites, function(x)
        !is.na(x) & !is.na(x[first]) & x != x[first]))
      bad_rows <- unique((which(changed) - 1L) %/% tt)
      block <- as.vector(t(outer(bad_rows * tt, seq_len(tt), `+`)))
      add("physical_site_alignment", "Primary periods in one native site row have different physical-site keys.", block)
    }
  }
  if (!is.null(expected_event_keys)) {
    if (!is.data.frame(expected_event_keys) || any(!event_id %in% names(expected_event_keys)) ||
        nrow(expected_event_keys) != ne || anyDuplicated(names(expected_event_keys)))
      stop("expected_event_keys must contain the event_id columns and one row per event.", call. = FALSE)
    reference <- key_text(expected_event_keys, event_id)
    aligned <- Reduce(`&`, Map(function(x, y) !is.na(x) & !is.na(y) & x == y, keys, reference))
    add("event_alignment", "Event keys differ from the independently supplied frame order.", which(!aligned))
    add("duplicate_reference_event", "The independent reference contains duplicate event identifiers.",
        which(duplicated(reference) | duplicated(reference, fromLast = TRUE)))
  }
  for (nm in c("siteCovs", "yearlySiteCovs", "obsCovs")) {
    covs <- methods::slot(object, nm)
    expected <- switch(nm, siteCovs = m, yearlySiteCovs = ne, obsCovs = ne * jr)
    if (is.null(covs)) next
    if (nrow(covs) != expected || anyDuplicated(names(covs))) {
      add("covariate_layout", paste(nm, "does not have the required native row layout and unique names."))
      next
    }
    shared <- intersect(names(covs), unique(c(event_id, site_id)))
    if (!length(shared)) next
    pos <- switch(nm, siteCovs = rep(seq_len(m), each = tt), yearlySiteCovs = seq_len(ne),
                  obsCovs = seq_len(ne * jr))
    event_pos <- if (nm == "obsCovs") rep(seq_len(ne), each = jr) else seq_len(ne)
    for (column in shared) {
      a <- as.character(covs[[column]][pos]); b <- as.character(event_keys[[column]][event_pos])
      matched <- !is.na(a) & !is.na(b) & a == b
      add("covariate_key_alignment", paste(nm, "column", column, "does not match event_keys."),
          unique(event_pos[!matched]))
    }
  }
  if (!identical(dim(object@y), dim(yr)) || !isTRUE(all.equal(object@y, yr, check.attributes = FALSE)))
    add("response_slot_alignment", "The inherited response and removal-count slots disagree.")
  db <- object@dist.breaks; pl <- object@period.lengths
  if (length(db) != jd + 1L || any(!is.finite(db)) || db[1L] != 0 || any(diff(db) <= 0))
    add("distance_breaks", "Distance breaks must start at zero, increase, and match the distance columns.")
  if (length(object@unitsIn) != 1L || !object@unitsIn %in% c("m", "km"))
    add("distance_units", "Native GDR distance units must be m or km.")
  if (length(pl) != jr || any(!is.finite(pl) | pl <= 0 | pl != round(pl)))
    add("interval_lengths", "Native removal intervals must have positive integer lengths in the stated unit.")
  d <- matrix(as.vector(t(yd)), ne, jd, byrow = TRUE)
  r <- matrix(as.vector(t(yr)), ne, jr, byrow = TRUE)
  for (nm in c("distance", "removal")) {
    x <- if (nm == "distance") d else r
    invalid <- is.nan(x) | (!is.na(x) & (!is.finite(x) | x < 0 | x != round(x)))
    add("invalid_counts", paste(nm, "counts must be finite nonnegative integers or missing."),
        which(rowSums(invalid) > 0))
  }
  dn <- rowSums(is.na(d)); rn <- rowSums(is.na(r))
  missing <- dn == jd & rn == jr
  complete <- dn == 0 & rn == 0
  add("partial_history", "Both paired margins must be fully observed or fully missing for an event.",
      which(!missing & !complete))
  mismatch <- complete & rowSums(d, na.rm = TRUE) != rowSums(r, na.rm = TRUE)
  add("unequal_totals", "Paired distance and removal margins describe different detected totals.", which(mismatch))
  if (is.null(effort)) {
    add("effort_unspecified", "Counts alone cannot distinguish unrecorded effort from a surveyed zero.",
        severity = "warning")
  } else {
    if (!is.numeric(effort) || length(effort) != ne)
      stop("effort must be a numeric vector with one value per event.", call. = FALSE)
    invalid <- !is.na(effort) & (!is.finite(effort) | effort < 0)
    add("invalid_effort", "Effort must be finite and nonnegative, or missing.", which(invalid))
    add("unsurveyed_counts", "Zero-effort events require missing counts, not recorded detections or zeros.",
        which(!is.na(effort) & effort == 0 & !missing))
    add("surveyed_missing", "Positive-effort events have fully missing paired histories.",
        which(!is.na(effort) & effort > 0 & missing))
    add("unknown_effort", "Survey effort is unknown for this event.", which(is.na(effort)), "warning")
  }
  finish(sum(missing))
}

#' @rdname validate_unmarked_data
#' @export
print.msr_unmarked_audit <- function(x, ...) {
  cat("GDR survey audit:", x$summary$events, "events;", x$summary$errors,
      "errors;", x$summary$warnings, "warnings\n")
  if (nrow(x$findings)) print(utils::head(x$findings, 10L), ...)
  if (nrow(x$findings) > 10L)
    cat("Showing the first 10 of", nrow(x$findings), "findings; see $findings for all records.\n")
  invisible(x)
}
