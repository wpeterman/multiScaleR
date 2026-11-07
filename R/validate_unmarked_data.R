#' Audit Survey Identities and Paired Distance-Removal Counts
#'
#' Check an explicit survey-event mapping and the paired counts in an
#' `unmarkedFrameGDR`. Findings identify records that need review. The function
#' does not combine surveys, reorder rows, remove records, or repair counts.
#'
#' @param object An `unmarkedFrameGDR` containing paired distance and removal
#'   counts. Other frame classes are not yet supported.
#' @param event_keys A data frame with one row per site and primary period,
#'   in native site-major order: all periods at the first site, then all periods
#'   at the second site. Supply original survey identifiers, independently of
#'   any artificial frame-row index.
#' @param event_id Character vector naming the columns of `event_keys` that
#'   together identify one survey event. Their combinations must be unique.
#' @param site_id Optional character vector naming physical-site columns in
#'   `event_keys`. Repeated sites across distinct events are allowed.
#' @param interval_units One nonempty character string giving the time unit of
#'   `object@period.lengths`, for example `"minutes"`. Removal probabilities
#'   are defined per this unit. The audit records the label without converting it.
#' @param effort Optional numeric vector in event order. Positive values mean
#'   surveyed events; zero means not surveyed; `NA` means effort is unknown.
#'   Unsurveyed events must contain missing counts rather than observed zeros.
#'   Without this vector, the audit cannot establish effort from counts alone.
#' @param expected_event_keys Optional independent data frame in the intended
#'   frame order, containing the `event_id` columns. Use it to detect a reordered
#'   or unmatched event mapping. Keys alone cannot prove their own alignment.
#' @return A list with class `msr_unmarked_audit`. `valid` is `TRUE` when no
#'   errors were found; warnings still require interpretation. `findings` is a
#'   data frame with severity, code, message, native event/row/primary positions,
#'   and the supplied key columns. Global findings have missing positions.
#'   `events` retains the supplied event mapping. `summary` records frame
#'   dimensions, error/warning counts, and the number of missing event histories.
#'   `protocol` records distance and removal units, breaks, and interval lengths.
#'   `effort` retains the supplied effort vector.
#'   Printing shows the error/warning counts and at most ten findings, returning
#'   the object invisibly. The complete findings table remains in the object.
#' @details
#' Distance and removal margins describe the same detected individuals. Equal
#' totals are necessary but cannot identify an event that was split into two
#' artificial surveys. Supply the original event keys so such duplicates can
#' be reported. Index-like key names generate a warning, not an inferred repair.
#'
#' Fully missing histories differ from surveyed events with zero detections.
#' Partly missing paired margins are flagged because their totals cannot be
#' compared reliably. Shared key columns in site, annual, or observation
#' covariates are checked against the supplied mapping in native row order.
#' Character and factor categories, including the literal string `"None"`,
#' are preserved. This audit does not decide how detections outside a survey
#' radius should be truncated or how raw bird histories should be aggregated.
#' @examples
#' frame <- unmarked::unmarkedFrameGDR(
#'   yDistance = matrix(c(1, 0, 0, 0), 2, byrow = TRUE),
#'   yRemoval = matrix(c(0, 1, 0, 0), 2, byrow = TRUE),
#'   dist.breaks = c(0, 25, 75), unitsIn = "m",
#'   period.lengths = c(1, 2))
#' keys <- data.frame(site = c("A", "A"), visit = c("spring", "summer"))
#' audit <- validate_unmarked_data(frame, keys, c("site", "visit"),
#'   site_id = "site", interval_units = "minutes", effort = c(3, 3))
#' audit$valid
#' audit$findings
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
#' @param x An object returned by `validate_unmarked_data()`.
#' @param ... Further arguments passed to the data-frame print method.
#' @export
print.msr_unmarked_audit <- function(x, ...) {
  cat("GDR survey audit:", x$summary$events, "events;", x$summary$errors,
      "errors;", x$summary$warnings, "warnings\n")
  if (nrow(x$findings)) print(utils::head(x$findings, 10L), ...)
  if (nrow(x$findings) > 10L)
    cat("Showing the first 10 of", nrow(x$findings), "findings; see $findings for all records.\n")
  invisible(x)
}
