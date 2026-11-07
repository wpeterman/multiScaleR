#' Prepare keyed native GDR interval covariates
#'
#' Add observation-interval predictors to a copy of an
#' `unmarkedFrameGDR`, retaining an explicit record of their alignment.
#'
#' @param object An `unmarkedFrameGDR` with paired histories and its original
#'   site, primary-period and observation covariates.
#' @param event_keys A data frame with one row per native site row and primary
#'   period, in native site-major order. Supply original survey identifiers.
#' @param event_id Character vector naming the columns of `event_keys` that
#'   uniquely identify each survey. Compound keys are supported.
#' @param interval_data A data frame containing the `event_id` columns, the
#'   interval index, and new observation predictors. It must contain exactly
#'   one row for every survey and removal interval. Rows may be shuffled;
#'   `alignment` records their explicit mapping to native order.
#' @param interval_units Nonempty time-unit label for the original removal
#'   durations, such as `"minutes"`. No time conversion is performed.
#' @param interval One character string naming the numeric interval-index
#'   column in `interval_data`. Its values must be integers from one through
#'   the number of removal intervals per primary period. Default: `"interval"`.
#'   This column is retained as an observation predictor.
#' @param site_id Optional physical-site columns in `event_keys`, passed to
#'   [validate_unmarked_data()]. Native frame rows need not be physical sites.
#' @param effort Optional event-order effort vector, passed to
#'   [validate_unmarked_data()]. Unknown, unsurveyed and surveyed-zero events
#'   retain their separate meanings.
#' @param expected_event_keys Optional independently supplied event mapping in
#'   native frame order, passed to [validate_unmarked_data()].
#' @return A list of class `msr_gdr_intervals` containing `data`, the updated
#'   native frame; `survey_audit`, the original survey audit; `alignment`, the
#'   native event/row/primary/interval positions and matching input row; and
#'   `protocol`, the units, durations, added columns and provider metadata.
#'   `fit_eligibility_established` is always `FALSE`. Structural errors stop
#'   without changing the input; interval errors carry the survey audit.
#'   Printing reports dimensions and whether the input order changed, returning
#'   the object invisibly. Full alignment and warning records remain available.
#' @details
#' Native GDR observation rows are ordered by frame row, primary period, then
#'   removal interval. The native observation-covariate setter does not verify
#'   this layout. This helper checks keys before using that setter. It rejects
#'   unknown, duplicate, missing or extra survey/interval combinations and any
#'   attempt to overwrite an existing observation column. Event identifier
#'   columns are used for alignment; all other supplied columns, including the
#'   interval index, are added to native observation covariates.
#'
#' The response matrices, existing covariates, distance geometry, durations,
#'   units, column classes, factor levels and contrasts remain unchanged.
#'   Missing predictor values are retained, and native model fitting still
#'   applies its own missingness rules. No records are aggregated or imputed.
#'   `.input_row`, `.input_row_name` and `.interval` are reserved alignment names
#'   in `event_keys`. Shared site/annual predictor names are recorded for
#'   review; scale-fitting guards against shadowed optimized site predictors
#'   remain in force. Use distinct names for fields with different meanings.
#'
#' Use native formulas such as `~ I(interval == 1)` or `~ factor(interval)`
#'   when the observation process justifies them. With equal unit intervals,
#'   hazards pi in the first interval and b thereafter reproduce the first-time
#'   law of a simple encounter Markov process regardless of persistence after
#'   an encounter. First times alone do not estimate that persistence or
#'   establish a finite corrected density. Unequal native durations retain
#'   their per-unit probability interpretation. Keep linked complete histories
#'   when assessing re-detections. This helper adds no likelihood, parameter,
#'   prior, uncertainty method or conditional-random prediction admission.
#' @examples
#' frame <- unmarked::unmarkedFrameGDR(
#'   yDistance = matrix(c(1, 0, 0, 1), 2, byrow = TRUE),
#'   yRemoval = matrix(c(1, 0, 0, 0, 1, 0), 2, byrow = TRUE),
#'   siteCovs = data.frame(site = factor(c("A", "B"))),
#'   dist.breaks = c(0, 25, 75), unitsIn = "m",
#'   period.lengths = c(1, 2, 1))
#' keys <- data.frame(site = c("A", "B"), survey = c("s1", "s2"))
#' bands <- data.frame(survey = rep(c("s2", "s1"), each = 3),
#'   interval = rep(c(3, 1, 2), 2),
#'   observer = factor(c("Ann", "Ann", "Bob", "Bob", "Ann", "Ann")))
#' prepared <- prepare_gdr_intervals(frame, keys, "survey", bands,
#'   interval_units = "minutes", site_id = "site")
#' unmarked::obsCovs(prepared$data)
#' prepared$alignment
#' @export
prepare_gdr_intervals <- function(object, event_keys, event_id, interval_data,
    interval_units, interval = "interval", site_id = NULL, effort = NULL,
    expected_event_keys = NULL) {
  survey <- validate_unmarked_data(object, event_keys, event_id, site_id,
    interval_units, effort, expected_event_keys)
  fail <- function(message, input_rows = integer()) {
    stop(structure(list(message = message, call = NULL, survey_audit = survey,
      input_rows = input_rows), class = c("msr_gdr_interval_error", "error", "condition")))
  }
  if (!isTRUE(survey$valid))
    fail("The original GDR survey audit has errors; inspect survey_audit$findings.")
  validate_character_scalar(interval, "interval")
  if (interval %in% event_id)
    fail("The interval index must have a different name from the event identifiers.")
  reserved <- c(".input_row", ".input_row_name", ".interval")
  if (any(names(event_keys) %in% reserved)) fail("event_keys uses reserved interval alignment names.")
  if (!is.data.frame(interval_data) || anyNA(names(interval_data)) ||
      anyDuplicated(names(interval_data)) || any(!nzchar(names(interval_data))) ||
      any(!c(event_id, interval) %in% names(interval_data)))
    fail("interval_data must contain uniquely named event identifiers and the interval index.")
  if (any(vapply(interval_data[event_id], function(x) !is.atomic(x) || !is.null(dim(x)), logical(1))))
    fail("Interval event identifiers must be atomic columns without matrix dimensions.")
  ne <- nrow(event_keys)
  tt <- object@numPrimary
  jr <- ncol(object@yRemoval) / tt
  total <- as.double(ne) * jr
  if (!is.finite(total) || total > .Machine$integer.max || nrow(interval_data) != total)
    fail("Supply exactly one row per survey event and native removal interval.")
  if (unmarked::obsNum(object) != ncol(object@yRemoval))
    fail("The frame observation mapping does not match its removal columns.")
  band <- interval_data[[interval]]
  bad <- if (!is.numeric(band) || !is.null(dim(band))) seq_len(nrow(interval_data)) else
    which(!is.finite(band) | band != floor(band) | band < 1 | band > jr)
  if (length(bad)) fail("Interval indices must be complete integers in the native removal range.", bad)
  text_keys <- lapply(interval_data[event_id], as.character)
  bad <- Reduce(`|`, lapply(text_keys, function(x) is.na(x) | !nzchar(trimws(x))))
  if (any(bad)) fail("Interval event identifiers must be complete.", which(bad))
  # Length prefixes preserve compound identifiers even when values contain
  # separators or non-ASCII text. Key conversion matches the survey audit.
  encode <- function(data) do.call(paste0, lapply(data, function(x) {
    x <- enc2utf8(as.character(x))
    paste0(nchar(x, type = "bytes"), ":", x)
  }))
  event <- match(encode(interval_data[event_id]), encode(event_keys[event_id]))
  if (anyNA(event)) fail("Interval rows contain unknown survey event identifiers.", which(is.na(event)))
  position <- (as.double(event) - 1) * jr + band
  duplicate <- duplicated(position) | duplicated(position, fromLast = TRUE)
  if (any(duplicate)) fail("Duplicate survey/interval combinations cannot be aligned.", which(duplicate))
  ordering <- match(seq_len(total), position)
  if (anyNA(ordering)) fail("Some required survey/interval combinations are missing.")
  added <- setdiff(names(interval_data), event_id)
  previous <- unmarked::obsCovs(object)
  overlap <- intersect(added, names(previous))
  if (length(overlap)) fail(paste0("Existing observation covariates cannot be overwritten: ",
    paste(overlap, collapse = ", "), "."))
  extra <- as.data.frame(interval_data[ordering, added, drop = FALSE])
  rownames(extra) <- if (is.null(previous)) as.character(seq_len(total)) else rownames(previous)
  covariates <- if (is.null(previous)) extra else data.frame(previous, extra, check.names = FALSE)
  out <- object
  unmarked::obsCovs(out) <- covariates
  event_order <- rep(seq_len(ne), each = jr)
  native <- data.frame(.event = event_order,
    .site_row = (event_order - 1L) %/% tt + 1L,
    .primary = (event_order - 1L) %% tt + 1L,
    .interval = rep(seq_len(jr), ne), .input_row = ordering,
    .input_row_name = rownames(interval_data)[ordering])
  alignment <- data.frame(native, event_keys[event_order, , drop = FALSE], check.names = FALSE)
  rownames(alignment) <- NULL
  shared <- intersect(added, union(names(unmarked::siteCovs(object)),
    names(unmarked::yearlySiteCovs(object))))
  structure(list(data = out, survey_audit = survey, alignment = alignment,
    protocol = c(survey$protocol, list(native_rows = nrow(object@yDistance),
      primary_periods = tt, events = ne, removal_intervals = jr,
      interval_column = interval, added_columns = added,
      shared_site_or_annual_names = shared, reordered = !identical(ordering, seq_len(total)),
      unmarked_version = as.character(utils::packageVersion("unmarked")),
      qualification = "native_data_alignment_only")),
    fit_eligibility_established = FALSE), class = "msr_gdr_intervals")
}

#' @rdname prepare_gdr_intervals
#' @param x A prepared `msr_gdr_intervals` object to print.
#' @param ... Unused.
#' @export
print.msr_gdr_intervals <- function(x, ...) {
  p <- x$protocol
  cat("Prepared GDR interval covariates:", p$native_rows, "native rows,",
    p$primary_periods, "primary periods,", p$removal_intervals, "intervals.\n")
  cat("Added", length(p$added_columns), "observation columns.\n")
  cat(if (isTRUE(p$reordered)) "Input order changed; inspect alignment.\n" else
    "Input already follows native interval order.\n")
  if (x$survey_audit$summary$warnings)
    cat("Survey audit warnings:", x$survey_audit$summary$warnings, "; inspect survey_audit.\n")
  invisible(x)
}
