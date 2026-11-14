#' Add Interval-Level Covariates to a Distance-Removal Frame
#'
#' Add covariates that can change from one removal interval to the next,
#' such as who was listening in each interval, to a copy of an
#' `unmarkedFrameGDR`. Rows are matched by survey identifier and interval
#' number, so the input can be in any order, and the alignment is recorded.
#'
#' @param object An `unmarkedFrameGDR` with its original counts and its site,
#'   yearly and observation covariates.
#' @param event_keys A data frame with one row per survey event, in the
#'   frame's order (all primary periods of the first frame row, then the
#'   second, and so on). Use the original survey identifiers. The names
#'   `.input_row`, `.input_row_name` and `.interval` are reserved.
#' @param event_id Character vector naming the column(s) of `event_keys` that
#'   identify each survey.
#' @param interval_data A data frame with the `event_id` columns, the interval
#'   number column and the new covariates. It needs exactly one row for every
#'   survey and removal interval. Rows can be in any order, and `alignment`
#'   records how they were matched.
#' @param interval_units Character string naming the time unit of the removal
#'   intervals, such as `"minutes"`. No time conversion is performed.
#' @param interval Character string, default `"interval"`. The name of the
#'   interval number column in `interval_data`. Its values must be whole
#'   numbers from 1 to the number of removal intervals per primary period. The
#'   column is also added to the observation covariates, so formulas can use
#'   it.
#' @param site_id Optional column(s) of `event_keys` that identify the
#'   physical site, passed to [validate_unmarked_data()]. Frame rows need not
#'   be physical sites.
#' @param effort Optional numeric vector of survey effort, one value per
#'   event, passed to [validate_unmarked_data()]. Positive values mark
#'   completed surveys, zero marks unsurveyed events and `NA` marks unknown
#'   effort.
#' @param expected_event_keys Optional independent event list in frame
#'   order, passed to [validate_unmarked_data()].
#' @param x A prepared `msr_gdr_intervals` object to print.
#' @param ... Not used.
#' @return A list of class `msr_gdr_intervals`:
#'   * `data`: a copy of `object` with the new observation covariates.
#'   * `survey_audit`: the result of [validate_unmarked_data()] for the frame.
#'   * `alignment`: one row per observation row in frame order, with the
#'     event (`.event`), frame row (`.site_row`), primary period
#'     (`.primary`), interval (`.interval`), the matching row of
#'     `interval_data` (`.input_row`, `.input_row_name`) and the `event_keys`
#'     columns.
#'   * `protocol`: the survey protocol from the audit (units, distance breaks,
#'     interval lengths), plus the frame dimensions, the `interval_column`,
#'     the `added_columns`, any added names that also occur among the site or
#'     yearly covariates (`shared_site_or_annual_names`), whether the input
#'     order changed (`reordered`) and the `unmarked` version.
#'   * `fit_eligibility_established`: always `FALSE`.
#'
#'   Printing reports the dimensions and whether the input order changed,
#'   and returns the object invisibly. If the survey audit finds an error, or
#'   the interval rows cannot be aligned, the function stops with a condition
#'   of class `msr_gdr_interval_error` carrying `survey_audit` and the
#'   offending `input_rows`. The input frame is never changed.
#' @details
#' # How rows are aligned
#'
#' `unmarked` orders observation rows by frame row, primary period and then
#' removal interval, and its own setter for observation covariates does not
#' check that order. This function checks the keys first. It rejects unknown,
#' duplicate, missing or extra survey and interval combinations, and any
#' attempt to overwrite an existing observation column. The event identifier
#' columns are used only for matching. Every other column of
#' `interval_data`, including the interval number, is added.
#'
#' The counts, existing covariates, distance bands, interval lengths, units,
#' column classes, factor levels and contrasts are unchanged. Missing values
#' are kept, and `unmarked` applies its own missing-data rules when fitting.
#' Nothing is aggregated or imputed. Added names that also occur among the
#' site or yearly covariates are recorded for review. Use distinct names for
#' fields with different meanings, because [multiScale_optim()] will not
#' estimate a scale for a site covariate whose name is shadowed this way.
#'
#' # Model the first interval
#'
#' Use formulas such as `~ I(interval == 1)` or `~ factor(interval)` when the
#' observation process justifies them. With equal unit intervals, hazards pi
#' in the first interval and b thereafter reproduce the first-time law of a
#' simple encounter Markov process regardless of persistence after an
#' encounter. First times alone do not estimate that persistence or
#' establish a finite corrected density. Unequal durations keep their
#' per-unit probability interpretation. Keep linked complete histories when
#' assessing re-detections. This function only arranges data. It adds no
#' model terms, parameters or uncertainty methods.
#'
#' For covariates that are constant within a survey, [expand_gdr_intervals()]
#' builds `interval_data` for you.
#' @seealso [expand_gdr_intervals()], [validate_unmarked_data()],
#'   [prepare_gdr_histories()]
#' @examples
#' frame <- unmarked::unmarkedFrameGDR(
#'   yDistance = matrix(c(1, 0, 0, 1), 2, byrow = TRUE),
#'   yRemoval = matrix(c(1, 0, 0, 0, 1, 0), 2, byrow = TRUE),
#'   siteCovs = data.frame(site = factor(c("A", "B"))),
#'   dist.breaks = c(0, 25, 75), unitsIn = "m",
#'   period.lengths = c(1, 2, 1))
#' keys <- data.frame(site = c("A", "B"), survey = c("s1", "s2"))
#' # One row per survey and interval, deliberately out of frame order
#' bands <- data.frame(survey = rep(c("s2", "s1"), each = 3),
#'   interval = rep(c(3, 1, 2), 2),
#'   observer = factor(c("Ann", "Ann", "Bob", "Bob", "Ann", "Ann")))
#' prepared <- prepare_gdr_intervals(frame, keys, "survey", bands,
#'   interval_units = "minutes", site_id = "site")
#' # Observer and interval number, now in frame order
#' unmarked::obsCovs(prepared$data)
#' # Which input row filled each observation row
#' prepared$alignment[, c(".event", ".interval", ".input_row", "survey")]
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
