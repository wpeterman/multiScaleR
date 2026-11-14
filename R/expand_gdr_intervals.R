#' Expand Site and Survey Covariates into GDR Intervals
#'
#' Copy chosen site covariates, or covariates recorded once per visit, into
#' the observation table of an `unmarkedFrameGDR`, which has one row per
#' survey and removal interval. A removal formula reads its covariates from
#' that table, so a covariate constant within a survey, such as the observer,
#' must be repeated there before you can use it in the removal formula. The
#' function also adds the interval number and returns a copy of the frame
#' with a record of how every row was aligned.
#'
#' @param object An `unmarkedFrameGDR` with its original counts and its site,
#'   yearly and observation covariates.
#' @param event_keys Data frame with one row per survey event, in the frame's
#'   order (all primary periods of the first frame row, then the second, and
#'   so on). Include completed surveys with zero counts. Use the original
#'   survey identifiers.
#' @param event_id Character vector naming the column(s) of `event_keys` that
#'   identify a survey.
#' @param interval_units Character string naming the time unit of the removal
#'   intervals, such as `"minutes"`. Interval lengths are kept, not converted.
#' @param site_covariates Character vector naming the site covariates
#'   (`siteCovs`) to repeat. The default is none. Name an element to give the
#'   copy a new name, for example
#'   `c(rem_forest = "forest", rem_observer = "Observer")`. Unnamed elements
#'   keep the source name. Nothing is inferred from model formulas, so list
#'   every covariate you need.
#' @param yearly_covariates Character vector naming covariates in the frame's
#'   yearly table (`yearlySiteCovs`, one row per primary period) to repeat,
#'   with the same naming convention. The default is none. Despite the name,
#'   these rows are primary periods (visits), not necessarily years.
#' @param interval Character string, default `"interval"`. The name of the
#'   new column holding the interval number (1, 2, ...). It must differ from
#'   the event keys and from the other new columns.
#' @param site_id Optional column(s) of `event_keys` that identify the
#'   physical site, passed to [prepare_gdr_intervals()] and
#'   [validate_unmarked_data()].
#' @param effort Optional numeric vector of survey effort, one value per
#'   event, passed to [validate_unmarked_data()]. Positive values mark
#'   completed surveys, zero marks unsurveyed events and `NA` marks unknown
#'   effort.
#' @param expected_event_keys Optional independent event list in frame
#'   order, passed to [validate_unmarked_data()].
#' @return A list of class `msr_gdr_intervals`, as from
#'   [prepare_gdr_intervals()]: `data` is the copied frame with the new
#'   observation columns, `survey_audit` is the survey audit, `alignment`
#'   records which input row fills each observation row, and `protocol`
#'   records units, interval lengths and the added columns. In addition,
#'   `protocol$expanded_covariates` is a data frame with one row per copied
#'   covariate: `source_table` (`"site"` or `"yearly"`), the `source` column
#'   and its `target` name in the observation table. With no covariates
#'   selected, only the interval column is added. `fit_eligibility_established`
#'   is always `FALSE`.
#' @details
#' # Why repeat covariates
#'
#' When the removal formula contains a random effect, `unmarked` evaluates
#' the random-effect part on the observation table, and that step also needs
#' the formula's fixed covariates there. Repeating only an observer label can
#' therefore leave a
#' formula such as `~ rem_forest + I(interval == 1) + (1 | rem_observer)`
#' without `rem_forest`. Select every survey-constant covariate the
#' observation formula needs.
#'
#' # How values are repeated
#'
#' A site covariate is repeated for every primary period and interval at that
#' site. A yearly covariate is repeated only for the intervals of its own
#' primary period. Column classes, matrix columns, factor levels, contrasts
#' and missing values are kept. Nothing is standardized, imputed or dropped,
#' and the original covariate tables, counts, distance bands and interval
#' lengths are unchanged. Existing observation columns cannot be overwritten.
#' For covariates that change from one interval to the next, use
#' [prepare_gdr_intervals()] directly.
#'
#' # Name copies by role
#'
#' Give a copy a new name when the same field plays different roles, for
#' example `forest` in the abundance formula and `rem_forest` in the removal
#' formula. A copy holds the values present when it was made. If you later
#' change the site covariate, for example while estimating its scale of
#' effect, the copy does not change with it. [multiScale_optim()] also stops
#' if a covariate whose scale it estimates shares its name with a column in
#' the observation or yearly tables, so new names avoid that conflict.
#'
#' Aligning the data does not show that a model fitted to it converges or
#' that its parameters are identifiable. See the "Add covariates for the time
#' intervals" section of
#' `vignette("v12-distance-removal", package = "multiScaleR")` for a worked
#' example.
#' @seealso [prepare_gdr_intervals()], [prepare_gdr_histories()],
#'   [validate_unmarked_data()]
#' @examples
#' # Two surveys, two distance bands and three removal intervals
#' frame <- unmarked::unmarkedFrameGDR(matrix(0, 2, 2), matrix(0, 2, 3),
#'   siteCovs = data.frame(forest = c(.2, .8), Observer = factor(c("A", "B"))),
#'   dist.breaks = c(0, 25, 75), unitsIn = "m", period.lengths = c(1, 2, 1))
#' keys <- data.frame(survey = c("s1", "s2"))
#' # Copy forest and observer into the interval table under new names
#' prepared <- expand_gdr_intervals(frame, keys, "survey", "minutes",
#'   site_covariates = c(rem_forest = "forest", rem_observer = "Observer"))
#' # One row per survey and interval, ready for a removal formula such as
#' # ~ rem_forest + I(interval == 1) + (1 | rem_observer)
#' unmarked::obsCovs(prepared$data)
#' prepared$protocol$expanded_covariates
#' @export
expand_gdr_intervals <- function(object, event_keys, event_id, interval_units,
    site_covariates = character(), yearly_covariates = character(),
    interval = "interval", site_id = NULL, effort = NULL,
    expected_event_keys = NULL) {
  fail <- function(message) stop(message, call. = FALSE)
  if (!methods::is(object, "unmarkedFrameGDR"))
    fail("object must be an unmarkedFrameGDR with its original survey layout.")
  validate_character_scalar(interval, "interval")
  if (!is.data.frame(event_keys) || !is.character(event_id) || !length(event_id) ||
      anyNA(event_id) || anyDuplicated(event_id) || any(!event_id %in% names(event_keys)))
    fail("Supply event_keys and explicit event_id columns in native event order.")
  primary <- object@numPrimary
  events <- as.double(nrow(object@yDistance)) * primary
  bands <- ncol(object@yRemoval) / primary
  total <- events * bands
  if (length(primary) != 1L || !is.finite(primary) || primary < 1 || primary != floor(primary) ||
      !is.finite(bands) || bands < 1 || bands != floor(bands) ||
      !is.finite(total) || total > .Machine$integer.max || nrow(event_keys) != events)
    fail("The frame and keys must define a complete bounded native survey/interval layout.")
  selection <- function(x, table, label) {
    if (!is.character(x) || anyNA(x) || any(!nzchar(x)) || any(!x %in% names(table)))
      fail(paste0(label, " must explicitly name existing covariate columns."))
    target <- names(x)
    if (is.null(target)) target <- unname(x)
    if (anyNA(target)) fail("Covariate aliases must be complete names.")
    target[!nzchar(target)] <- unname(x)[!nzchar(target)]
    data.frame(source_table = rep(label, length(x)), source = unname(x),
      target = target, stringsAsFactors = FALSE)
  }
  site <- unmarked::siteCovs(object)
  yearly <- unmarked::yearlySiteCovs(object)
  map <- rbind(selection(site_covariates, site, "site"),
    selection(yearly_covariates, yearly, "yearly"))
  targets <- c(interval, map$target)
  if (anyDuplicated(targets) || any(targets %in% event_id))
    fail("Expanded aliases, event keys and the interval index must have distinct names.")
  if ((length(site_covariates) && nrow(site) != nrow(object@yDistance)) ||
      (length(yearly_covariates) && nrow(yearly) != events))
    fail("Selected source tables must retain their complete native row layouts.")
  event <- rep(seq_len(events), each = bands)
  interval_data <- event_keys[event, event_id, drop = FALSE]
  rownames(interval_data) <- NULL
  interval_data[[interval]] <- rep(seq_len(bands), events)
  # Data-frame subsetting retains matrix columns and trained factor attributes.
  for (table in c("site", "yearly")) {
    selected <- map[map$source_table == table, , drop = FALSE]
    if (!nrow(selected)) next
    source <- if (table == "site") site else yearly
    rows <- if (table == "site") (event - 1L) %/% primary + 1L else event
    values <- source[rows, selected$source, drop = FALSE]
    names(values) <- selected$target
    rownames(values) <- NULL
    interval_data <- data.frame(interval_data, values, check.names = FALSE)
  }
  out <- prepare_gdr_intervals(object, event_keys, event_id, interval_data,
    interval_units, interval, site_id, effort, expected_event_keys)
  out$protocol$expanded_covariates <- map
  out
}
