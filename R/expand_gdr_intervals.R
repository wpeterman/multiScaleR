#' Expand Site and Survey Covariates into GDR Intervals
#'
#' Repeat explicitly selected covariates into native removal-interval rows and
#' prepare a copy of the GDR frame with its keyed alignment and source aliases.
#' @param object An `unmarkedFrameGDR` with its original paired response,
#'   site, primary-survey and observation covariates.
#' @param event_keys Data frame in native site-major, primary-survey order,
#'   including completed zero surveys. Supply original survey identifiers.
#' @param event_id Character vector of columns identifying a survey in
#'   `event_keys`. Compound keys are supported.
#' @param interval_units Nonempty label for original removal-duration units,
#'   such as `"minutes"`. Durations are retained without conversion.
#' @param site_covariates Character vector naming static site covariates to
#'   repeat. Default: none. Optional names specify new observation-column
#'   aliases, for example `c(rem_forest = "forest", rem_observer = "Observer")`.
#'   Unnamed entries retain their source names. Selection is explicit; no
#'   formula variables or observer roles are inferred.
#' @param yearly_covariates Character vector naming native primary-survey
#'   covariates to repeat, with the same optional alias convention. Default:
#'   none. These rows refer to native primary periods, not necessarily years.
#' @param interval One string naming the retained numeric removal-interval
#'   index. Default: `"interval"`. It must not collide with an event key or
#'   another requested observation column.
#' @param site_id Optional physical-site identifier columns in `event_keys`,
#'   passed to [prepare_gdr_intervals()].
#' @param effort Optional event-order effort vector, passed to
#'   [prepare_gdr_intervals()]. Surveyed zeros and unsurveyed missing responses
#'   keep their original meanings.
#' @param expected_event_keys Optional independently supplied event mapping
#'   in native order, passed to [prepare_gdr_intervals()].
#' @return An `msr_gdr_intervals` list as returned by [prepare_gdr_intervals()],
#'   with the copied frame in `data`, survey audit, complete alignment and
#'   protocol. `protocol$expanded_covariates` is a data frame with `source_table`
#'   (`"site"` or `"yearly"`), original `source` name and observation-column
#'   `target` alias. Its rows identify every repeated field. An empty selection
#'   adds only the interval index. `fit_eligibility_established` stays `FALSE`.
#' @details
#' The native removal random-term model frame can evaluate fixed covariates
#'   from the observation table as well as its grouping factor. Repeating only
#'   an observer label can therefore leave a formula such as
#'   `~ rem_forest + I(interval == 1) + (1 | rem_observer)` incomplete. Select
#'   all survey-constant fields needed by that observation formula.
#'
#' Static site fields repeat across every primary survey and interval for the
#'   site. Primary-survey fields repeat only within their own survey. Original
#'   column classes, matrix columns, factor levels, contrasts and missing values
#'   are preserved. No model is fitted, no predictor is standardized and no
#'   missingness decision is made here. Original covariate tables, responses,
#'   geometry and durations remain unchanged.
#'
#' Use aliases for different component roles. For example, `forest` in abundance
#'   and `rem_forest` in removal are distinct inputs. An alias copies the values
#'   present at preparation; it does not keep later changes synchronized. The
#'   existing guard against optimizing a static site name shadowed by an
#'   observation or primary-survey column remains in force. Existing observation
#'   columns cannot be overwritten. Interval-varying measurements can instead
#'   be supplied directly through [prepare_gdr_intervals()].
#'
#' Data alignment does not qualify native TMB controls, likelihood-only fitting,
#'   latent integration, start agreement, score, covariance or identification.
#'   Priors, parameter counts, diagnostics and uncertainty methods are unchanged.
#' @seealso [prepare_gdr_intervals()], [prepare_gdr_histories()],
#'   [validate_unmarked_data()]
#' @examples
#' frame <- unmarked::unmarkedFrameGDR(matrix(0, 2, 2), matrix(0, 2, 3),
#'   siteCovs = data.frame(forest = c(.2, .8), Observer = factor(c("A", "B"))),
#'   dist.breaks = c(0, 25, 75), unitsIn = "m", period.lengths = c(1, 2, 1))
#' keys <- data.frame(survey = c("s1", "s2"))
#' prepared <- expand_gdr_intervals(frame, keys, "survey", "minutes",
#'   site_covariates = c(rem_forest = "forest", rem_observer = "Observer"))
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
