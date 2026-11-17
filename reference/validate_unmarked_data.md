# Audit Survey Identities and Paired Distance-Removal Counts

Check an `unmarkedFrameGDR` against a table of original survey
identifiers before you fit a model. The audit looks for split or
duplicated surveys, covariates out of line with the surveys, distance
and removal counts that disagree, and counts that conflict with survey
effort. It only reports. It does not combine surveys, reorder rows,
remove records or repair counts.

## Usage

``` r
validate_unmarked_data(
  object,
  event_keys,
  event_id,
  site_id = NULL,
  interval_units,
  effort = NULL,
  expected_event_keys = NULL
)

# S3 method for class 'msr_unmarked_audit'
print(x, ...)
```

## Arguments

- object:

  An `unmarkedFrameGDR` holding paired distance and removal counts for
  point surveys. Other frame classes are not supported.

- event_keys:

  A data frame with one row per survey event, in the frame's order. An
  event is one primary period at one frame row, so with one visit per
  point there is one row per frame row. With several primary periods,
  list all periods of the first frame row, then all periods of the
  second, and so on. Use the original survey identifiers, such as a
  point name and date, rather than a row number.

- event_id:

  Character vector naming the columns of `event_keys` that together
  identify one survey event. Their combinations should be unique;
  repeated combinations are reported.

- site_id:

  Optional character vector naming the columns of `event_keys` that
  identify the physical site. A site may appear in several events
  (repeat visits). The default `NULL` skips the site checks.

- interval_units:

  Character string naming the time unit of `object@period.lengths`, for
  example `"minutes"`. Removal probabilities are per unit of this time.
  The label is recorded, not converted.

- effort:

  Optional numeric vector with one value per event, in event order, such
  as minutes surveyed. A positive value marks a completed survey, zero
  marks an event that was not surveyed (its counts must be `NA`), and
  `NA` marks unknown effort. With the default `NULL`, the audit adds a
  warning, because counts alone cannot separate a completed survey with
  no detections from one that never happened.

- expected_event_keys:

  Optional data frame, built independently of `event_keys` (for example
  from the field schedule), with the `event_id` columns in the order the
  frame should have. Rows whose keys differ are reported. A key table
  cannot show by itself that it is in the right order.

- x:

  An object returned by `validate_unmarked_data()`.

- ...:

  Further arguments passed to the data-frame print method.

## Value

A list of class `msr_unmarked_audit`:

- `valid`: `TRUE` when no errors were found. Warnings can remain, so
  read them too.

- `findings`: a data frame with one row per finding. Columns are
  `severity` (`"error"` or `"warning"`), `code` (see "Finding codes"),
  `message`, the event's position in frame order (`.event`), its frame
  row (`.site_row`) and primary period (`.primary`), followed by the
  `event_keys` columns for that event. Findings about the whole frame
  have missing positions.

- `events`: `event_keys` as supplied.

- `summary`: the number of frame rows (`sites`), `primary_periods`,
  `events`, `errors` and `warnings`, and `missing_histories`, the number
  of events with both count tables fully missing.

- `protocol`: `distance_units`, `interval_units`, `distance_breaks` and
  `interval_lengths` from the frame.

- `effort`: the supplied effort vector.

Printing shows the error and warning counts and at most ten findings,
and returns the object invisibly. The full table stays in `$findings`.

## Identify survey events and physical sites

A physical site can be visited more than once. Give each visit its own
event identifier, and keep the site in a separate column named by
`site_id`. A row number from a spreadsheet cannot show that a detection
record and a zero-detection record belong to the same survey, so resolve
survey identities when you build the frame from the original records,
for example with
[`prepare_gdr_histories()`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_histories.md).

The distance and removal tables count the same birds, so their totals
must match for each event. Matching totals are necessary but cannot
reveal a survey split across two frame rows. Supplying the original
event keys lets the audit report such a split as `duplicate_event`. Key
names that look like row indices (such as `index`, `row_id` or `X`) give
a warning. When an `event_id` or `site_id` column also appears in the
site, yearly or observation covariates, its values are compared with
`event_keys` row by row. Without `expected_event_keys` or such shared
columns, unique keys alone cannot prove that the rows are aligned.

## Record survey effort

A survey that took place and detected nothing is an observed zero. A
survey that never took place should have fully missing counts in both
tables. Supply `effort` so the audit can tell these apart. Counts that
are missing in only one table, or in only some columns, are reported
because their totals cannot be compared. Character and factor values,
including the literal string `"None"`, are kept as values and not
treated as missing.

## Finding codes

|  |  |  |
|----|----|----|
| Code | Severity | Meaning |
| `frame_layout` | error | The count tables do not have matching rows and complete primary-period blocks. Other checks are skipped |
| `survey_type` | error | The frame is not a point survey |
| `event_layout` | error | `event_keys` does not have one row per event. Other checks are skipped |
| `missing_event_key` | error | An event identifier is missing or blank |
| `duplicate_event` | error | One event key appears in more than one event, as when a survey is split |
| `index_like_event_key` | warning | An `event_id` column name looks like a row index |
| `missing_site_key` | error | A site identifier is missing or blank |
| `physical_site_alignment` | error | Primary periods in one frame row have different site keys |
| `event_alignment` | error | The key differs from `expected_event_keys` in that row |
| `duplicate_reference_event` | error | `expected_event_keys` repeats an event |
| `covariate_layout` | error | A covariate table has the wrong number of rows or repeated names |
| `covariate_key_alignment` | error | A key column in a covariate table differs from `event_keys` |
| `response_slot_alignment` | error | The frame's response and removal tables disagree |
| `distance_breaks` | error | Distance breaks do not start at zero, increase, and match the distance columns |
| `distance_units` | error | Distance units are not `"m"` or `"km"` |
| `interval_lengths` | error | Interval lengths are not positive whole numbers, one per removal column |
| `invalid_counts` | error | A count is negative, fractional, infinite or `NaN` |
| `partial_history` | error | An event's counts are only partly missing |
| `unequal_totals` | error | Distance and removal totals differ for the event |
| `effort_unspecified` | warning | No `effort` was supplied |
| `invalid_effort` | error | Effort is negative or infinite |
| `unsurveyed_counts` | error | An event with zero effort has counts, including zeros |
| `surveyed_missing` | error | An event with positive effort has fully missing counts |
| `unknown_effort` | warning | Effort is `NA` for the event |

The audit does not decide how to treat a bird first detected beyond the
survey radius and later inside it, or how to turn raw bird records into
counts.
[`prepare_gdr_histories()`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_histories.md)
handles that step. See the "Check that the records line up with the
surveys" section of
[`vignette("v12-distance-removal", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v12-distance-removal.md)
for the audit in a full workflow.

## See also

[`prepare_gdr_histories()`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_histories.md),
[`prepare_gdr_intervals()`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_intervals.md),
[`expand_gdr_intervals()`](https://wpeterman.github.io/multiScaleR/reference/expand_gdr_intervals.md),
[`sim_dat_gdr()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat_gdr.md)

## Examples

``` r
# Two visits to site A, one frame row each, with matching totals
frame <- unmarked::unmarkedFrameGDR(
  yDistance = matrix(c(1, 0, 0, 0), 2, byrow = TRUE),
  yRemoval = matrix(c(0, 1, 0, 0), 2, byrow = TRUE),
  dist.breaks = c(0, 25, 75), unitsIn = "m",
  period.lengths = c(1, 2))
keys <- data.frame(site = c("A", "A"), visit = c("spring", "summer"))
audit <- validate_unmarked_data(frame, keys, c("site", "visit"),
  site_id = "site", interval_units = "minutes", effort = c(3, 3))
audit$valid      # TRUE: no errors and, with effort supplied, no warnings
#> [1] TRUE
audit
#> GDR survey audit: 2 events; 0 errors; 0 warnings

# Give both rows the same key, as if one survey had been split in two.
# Their totals still match, but the event mapping is now invalid.
split_keys <- keys
split_keys$visit[2] <- "spring"
split_audit <- validate_unmarked_data(frame, split_keys, c("site", "visit"),
  site_id = "site", interval_units = "minutes", effort = c(3, 3))
split_audit$valid
#> [1] FALSE
split_audit$findings[, c("severity", "code", ".event", "visit")]
#>   severity            code .event  visit
#> 1    error duplicate_event      1 spring
#> 2    error duplicate_event      2 spring
```
