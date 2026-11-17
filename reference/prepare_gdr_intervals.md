# Add Interval-Level Covariates to a Distance-Removal Frame

Add covariates that can change from one removal interval to the next,
such as who was listening in each interval, to a copy of an
`unmarkedFrameGDR`. Rows are matched by survey identifier and interval
number, so the input can be in any order, and the alignment is recorded.

## Usage

``` r
prepare_gdr_intervals(
  object,
  event_keys,
  event_id,
  interval_data,
  interval_units,
  interval = "interval",
  site_id = NULL,
  effort = NULL,
  expected_event_keys = NULL
)

# S3 method for class 'msr_gdr_intervals'
print(x, ...)
```

## Arguments

- object:

  An `unmarkedFrameGDR` with its original counts and its site, yearly
  and observation covariates.

- event_keys:

  A data frame with one row per survey event, in the frame's order (all
  primary periods of the first frame row, then the second, and so on).
  Use the original survey identifiers. The names `.input_row`,
  `.input_row_name` and `.interval` are reserved.

- event_id:

  Character vector naming the column(s) of `event_keys` that identify
  each survey.

- interval_data:

  A data frame with the `event_id` columns, the interval number column
  and the new covariates. It needs exactly one row for every survey and
  removal interval. Rows can be in any order, and `alignment` records
  how they were matched.

- interval_units:

  Character string naming the time unit of the removal intervals, such
  as `"minutes"`. No time conversion is performed.

- interval:

  Character string, default `"interval"`. The name of the interval
  number column in `interval_data`. Its values must be whole numbers
  from 1 to the number of removal intervals per primary period. The
  column is also added to the observation covariates, so formulas can
  use it.

- site_id:

  Optional column(s) of `event_keys` that identify the physical site,
  passed to
  [`validate_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md).
  Frame rows need not be physical sites.

- effort:

  Optional numeric vector of survey effort, one value per event, passed
  to
  [`validate_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md).
  Positive values mark completed surveys, zero marks unsurveyed events
  and `NA` marks unknown effort.

- expected_event_keys:

  Optional independent event list in frame order, passed to
  [`validate_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md).

- x:

  A prepared `msr_gdr_intervals` object to print.

- ...:

  Not used.

## Value

A list of class `msr_gdr_intervals`:

- `data`: a copy of `object` with the new observation covariates.

- `survey_audit`: the result of
  [`validate_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md)
  for the frame.

- `alignment`: one row per observation row in frame order, with the
  event (`.event`), frame row (`.site_row`), primary period
  (`.primary`), interval (`.interval`), the matching row of
  `interval_data` (`.input_row`, `.input_row_name`) and the `event_keys`
  columns.

- `protocol`: the survey protocol from the audit (units, distance
  breaks, interval lengths), plus the frame dimensions, the
  `interval_column`, the `added_columns`, any added names that also
  occur among the site or yearly covariates
  (`shared_site_or_annual_names`), whether the input order changed
  (`reordered`) and the `unmarked` version.

- `fit_eligibility_established`: always `FALSE`.

Printing reports the dimensions and whether the input order changed, and
returns the object invisibly. If the survey audit finds an error, or the
interval rows cannot be aligned, the function stops with a condition of
class `msr_gdr_interval_error` carrying `survey_audit` and the offending
`input_rows`. The input frame is never changed.

## How rows are aligned

`unmarked` orders observation rows by frame row, primary period and then
removal interval, and its own setter for observation covariates does not
check that order. This function checks the keys first. It rejects
unknown, duplicate, missing or extra survey and interval combinations,
and any attempt to overwrite an existing observation column. The event
identifier columns are used only for matching. Every other column of
`interval_data`, including the interval number, is added.

The counts, existing covariates, distance bands, interval lengths,
units, column classes, factor levels and contrasts are unchanged.
Missing values are kept, and `unmarked` applies its own missing-data
rules when fitting. Nothing is aggregated or imputed. Added names that
also occur among the site or yearly covariates are recorded for review.
Use distinct names for fields with different meanings, because
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
will not estimate a scale for a site covariate whose name is shadowed
this way.

## Model the first interval

Use formulas such as `~ I(interval == 1)` or `~ factor(interval)` when
the observation process justifies them. With equal unit intervals,
hazards pi in the first interval and b thereafter reproduce the
first-time law of a simple encounter Markov process regardless of
persistence after an encounter. First times alone do not estimate that
persistence or establish a finite corrected density. Unequal durations
keep their per-unit probability interpretation. Keep linked complete
histories when assessing re-detections. This function only arranges
data. It adds no model terms, parameters or uncertainty methods.

For covariates that are constant within a survey,
[`expand_gdr_intervals()`](https://wpeterman.github.io/multiScaleR/reference/expand_gdr_intervals.md)
builds `interval_data` for you.

## See also

[`expand_gdr_intervals()`](https://wpeterman.github.io/multiScaleR/reference/expand_gdr_intervals.md),
[`validate_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md),
[`prepare_gdr_histories()`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_histories.md)

## Examples

``` r
frame <- unmarked::unmarkedFrameGDR(
  yDistance = matrix(c(1, 0, 0, 1), 2, byrow = TRUE),
  yRemoval = matrix(c(1, 0, 0, 0, 1, 0), 2, byrow = TRUE),
  siteCovs = data.frame(site = factor(c("A", "B"))),
  dist.breaks = c(0, 25, 75), unitsIn = "m",
  period.lengths = c(1, 2, 1))
keys <- data.frame(site = c("A", "B"), survey = c("s1", "s2"))
# One row per survey and interval, deliberately out of frame order
bands <- data.frame(survey = rep(c("s2", "s1"), each = 3),
  interval = rep(c(3, 1, 2), 2),
  observer = factor(c("Ann", "Ann", "Bob", "Bob", "Ann", "Ann")))
prepared <- prepare_gdr_intervals(frame, keys, "survey", bands,
  interval_units = "minutes", site_id = "site")
# Observer and interval number, now in frame order
unmarked::obsCovs(prepared$data)
#>   interval observer
#> 1        1      Ann
#> 2        2      Ann
#> 3        3      Bob
#> 4        1      Ann
#> 5        2      Bob
#> 6        3      Ann
# Which input row filled each observation row
prepared$alignment[, c(".event", ".interval", ".input_row", "survey")]
#>   .event .interval .input_row survey
#> 1      1         1          5     s1
#> 2      1         2          6     s1
#> 3      1         3          4     s1
#> 4      2         1          2     s2
#> 5      2         2          3     s2
#> 6      2         3          1     s2
```
