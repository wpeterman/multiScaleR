# Convert Linked Bird Records into Paired Distance-Removal Counts

Turn a matrix of bird records, one row per individual at a survey, into
the distance and removal count tables of a distance-removal model. The
function finds each bird's first detection across all distance bands,
including bands beyond the survey radius, before it applies the radius.
It keeps every source record, reports what happened to each bird, and
never deletes or repairs a record.

## Usage

``` r
prepare_gdr_histories(
  histories,
  history_keys,
  event_keys,
  event_id,
  individual_id,
  cell_map,
  dist_breaks,
  radius,
  interval_lengths,
  interval_units,
  units_in,
  surveyed,
  expected_event_keys = NULL,
  native_frame = NULL
)

# S3 method for class 'msr_gdr_histories'
print(x, ...)
```

## Arguments

- histories:

  Numeric or logical matrix with one row per individual bird at a survey
  and one named column per source cell (one distance band in one time
  interval). Enter 1 where the bird was recorded and 0 elsewhere. Any
  other value, including `NA`, is kept as it is and reported in
  `findings`; it is never treated as 0.

- history_keys:

  Data frame with one row per row of `histories`, in the same order,
  holding the original survey and individual identifiers and any other
  notes you want to keep.

- event_keys:

  Data frame with one row per survey event, in the order the count
  tables should have. Include completed surveys that recorded no birds.
  Event keys must be unique.

- event_id:

  Character vector naming the column(s) that identify a survey in both
  `event_keys` and `history_keys`.

- individual_id:

  Character vector naming the column(s) of `history_keys` that identify
  a bird within a survey. Each survey and individual combination must be
  unique.

- cell_map:

  Data frame with columns `column` (a column name of `histories`),
  `interval` and `distance_bin` (positive whole numbers). It must map
  every column of `histories` to exactly one combination of interval and
  distance band, and every combination must appear once. Columns of
  `histories` can be in any order.

- dist_breaks:

  Increasing distance boundaries of the source bands, starting at zero,
  in `units_in`. Include bands beyond the survey radius when the records
  have them. Only the last boundary may be `Inf`.

- radius:

  Survey radius, in `units_in`. It must equal one of `dist_breaks`,
  because a band cannot be split without finer distances.

- interval_lengths:

  Positive whole numbers giving the duration of each time interval, in
  `interval_units`. Intervals may differ in length.

- interval_units:

  Character string naming the time unit, such as `"seconds"`. It is
  recorded, not converted.

- units_in:

  Distance units of `dist_breaks` and `radius`, `"m"` or `"km"`.

- surveyed:

  Logical vector with one value per event, in event order, with no `NA`:
  `TRUE` for a completed survey and `FALSE` for an event that was not
  surveyed. A recorded bird at an unsurveyed event is an error.

- expected_event_keys:

  Optional data frame, built independently, with the `event_id` columns
  in the intended event order. The function stops if `event_keys`
  differs from it.

- native_frame:

  Optional single-visit `unmarkedFrameGDR` built earlier from these
  records. It is audited with
  [`validate_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md),
  and the function stops unless its two count tables, its distance
  boundaries up to the radius and its interval lengths match the
  prepared results exactly. Distances in meters and kilometers are
  compared in meters, without changing either input.

- x:

  A prepared `msr_gdr_histories` object.

- ...:

  Not used.

## Value

A list of class `msr_gdr_histories`:

- `histories`, `history_keys`, `event_keys` and `cell_map`: the inputs,
  unchanged.

- `first`: one row per bird, with its row in `histories`
  (`.history_row`), its event (`.event`), the first interval with a
  recorded mark (`candidate_interval`), the resolved first detection
  (`first_interval`, `first_distance_bin`; `NA` when unresolved),
  whether the bird is counted (`primary_eligible`), whether every cell
  holds a valid 0 or 1 with at most one band per interval
  (`history_complete`), whether both hold (`complete_history_eligible`)
  and its `status` (see "Read the status of each bird").

- `events`: one row per event, with `surveyed`, `response_complete`
  (whether its counts are observed), `known_first_records` (birds
  counted) and `unresolved_first_histories`.

- `findings`: one row per problem, with the history row, the cell
  (`.cell`, `.column`, `interval`, `distance_bin`), a `code` and the
  original `value`. Findings about a whole history have missing cell
  fields.

- `first_cells`: counts of first detections by event (rows) and cell
  inside the radius (columns), ordered by distance band, then interval.

- `native_cells`: the rows of `cell_map` for those columns, in the same
  order.

- `yDistance` and `yRemoval`: the two count tables, events by distance
  band and events by interval, ready for
  [`unmarked::unmarkedFrameGDR()`](https://ecoverseR.github.io/unmarked/reference/unmarkedFrameGDR.html).

- `protocol`: the source and in-radius distance breaks, the radius,
  units, interval lengths and the column order of `first_cells`.

- `surveyed`, `survey_audit` (the audit of `native_frame`, or `NULL`),
  `fingerprint` (a SHA-256 of the inputs) and `summary` (counts of
  histories, events, counted birds, complete histories, invalid cells,
  unresolved events and birds first detected outside the radius).

- `fit_eligibility_established` and
  `density_identification_established`: always `FALSE`. Preparing counts
  does not show that a model fitted to them is valid.

Printing reports the summary counts and returns the object invisibly.
Errors in the inputs stop with a condition of class
`msr_gdr_history_error` whose `history_rows` element lists the rows
involved, when they are known.

## Find each bird's first detection

The first detection is found across all source bands before the radius
is applied. A bird first recorded beyond the radius stays out of the
counts, even if a later mark places it near the point. If a bird has
marks in more than one distance band in the same interval, its distance
in that interval is not unique. Such marks are kept for review and never
collapsed into one band. They may call for a model of movement or of the
recording process.

A bird is counted when its first detection is known and inside the
radius, even if a later mark is malformed. That bird is then left out of
the complete-history subset. When a malformed or multi-band mark falls
at or before a bird's first detection, the first detection is
unresolved, as it is for a history with a malformed mark and no
detection. The whole event's counts then become missing, so that an
incomplete count cannot pass as an observed zero. The known partial
count stays in `events$known_first_records`, and the source rows stay in
`histories`.

## Read the status of each bird

|  |  |
|----|----|
| `status` | Meaning |
| `retained_first` | Counted. The history is complete |
| `retained_first_incomplete_history` | Counted, but a later mark is malformed or spans two bands |
| `first_outside_radius` | Not counted. The bird was first detected beyond the radius |
| `unresolved_first` | The first detection cannot be determined. The event's counts are missing |
| `no_record` | No detections were recorded, and every cell is a valid 0 |
| `unsurveyed_history` | The history belongs to an event that was not surveyed |

`findings` uses the codes `invalid_mark` (a cell that is not 0 or 1,
with its `value`), `multiple_distance_bins` (marks in more than one band
in an interval), and `unresolved_first`, `first_outside_radius` and
`unsurveyed_history` as in the table.

## Read the count tables

Summing the columns of `first_cells` over intervals gives `yDistance`,
and summing over distance bands gives `yRemoval`, so the two tables
always have equal row totals. A completed survey with no counted birds
has a row of zeros. An unsurveyed event, or an event with an unresolved
bird, has a row of `NA`.

## Keep the other responses separate

The paired first detections are one response. Complete linked histories
of the counted birds, re-detections, and marks beyond the radius are
others. Their likelihoods have different constants and conditioning, so
they cannot be compared in one AIC table. The time intervals within one
count are also not independent visits. Keep the source records for any
separate analysis. Better records of the detected birds do not by
themselves identify the birds that were never detected.

See the "Prepare your own records" section of
[`vignette("v12-distance-removal", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v12-distance-removal.md)
for a worked example.

## See also

[`validate_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md),
[`prepare_gdr_intervals()`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_intervals.md),
[`expand_gdr_intervals()`](https://wpeterman.github.io/multiScaleR/reference/expand_gdr_intervals.md)

## Examples

``` r
# Source cells: 2 time intervals by 3 distance bands (0-25 m, 25-75 m and
# beyond the 75 m radius)
map <- expand.grid(interval = 1:2, distance_bin = 1:3)
map$column <- paste0("cell", 1:6)
h <- matrix(0, 2, 6, dimnames = list(NULL, map$column))
h[1, c(5, 2)] <- 1 # Beyond 75 m first, then near: not counted
h[2, 1] <- 1       # Interval 1, 0-25 m: counted
h[2, 2] <- 2       # Malformed later mark: the valid first still counts
prepared <- prepare_gdr_histories(h,
  history_keys = data.frame(survey = c("a", "a"), bird = c("one", "two")),
  event_keys = data.frame(survey = c("a", "b", "c")), event_id = "survey",
  individual_id = "bird", cell_map = map, dist_breaks = c(0, 25, 75, Inf),
  radius = 75, interval_lengths = c(75, 75), interval_units = "seconds",
  units_in = "m", surveyed = c(TRUE, TRUE, FALSE))
# What happened to each bird
prepared$first[, c("first_interval", "first_distance_bin", "status")]
#>   first_interval first_distance_bin                            status
#> 1              1                  3              first_outside_radius
#> 2              1                  1 retained_first_incomplete_history
# Survey b is an observed zero. Survey c was not surveyed, so it is NA
prepared$yDistance
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    0
#> [3,]   NA   NA
prepared$yRemoval
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    0
#> [3,]   NA   NA
# Both tables are totals of the same first-detection cells
prepared$first_cells
#>      cell1 cell2 cell3 cell4
#> [1,]     1     0     0     0
#> [2,]     0     0     0     0
#> [3,]    NA    NA    NA    NA
```
