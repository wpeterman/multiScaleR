# Estimate Memory Requirements for Multiscale Optimization

Estimates how much memory a
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
object (and optionally a fitted starting model) needs during
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
and recommends a conservative maximum number of parallel workers
(`n_cores`). Use it before running a large or slow model in parallel,
because each worker holds its own copy of the inputs.

## Usage

``` r
estimate_multiscale_ram(
  kernel_inputs,
  fitted_mod = NULL,
  join_by = NULL,
  refit_fn = NULL,
  n_cores = NULL,
  PSOCK = (.Platform$OS.type != "unix"),
  safety_factor = 1.5,
  ram_fraction = 0.75,
  leave_free = 2L
)

# S3 method for class 'multiScaleR_ram'
print(x, ...)
```

## Arguments

- kernel_inputs:

  A `multiScaleR_data` object created by
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  or
  [`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md).

- fitted_mod:

  Optional fitted starting model, as you would pass to
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
  When supplied, the estimate also includes the optimization context
  sent to each worker. Without it, the estimate can understate the real
  payload. Default: `NULL`.

- join_by:

  Optional data frame passed through to
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  for `unmarked` workflows. Default: `NULL`.

- refit_fn:

  Optional custom refit function passed through to
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
  Default: `NULL`.

- n_cores:

  Optional positive integer. Worker count at which to report the peak
  memory estimate. When omitted, the recommended maximum is used.

- PSOCK:

  Logical. Budget workers as PSOCK (socket) workers. Default: `TRUE` on
  Windows and `FALSE` on Unix-like systems. The estimate is the same
  either way, because forked workers are also budgeted as though each
  holds its own copy.

- safety_factor:

  Numeric scalar, at least 1. Multiplier on each worker's payload to
  allow for temporary allocations during optimization. Larger values
  give fewer recommended workers. Default: `1.5`.

- ram_fraction:

  Numeric scalar between 0 and 1. Share of total system RAM treated as
  usable. The rest is left for the operating system and other R objects.
  Default: `0.75`.

- leave_free:

  Non-negative integer. Number of physical CPU cores to leave unused
  when computing the CPU limit. Default: `2`.

- x:

  A `multiScaleR_ram` object returned by `estimate_multiscale_ram()`.

- ...:

  Not used.

## Value

`estimate_multiscale_ram()` returns a list of class `multiScaleR_ram`
with these elements:

- `requested_n_cores`: the worker count used for the peak estimate
  (`n_cores`, or the recommendation when `n_cores` is omitted).

- `recommended_n_cores`: the smaller of the CPU and RAM limits, or the
  CPU limit when RAM could not be detected. Never less than 1.

- `max_cores_by_cpu`: physical cores minus `leave_free`, at least 1.

- `max_cores_by_ram`: the number of workers that fit in the usable RAM,
  at least 1, or `NA` when total RAM could not be detected.

- `backend`: `"PSOCK"` or `"fork"`, the assumed parallel backend.

- `system`: list with the detected `physical_cores` and `logical_cores`,
  `leave_free`, `total_ram_bytes`, `usable_ram_fraction`,
  `usable_ram_bytes`, `safety_factor`, and readable versions of the RAM
  values (`total_ram_pretty`, `usable_ram_pretty`).

- `point_summary`: one-row data frame with the number of points, the
  mean, median, and maximum number of extracted cells per point (`NA` in
  lean mode), the number of source layers, whether cell-level data are
  stored (`cell_data_stored`), and the number of distance bins
  (`nbins`).

- `component_bytes`: data frame with columns `component`, `bytes`, and
  `pretty` (a readable size) for each stored object
  (`kernel_inputs_total`, `kernel_dat`, `d_list`, `raw_cov`, `binned`,
  `join_by`, `fitted_mod`, `opt_context`), the two bundles described in
  Details (`worker_bundle` and `main_bundle`), and
  `peak_parallel_estimate`.

- `notes`: character vector of interpretation notes and detection
  warnings.

[`print()`](https://rdrr.io/r/base/print.html) displays the main limits,
the detected hardware, the peak estimate, and the notes, and returns `x`
invisibly.

## How the estimate is built

Sizes come from
[`utils::object.size()`](https://rdrr.io/r/utils/object.size.html), so
they describe stored R objects, not the transient memory used by
`exactextractr`, model refits, or linear algebra libraries. The *main
bundle* is what the main R session holds (the kernel inputs, `join_by`,
the optimization context, and the model). The *worker bundle* is what
each parallel worker receives (the per-cell distances and values,
`join_by`, the optimization context, and the model). The limits and the
peak estimate follow from these sizes.

- Usable RAM = total RAM x `ram_fraction`.

- RAM limit = (usable RAM - main bundle) / (worker bundle x
  `safety_factor`), rounded down.

- Peak estimate = main bundle + `n_cores` x worker bundle x
  `safety_factor`.

Even where forking lets workers share memory, each worker is budgeted as
if it held its own copy. Hardware measurements that succeed are kept for
the rest of the R session. A failed RAM lookup is retried on the next
call.

## Read the report

Compare the worker count you plan to use with `max_cores_by_ram`, and
read the peak estimate at that count. These are conservative estimates,
not measured memory use. The recommendation is never below one worker,
but that does not guarantee a serial run fits in memory. A run at or
below the limit can still run out of memory during a model refit. If RAM
detection or the estimate fails, no limit can be enforced.

## Reduce the memory needed

When every covariate is a kernel mean, rebuild the original
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
call with `bin = TRUE, store_cell_data = FALSE`, keeping its other
settings, and check the estimate again. This lean mode keeps only the
distance-binned summaries that optimization uses. It is not available
for landscape or surface metrics, which need the cell-level data.
Otherwise, use fewer workers or `n_cores = 1`. Do not delete components
of the object by hand, and do not shrink a scientifically chosen `max_D`
only to pass the check.

## How multiScale_optim() uses the estimate

With `n_cores > 1`,
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
calls this function with its default settings before starting workers.
It stops with an error when `n_cores` exceeds `max_cores_by_ram`. It
does not lower `n_cores` for you. When RAM cannot be detected or the
estimate fails, the check does not block the run.

## See also

[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md).
The section "Estimate the scales of effect" in
[`vignette("v02-user-guide", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v02-user-guide.md)
shows the check on a worked example.

## Examples

``` r
data("pts")
data("count_data")
hab <- terra::rast(system.file("extdata", "hab.tif",
                               package = "multiScaleR"))
kernel_inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
                             kernel = "gaussian", verbose = FALSE)
mod <- glm(y ~ hab, family = poisson,
           data = data.frame(count_data["y"], kernel_inputs$kernel_dat))

# Recommended worker limits and the peak memory estimate for two workers
ram <- estimate_multiscale_ram(kernel_inputs, fitted_mod = mod, n_cores = 2)
ram
#> 
#> Estimated parallel RAM budget for `multiScaleR`
#> Backend assumption: fork
#> Recommended max workers: 2
#> CPU limit after reserve: 2
#> RAM limit: 2975
#> Physical cores detected: 4
#> Logical cores detected: 4
#> Total RAM detected: 15.6 Gb
#> Usable RAM budget: 11.7 Gb
#> Peak estimate at 2 worker(s): 11.1 Mb
#> 
#> Notes:
#> - Parallel RAM budgeting remains conservative even when forked workers may share memory.

# Estimated size of each component, largest first
ram$component_bytes[order(-ram$component_bytes$bytes), ]
#>                 component    bytes   pretty
#> 11 peak_parallel_estimate 11591216  11.1 Mb
#> 10            main_bundle  3143816     3 Mb
#> 9           worker_bundle  2815800   2.7 Mb
#> 1     kernel_inputs_total  2547720   2.4 Mb
#> 3                  d_list  1243160   1.2 Mb
#> 4                 raw_cov   976288 953.4 Kb
#> 8             opt_context   461744 450.9 Kb
#> 5                  binned   316504 309.1 Kb
#> 7              fitted_mod   133832 130.7 Kb
#> 2              kernel_dat     6120     6 Kb
#> 6                 join_by        0  0 bytes
```
