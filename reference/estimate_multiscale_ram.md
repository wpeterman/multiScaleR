# Estimate Memory Requirements for Multiscale Optimization

Estimate the memory footprint of a
[`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
object and provide a conservative recommendation for parallel worker
counts. The helper inspects the stored kernel inputs, optionally builds
the same optimization context used by
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
detects the number of physical CPU cores, and estimates a maximum
`n_cores` value while leaving a user-defined number of cores free.

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
```

## Arguments

- kernel_inputs:

  A `"multiScaleR_data"` object created by
  [`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md).

- fitted_mod:

  Optional fitted model object. When supplied, the helper also measures
  the optimization context that would be serialized to parallel workers
  during
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

- join_by:

  Optional data frame passed through to
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  for `unmarked` workflows.

- refit_fn:

  Optional custom refit function passed through to
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

- n_cores:

  Optional positive integer. When supplied, the helper also reports the
  estimated peak memory for that specific worker count. When omitted,
  the reported peak estimate uses the recommended maximum worker count.

- PSOCK:

  Logical. Should parallel runs be budgeted as PSOCK workers. Default:
  `TRUE` on non-Unix platforms and `FALSE` on Unix-like platforms. Even
  on systems that support forking, the recommendation uses a
  conservative per-worker duplication assumption.

- safety_factor:

  Numeric scalar \\\ge 1\\. Multiplier applied to the per-worker payload
  to allow for temporary allocations during optimization. Default:
  `1.5`.

- ram_fraction:

  Numeric scalar in `(0, 1]`. Fraction of detected total system RAM
  considered usable for this run. The remainder is reserved as headroom
  for the OS and other R objects. Default: `0.75`.

- leave_free:

  Integer scalar \\\ge 0\\. Number of physical CPU cores to leave unused
  when computing the recommended maximum worker count. Default: `2`.

## Value

A list of class `"multiScaleR_ram"` with:

- `requested_n_cores`:

  Worker count used for the `peak_parallel_estimate`.

- `recommended_n_cores`:

  Conservative recommended maximum worker count after applying both CPU
  and RAM limits.

- `max_cores_by_cpu`:

  Maximum worker count allowed by detected physical CPU cores after
  leaving `leave_free` cores unused.

- `max_cores_by_ram`:

  Maximum worker count allowed by the RAM budget, or `NA` if total
  system RAM could not be detected.

- `backend`:

  Character string identifying the assumed parallel backend budget:
  `"PSOCK"` or `"fork"`.

- `system`:

  Named list describing detected physical cores, logical cores, total
  RAM, usable RAM, and the reservation settings.

- `point_summary`:

  One-row data frame summarizing the number of points, extracted cells
  per point, source raster layer count, whether cell-level data were
  retained (`cell_data_stored`), and the number of distance bins
  (`nbins`).

- `component_bytes`:

  Data frame reporting estimated byte sizes for the major stored objects
  and worker bundles.

- `notes`:

  Character vector with interpretation notes and any detection warnings.

## Details

The estimate is advisory rather than exact. It is based on
[`utils::object.size()`](https://rdrr.io/r/utils/object.size.html) and
therefore reflects the size of stored R objects, not all transient
allocations made by `exactextractr`, model refits, or BLAS/LAPACK.

The recommendation is intentionally conservative. Even on Unix-like
systems where forking can share memory copy-on-write, worker counts are
budgeted as though each worker may require its own copy of the
optimization payload.

If `fitted_mod` is omitted, the RAM recommendation is based on
`kernel_inputs` alone and may understate the true optimization payload.

[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
stops before starting parallel workers when `n_cores > 1` exceeds the
RAM-based limit under the default settings of this helper. Use
`n_cores = 1` for serial execution, or rebuild
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
inputs with `bin = TRUE` and `store_cell_data = FALSE` when all
covariates are kernel-weighted. Landscape and surface metrics require
cell-level data and cannot use this lean mode. Recheck the estimate
after rebuilding inputs.

Worker recommendations have a minimum of one; this does not guarantee
that serial execution fits in memory. The budget check cannot enforce a
limit when RAM detection or estimation fails, and a run below the
estimated limit can still exhaust memory during a model refit.

## Examples

``` r
library(terra)

pts <- vect(cbind(c(3, 5, 7),
                  c(7, 5, 3)))
r <- rast(matrix(rnorm(100), nrow = 10))
names(r) <- "hab"

kernel_inputs <- kernel_prep(
  pts = pts,
  raster_stack = r,
  max_D = 2,
  kernel = "gaussian",
  verbose = FALSE
)

estimate_multiscale_ram(kernel_inputs)
#> 
#> Estimated parallel RAM budget for `multiScaleR`
#> Backend assumption: fork
#> Recommended max workers: 2
#> CPU limit after reserve: 2
#> RAM limit: 2053848
#> Physical cores detected: 4
#> Logical cores detected: 4
#> Total RAM detected: 15.6 Gb
#> Usable RAM budget: 11.7 Gb
#> Peak estimate at 2 worker(s): 38.5 Kb
#> 
#> Notes:
#> - No fitted model was supplied, so the recommendation is based on `kernel_inputs` alone.
#> - Parallel RAM budgeting remains conservative even when forked workers may share memory.
```
