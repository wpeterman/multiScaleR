# Simulate spatially autocorrelated raster surfaces

Creates a `SpatRaster` stack of four simulated landscape surfaces for
use with
[`sim_dat`](https://wpeterman.github.io/multiScaleR/reference/sim_dat.md)
and
[`sim_dat_unmarked`](https://wpeterman.github.io/multiScaleR/reference/sim_dat_unmarked.md).
The stack contains two binary (0/1) and two continuous (0–1) surfaces
with differing spatial autocorrelation ranges, allowing simulation of
multiscale ecological processes.

## Usage

``` r
sim_rast(
  dim = 100,
  resolution = 10,
  autocorr_range1 = NULL,
  autocorr_range2 = NULL,
  sill = 10,
  plot = FALSE,
  user_seed = NULL,
  ...
)
```

## Arguments

- dim:

  Positive integer. Number of cells on each side of the square raster.
  The output will be a `dim x dim` grid. Default: `100`.

- resolution:

  Positive numeric. Cell size (edge length) in map units. The raster
  extent will be `dim * resolution` in each direction. Default: `10`.

- autocorr_range1:

  Optional positive numeric. Spatial autocorrelation range (in map
  cells) for surfaces 1 and 3 (`bin1`, `cont1`). Controls the decay rate
  of the exponential covariance: larger values produce smoother, more
  broadly correlated patterns. If `NULL` (default), set to 5% of `dim`
  (fine-scale autocorrelation).

- autocorr_range2:

  Optional positive numeric. Autocorrelation range for surfaces 2 and 4
  (`bin2`, `cont2`). If `NULL` (default), set to 25% of `dim`
  (broad-scale autocorrelation).

- sill:

  Positive numeric. Partial sill (variance) of the Gaussian random field
  used for all four surfaces. Default: `10`.

- plot:

  Logical. If `TRUE`, the raster stack is plotted using
  [`terra::plot`](https://rspatial.github.io/terra/reference/plot.html)
  after simulation. Default: `FALSE`.

- user_seed:

  Optional integer seed for reproducibility. Different transformations
  of `user_seed` are applied to each of the four surfaces so they are
  independent but fully reproducible. Default: `NULL`.

- ...:

  Additional arguments. Not currently used.

## Value

A `SpatRaster` with four named layers:

- `bin1`:

  Binary (0/1) surface with fine-scale autocorrelation
  (`autocorr_range1`). Values of 1 where the underlying continuous field
  exceeds its 55th percentile.

- `bin2`:

  Binary (0/1) surface with broad-scale autocorrelation
  (`autocorr_range2`). Values of 1 where the underlying continuous field
  falls below its 40th percentile.

- `cont1`:

  Continuous surface rescaled to \[0, 1\] with fine-scale
  autocorrelation (75% of `autocorr_range1`).

- `cont2`:

  Continuous surface rescaled to \[0, 1\] with broad-scale
  autocorrelation (125% of `autocorr_range2`).

All layers span the same spatial extent:
`[0, dim * resolution] x [0, dim * resolution]`.

## Details

Each surface is generated as a Gaussian random field using a fast
Fourier transform (FFT) circulant embedding approach with an exponential
covariance function. The underlying continuous fields are normalized to
\[0, 1\] before thresholding (binary surfaces) or returning directly
(continuous surfaces).

The two autocorrelation ranges allow simulation of covariates that
operate at different spatial scales, a common scenario in landscape
ecology where some resources are patchily distributed at fine scales and
others vary broadly across the study area.

When `user_seed` is provided, independent but reproducible seeds are
derived for each surface as multiples of `user_seed`.
