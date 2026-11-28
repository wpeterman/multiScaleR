# Convert a sampling geometry to a footprint SD

**Experimental.** Returns the per-axis standard deviation, in map units,
of the area over which a response is measured.
[`multiScale_joint()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_joint.md)
uses it to relate the apparent scale of effect in a data set to the
ecological scale (`scales = "footprint"`).

## Usage

``` r
footprint_sd(
  type = c("point", "disk", "square", "rectangle", "halfnormal"),
  size = NULL,
  width = NULL
)
```

## Arguments

- type:

  Sampling geometry: `"point"` (no footprint), `"disk"` (counts or
  detections within a fixed radius), `"square"` (a square plot or grid
  cell), `"rectangle"` (a plot or transect strip), or `"halfnormal"`
  (animals detected with a half-normal detection function, as in
  distance sampling).

- size:

  Radius for `"disk"`, side length for `"square"`, length for
  `"rectangle"`, and the detection-function scale for `"halfnormal"`, in
  map units. Ignored for `"point"`.

- width:

  Width for `"rectangle"`, in map units.

## Value

A numeric vector of footprint SDs in map units.

## Details

The SD is the square root of the per-axis variance of a point drawn at
random from the footprint, weighted by detection where that applies:
`r / 2` for a disk of radius `r`, `L / sqrt(12)` for a square of side
`L`, `sqrt((L^2 + W^2) / 24)` for an `L` by `W` rectangle (its two axes
averaged, since the Gaussian kernel is isotropic), and the detection
scale itself for a half-normal detection function, whose detection
weights over the plane form a bivariate normal distribution.

## See also

[`multiScale_joint()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_joint.md).
[`vignette("v16-joint-data-sets", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v16-joint-data-sets.md)
uses a footprint for acoustic recorders.

## Examples

``` r
footprint_sd("disk", 100)          # 100 m fixed-radius point count: 50
#> [1] 50
footprint_sd("halfnormal", 60)     # distance sampling, sigma = 60 m: 60
#> [1] 60
footprint_sd("rectangle", 200, 20) # 200 m by 20 m transect strip
#> [1] 41.02845
```
