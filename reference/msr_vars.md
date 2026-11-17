# Define Model Covariates from Raster Layers

These helpers declare how each model covariate is computed from a source
raster layer, and whether its scale is estimated or fixed. Use
`kernel_var()` for a kernel-weighted mean, `landscape_var()` for a
landscape metric of a categorical raster, and `surface_var()` for a
texture metric of a continuous raster. Collect the specifications with
`msr_vars()` and pass the result to the `scale_vars` argument of
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
and
[`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md).

The helpers are optional. Without `scale_vars`,
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
turns each raster layer into one kernel-weighted mean covariate with the
layer's name.

## Usage

``` r
msr_vars(...)

kernel_var(source)

landscape_var(
  source,
  metric,
  radius = NULL,
  base = exp(1),
  classes_max = NULL,
  class = NULL
)

surface_var(source, metric, radius = NULL, weighted = FALSE)

# S3 method for class 'multiScaleR_vars'
print(x, ...)
```

## Arguments

- ...:

  For `msr_vars()`, one or more named specifications created with
  `kernel_var()`, `landscape_var()`, or `surface_var()`. Each name
  becomes the covariate's column name in the model data. Names must be
  unique, and `"x"` is reserved for internal coordinate data. At least
  one specification is required. For
  [`print()`](https://rdrr.io/r/base/print.html), further arguments
  passed to
  [`base::print.data.frame()`](https://rdrr.io/r/base/print.dataframe.html).

- source:

  Character scalar. Name of the source layer in the `raster_stack` given
  to
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md).
  It must match a layer name exactly. Several covariates can share one
  source.

- metric:

  Character scalar. Code of the metric to compute within the buffer,
  matched case-insensitively. `landscape_var()` accepts the 25 landscape
  metrics and `surface_var()` the six surface metrics defined in
  Details.

- radius:

  Optional positive numeric scalar. Buffer radius in map units. With the
  default `NULL`, the radius is estimated with the model. It is a hard
  radius, so it has no Hessian standard error (see
  [multiScaleR-kernels](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md)).
  When supplied, the metric is computed once at this radius and enters
  the model as an ordinary covariate, with no scale to estimate. Keep it
  within the `max_D` of
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md),
  because cells beyond `max_D` are never extracted.

- base:

  Positive numeric scalar, not equal to 1. Base of the logarithm used by
  Shannon diversity (`"shdi"`) and the information-theory metrics
  (`"ent"`, `"condent"`, `"joinent"`, `"mutinf"`). Default: `exp(1)`, so
  entropies are in nats. Use `2` for bits, which reproduces FRAGSTATS
  and the `landscapemetrics` package for the information-theory metrics,
  or `10` for Hartleys. The base cancels in the ratio metrics (`"shei"`,
  `"contag"`, `"iji"`, `"relmutinf"`). The modified Simpson metrics
  (`"msidi"`, `"msiei"`) always use the natural log, as FRAGSTATS and
  `landscapemetrics` define them. For these and for metrics that take no
  logarithm, `landscape_var()` warns and ignores a non-default `base`.

- classes_max:

  Positive numeric scalar. Number of classes that could occur in the
  landscape. Required for relative patch richness (`"rpr"`), where it is
  the denominator, and unused by every other metric. The observed class
  count is never substituted, so `landscape_var()` stops if `"rpr"` is
  requested without it.

- class:

  Optional integer-like scalar. The focal class code for the class-level
  metrics (`"clumpy"`, `"pland"`, `"ca"`), which describe one class
  rather than the whole buffer. It must match one of the class codes in
  the source raster. It is required for those three metrics and an error
  for every other metric.

- weighted:

  Logical, for `surface_var()` only. With `FALSE` (default), the metric
  summarizes every cell within a hard-edged circle, and the estimated
  scale is that radius. With `TRUE`, each cell is weighted by the
  model's distance kernel (the `kernel` given to
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md))
  and the estimated scale is the kernel's sigma, as for a
  kernel-weighted mean. Gaussian, negative exponential, and exponential
  power weights change smoothly with sigma, so sigma gets a standard
  error. With `kernel = "fixed"` the weights are 1 or 0, so the
  covariate is again a hard radius without a standard error. Weighting
  applies only to the distribution metrics (`"sa"`, `"sq"`, `"ssk"`,
  `"sku"`), not to `"sdq"` or `"sdr"`, and requires `radius = NULL`.

- x:

  A `multiScaleR_vars` object returned by `msr_vars()`.

## Value

`msr_vars()` returns a data frame of class
`c("multiScaleR_vars", "data.frame")` with one row per covariate and
these columns:

- `covariate`: the name given in `...`.

- `source`: the source raster layer.

- `type`: `"kernel"`, `"landscape"`, or `"surface"`.

- `metric`: the metric code, or `NA` for a kernel mean.

- `radius`: the fixed radius in map units, or `NA` when the scale is
  estimated.

- `optimize`: `TRUE` when the scale is estimated with the model.

- `weighted`: `TRUE` for a kernel-weighted surface metric.

- `base`: the logarithm base actually used.

- `classes_max`: the denominator of `"rpr"`, or `NA`.

- `class`: the focal class of a class-level metric, or `NA`.

`kernel_var()`, `landscape_var()`, and `surface_var()` each return a
list of class `"multiScaleR_var"` holding one specification. Pass it to
`msr_vars()`. It has no use on its own.

[`print()`](https://rdrr.io/r/base/print.html) displays the
specification table and returns `x` invisibly.

## Choose a covariate type

|  |  |  |  |
|----|----|----|----|
| Helper | Source raster | What it summarizes | Scale |
| `kernel_var()` | Any numeric layer | Kernel-weighted mean of the cell values | Sigma, estimated, with a standard error |
| `landscape_var()` | Categorical, with integer class codes | Composition, edges, aggregation, information theory, or one focal class | Hard radius, estimated without a standard error, or fixed |
| `surface_var()` | Continuous, such as elevation or canopy height | Variation, distribution shape, slope, or relief | Hard radius (estimated or fixed), or sigma with `weighted = TRUE` |

A kernel mean describes how much of something surrounds a point.
Landscape and surface metrics describe how it is arranged or how much it
varies, so they complement the mean rather than replace it. Note that
for a metric, changing the radius changes both the spatial scale and the
statistical quantity summarized. Landscape metrics require finite,
integer-like class codes, and both metric types require square raster
cells.

Whether a scale is a smooth kernel width or a hard radius decides how it
can be estimated and reported.
[multiScaleR-kernels](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md)
explains the difference, and why a hard radius has a stepped likelihood
and no Hessian standard error.

## Combine covariate types

One `msr_vars()` call can mix all three types, including several
covariates from the same source layer. For example, a kernel-weighted
mean of forest cover and a fixed-radius edge density of the same layer:

    vars <- msr_vars(
      forest_mean  = kernel_var("forest"),
      forest_ed500 = landscape_var("forest", metric = "ed", radius = 500)
    )
    kernel_inputs <- kernel_prep(pts, rasters, max_D = 1000,
                                 scale_vars = vars)

A covariate with a fixed `radius` is computed once and is not
re-evaluated during optimization, which can save considerable time.
Under `kernel = "expow"`, a hard-radius covariate does not use the
kernel shape. Its reported shape and shape standard error are `NA`, and
its unused shape position is left out of model-selection parameter
counts.

Not every fitting route accepts every covariate type. Kernel means,
weighted `"sq"`, and fixed-radius metrics work with every route,
including the analytic backends of
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
and
[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md),
within each route's kernel restrictions. Those routes need a Gaussian or
negative exponential kernel for weighted `"sq"`. An estimated hard
radius, and any weighted metric other than `"sq"`, needs the default
refit backend. See
[multiScaleR-fitting-methods](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-fitting-methods.md)
and
[`vignette("v07-fitting-methods", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v07-fitting-methods.md)
for the full feature table.

## Landscape metric definitions

Each landscape metric summarizes the cells whose centers lie within the
radius. In the formulas, `n` is the number of valid (non-missing) cells,
`n_i` the number in class `i`, `p_i = n_i / n`, `m` the number of
classes present, `res` the cell size in map units, and
`A = n * res^2 / 10000` the buffer area (hectares when map units are
meters). `g_ij` counts pairs of cells that share a side, have classes
`i` and `j`, and both lie in the buffer. Each pair is counted in both
directions, and `q_ij = g_ij / sum(g)`. Logarithms use `base` where
marked, and the natural log otherwise.

**Composition** (how many classes, and how evenly they share the
buffer):

|  |  |  |  |
|----|----|----|----|
| Code | Name | Definition | Units or range |
| `shdi` | Shannon diversity | `-sum(p_i * log(p_i))`, with `base` | Log units, 0 or more |
| `shei` | Shannon evenness | `shdi / log(m)`, with `base`. 0 when `m = 1` | 0 to 1 |
| `sidi` | Simpson diversity | `1 - sum(p_i^2)` | 0 to 1 |
| `siei` | Simpson evenness | `sidi / (1 - 1 / m)`. `NA` when `m = 1` | 0 to 1 |
| `msidi` | Modified Simpson diversity | `-log(sum(p_i^2))` | Natural log units, 0 or more |
| `msiei` | Modified Simpson evenness | `msidi / log(m)`. Undefined when `m = 1` | 0 to 1 |
| `pr` | Patch richness | `m` | Number of classes |
| `prd` | Patch richness density | `m / A * 100` | Classes per 100 ha |
| `rpr` | Relative patch richness | `m / classes_max * 100` | Percent |
| `ta` | Total area | `A` | Hectares |

**Edge** (how much boundary separates different classes). Here `E` is
the number of cell sides shared by two buffer cells of different
classes, and `B` the number of cell sides on the buffer's outer boundary
or next to a missing cell:

|  |  |  |  |
|----|----|----|----|
| Code | Name | Definition | Units or range |
| `te` | Total edge | `E * res` | Map units (meters) |
| `ed` | Edge density | `E * res / A` | Meters per hectare |
| `lsi` | Landscape shape index | `(E + B) / Pmin(n)`, where `Pmin(n)` is the smallest possible perimeter, in cell sides, of `n` cells | Ratio, 1 or more |

**Aggregation** (whether cells of a class sit together or are mixed with
other classes):

|  |  |  |  |
|----|----|----|----|
| Code | Name | Definition | Units or range |
| `pladj` | Proportion of like adjacencies | `sum(g_ii) / sum(g) * 100` | Percent |
| `ai` | Aggregation index | `sum(p_i * (g_ii / 2) / gmax_i) * 100`, where `gmax_i` is the most like adjacencies `n_i` cells can have | Percent |
| `contag` | Contagion | `(1 + sum(q_ij * log(q_ij)) / (2 * log(m))) * 100`. `NA` when `m < 2` | 0 (interspersed) to 100 (one class) |
| `iji` | Interspersion and juxtaposition | `-sum(e * log(e)) / log(m * (m - 1) / 2) * 100`, where `e` is each class pair's share of the adjacencies between different classes. `NA` when `m < 3` or no such adjacencies exist | Near 0 (one pair dominates) to 100 (all pairs equally adjacent) |

**Information theory** (pattern complexity, split into composition and
configuration, following Nowosad and Stepinski 2019). With
`q_i = sum_j(q_ij)`:

|  |  |  |  |
|----|----|----|----|
| Code | Name | Definition | Units or range |
| `ent` | Marginal entropy, H(x) | `-sum(q_i * log(q_i))`, with `base` | Log units. Rises with the number and evenness of classes |
| `joinent` | Joint entropy, H(x, y) | `-sum(q_ij * log(q_ij))`, with `base` | Log units. Total complexity of the pattern |
| `condent` | Conditional entropy, H(y given x) | `joinent - ent` | Log units. Uncertainty about a neighbor's class given a cell's class |
| `mutinf` | Mutual information, I(y, x) | `ent - condent` | Log units. How much a cell's class tells you about its neighbor's |
| `relmutinf` | Relative mutual information | `mutinf / ent`, and 1 when `mutinf = 0` | 0 to 1, unaffected by `base` |

Both clumping and regular alternation make neighbors predictable and
raise `mutinf`. Its size also depends on `ent`, so it is not a
composition-free aggregation index. `relmutinf` rescales it so that
configuration can be compared across buffers with different composition.

**Class level** (one focal class `k`, set with `class`):

|  |  |  |  |
|----|----|----|----|
| Code | Name | Definition | Units or range |
| `pland` | Percentage of landscape | `n_k / n * 100`. 0 when the class is absent | Percent |
| `ca` | Class area | `n_k * res^2 / 10000`. 0 when the class is absent | Hectares |
| `clumpy` | Clumpiness index | With `G = g_kk / (4 * n_k - Pmin(n_k))` and `P = n_k / n`: `(G - P) / (1 - P)` when `G >= P` or `P >= 0.5`, otherwise `(G - P) / P`, the FRAGSTATS definition. `NA` when the class is absent or fills the buffer | Unitless, from -1 (maximally dispersed) through 0 (random) to 1 (maximally clumped). `landscapemetrics` 2.2.1 divides by `-P` in the second case, which reverses the sign of dispersed classes |

## Surface metric definitions

Each surface metric summarizes the values `z` of the `n` valid cells
whose centers lie within the radius. `s` is the sample standard
deviation (`"sq"`) and `res` the cell size in map units.

|  |  |  |  |
|----|----|----|----|
| Code | Name | Definition | Units |
| `sa` | Average roughness | `sum(abs(z - mean(z))) / n` | Raster units |
| `sq` | Root mean square roughness | `sqrt(sum((z - mean(z))^2) / (n - 1))`. `NA` when `n < 2` | Raster units |
| `ssk` | Skewness | `sqrt(n * (n - 1)) / (n - 2) * (sum((z - mean(z))^3) / n) / s^3`. `NA` when `n < 3` or `s = 0` | Unitless. Positive when a few high values stand out |
| `sku` | Excess kurtosis | `(sum((z - mean(z))^4) / n) / s^4 - 3`. `NA` when `s = 0` | Unitless. 0 for normally distributed values |
| `sdq` | Root mean square slope | Square root of the mean of `((z - z_right)^2 + (z - z_below)^2) / res^2` over cells whose right and lower neighbors are both in the buffer | Raster units per map unit |
| `sdr` | Surface area ratio | `(sum(a) / (k * res^2) - 1) * 100`, where `a` is the three-dimensional area of each of the `k` blocks of 2 x 2 cell centers in the buffer | Percent |

`"sa"`, `"sq"`, `"ssk"`, and `"sku"` match the `geodiv` package. `"sdq"`
is a true slope (rise over run).
[`geodiv::sdq()`](https://rdrr.io/pkg/geodiv/man/sdq.html) assumes a
cell spacing of 1, so it equals `"sdq"` multiplied by `res`. `"sdr"`
uses the real values and cell size, so a plane tilted at slope `t`
returns `(sqrt(1 + t^2) - 1) * 100`.
[`geodiv::sdr()`](https://rdrr.io/pkg/geodiv/man/sdr.html) rescales the
surface to a unit range first, so its values differ. Read `"sdq"` and
`"sdr"` as a slope and an area ratio only when the raster values and the
map share units, as elevation in meters does.

## Kernel-weighted surface metrics

With `weighted = TRUE`, each cell gets the kernel weight `w = w(d)` of
its distance `d` from the point (see
[multiScaleR-kernels](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md)),
using the cells extracted within `max_D`. With the weighted mean
`mu = sum(w * z) / sum(w)`:

|       |                                                                 |
|-------|-----------------------------------------------------------------|
| Code  | Weighted definition                                             |
| `sa`  | `sum(w * abs(z - mu)) / sum(w)`                                 |
| `sq`  | `sqrt(sum(w * (z - mu)^2) / sum(w))`, and 0 for constant values |
| `ssk` | `(sum(w * (z - mu)^3) / sum(w)) / sq^3`. `NA` when `sq = 0`     |
| `sku` | `(sum(w * (z - mu)^4) / sum(w)) / sq^4 - 3`. `NA` when `sq = 0` |

Weighted `"sq"` divides by the sum of weights, so it is a
population-style spread and differs from the `n - 1` form of the
unweighted metric. The weighted skewness and kurtosis use plain
standardized moments, without the small-sample adjustment of the
unweighted `"ssk"`. Weighted `"sq"` is smooth in sigma except where the
weighted spread is zero. Absolute deviations (`"sa"`) and the
standardized higher moments have additional points where they are not
smooth or not defined. Within one evaluation, weighted metrics reuse
weights and centered moments only when their source, scale, shape, and
valid cells match exactly. Their scales are still estimated separately.
Under `kernel = "fixed"`, the fitted scale is a hard cutoff radius. It
has no standard error or Wald limits, although the weighted definitions
above still apply.

## Check a metric before you optimize it

Before estimating the radius of `"iji"`, `"siei"`, or `"msiei"`, inspect
a pilot
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
result or a few fixed radii. These metrics can be undefined or nearly
constant when buffers hold too few classes, too few adjacencies between
different classes, or nearly identical evenness at every point. `"iji"`
in particular returns `NA` with fewer than three classes, and small
buffers can make it `NA` everywhere. If so, use a larger radius,
increase `max_D`, or choose a more stable diversity, entropy, or edge
metric. A candidate scale that leaves a covariate undefined at a modeled
observation is rejected during optimization, so every candidate is
compared on the same observations.

The adjacency and information-theory metrics accept a limited number of
classes, which guards against passing a continuous raster by mistake and
keeps raster projection times bounded. The limits are 200 classes for
the like-adjacency metrics (`"ai"`, `"pladj"`, `"clumpy"`) and 50 for
the metrics built from the full class-pair table (`"contag"`, `"iji"`,
and the information-theory metrics), whose cost grows with the square of
the class count.

During optimization, the adjacency metrics reuse sorted cell and
neighbor-pair distances instead of rescanning every buffered cell at
each radius, and metrics on the same source and radius share their
counts. This saves time at the cost of retained memory, which
[`estimate_multiscale_ram()`](https://wpeterman.github.io/multiScaleR/reference/estimate_multiscale_ram.md)
includes. It does not change radius boundaries or metric definitions,
and it does not make a hard radius smooth or remove concerns about local
optima.

## Map a metric

[`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md)
computes the same covariates for every raster cell. Projection of
`"sq"`, `"ssk"`, `"sku"`, `"sdq"`, `"sdr"`, and the landscape metrics
uses fast Fourier transform (FFT) convolution by default. Unweighted
`"sa"` uses an exact compiled moving-window pass instead, because
average roughness cannot be written as a convolution. `"sdq"` and
`"sdr"` use each cell's neighbors, so their projected values agree with
point extraction to within a small tolerance at the buffer boundary
rather than exactly. Set `fft = FALSE` in
[`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md)
to use direct moving-window calculations instead.

## References

Nowosad, J., & Stepinski, T. F. (2019). Information theory as a
consistent framework for quantification and classification of landscape
patterns. *Landscape Ecology*, 34(9), 2091-2101.
[doi:10.1007/s10980-019-00830-x](https://doi.org/10.1007/s10980-019-00830-x)

## See also

[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md),
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
[`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md),
and
[multiScaleR-kernels](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md).
For worked examples, see
[`vignette("v04-landscape-metrics", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v04-landscape-metrics.md)
for the landscape metrics and
[`vignette("v05-surface-metrics", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v05-surface-metrics.md)
for the surface metrics.

## Examples

``` r
## Kernel-weighted mean only (the same as omitting scale_vars)
vars <- msr_vars(forest_prop = kernel_var("forest"))
vars
#> multiScaleR variable specifications:
#>    covariate source   type metric radius optimize weighted     base classes_max
#>  forest_prop forest kernel   <NA>     NA     TRUE    FALSE 2.718282          NA
#>  class
#>     NA

## Combine kernel, landscape, and surface covariates. The printed table
## shows which scales are estimated (optimize = TRUE) and which are fixed.
vars <- msr_vars(
  forest_prop    = kernel_var("forest"),
  forest_ed      = landscape_var("forest", metric = "ed"),
  cover_shdi_500 = landscape_var("landcover", metric = "shdi", radius = 500),
  elev_sd        = surface_var("elevation", metric = "sq")
)
print(vars)
#> multiScaleR variable specifications:
#>       covariate    source      type metric radius optimize weighted     base
#>     forest_prop    forest    kernel   <NA>     NA     TRUE    FALSE 2.718282
#>       forest_ed    forest landscape     ed     NA     TRUE    FALSE 2.718282
#>  cover_shdi_500 landcover landscape   shdi    500    FALSE    FALSE 2.718282
#>         elev_sd elevation   surface     sq     NA     TRUE    FALSE 2.718282
#>  classes_max class
#>           NA    NA
#>           NA    NA
#>           NA    NA
#>           NA    NA

## Shannon diversity in bits within a fixed 250 m radius
landscape_var("landcover", metric = "shdi", radius = 250, base = 2)
#> $type
#> [1] "landscape"
#> 
#> $source
#> [1] "landcover"
#> 
#> $metric
#> [1] "shdi"
#> 
#> $radius
#> [1] 250
#> 
#> $optimize
#> [1] FALSE
#> 
#> $weighted
#> [1] FALSE
#> 
#> $base
#> [1] 2
#> 
#> $classes_max
#> [1] NA
#> 
#> $class
#> [1] NA
#> 
#> attr(,"class")
#> [1] "multiScaleR_var"

## Interspersion of the whole buffer needs no focal class
landscape_var("landcover", metric = "iji", radius = 500)
#> $type
#> [1] "landscape"
#> 
#> $source
#> [1] "landcover"
#> 
#> $metric
#> [1] "iji"
#> 
#> $radius
#> [1] 500
#> 
#> $optimize
#> [1] FALSE
#> 
#> $weighted
#> [1] FALSE
#> 
#> $base
#> [1] 2.718282
#> 
#> $classes_max
#> [1] NA
#> 
#> $class
#> [1] NA
#> 
#> attr(,"class")
#> [1] "multiScaleR_var"

## Class-level metrics need `class`, here clumpiness and cover of class 3
landscape_var("landcover", metric = "clumpy", radius = 500, class = 3)
#> $type
#> [1] "landscape"
#> 
#> $source
#> [1] "landcover"
#> 
#> $metric
#> [1] "clumpy"
#> 
#> $radius
#> [1] 500
#> 
#> $optimize
#> [1] FALSE
#> 
#> $weighted
#> [1] FALSE
#> 
#> $base
#> [1] 2.718282
#> 
#> $classes_max
#> [1] NA
#> 
#> $class
#> [1] 3
#> 
#> attr(,"class")
#> [1] "multiScaleR_var"
landscape_var("landcover", metric = "pland", radius = 500, class = 3)
#> $type
#> [1] "landscape"
#> 
#> $source
#> [1] "landcover"
#> 
#> $metric
#> [1] "pland"
#> 
#> $radius
#> [1] 500
#> 
#> $optimize
#> [1] FALSE
#> 
#> $weighted
#> [1] FALSE
#> 
#> $base
#> [1] 2.718282
#> 
#> $classes_max
#> [1] NA
#> 
#> $class
#> [1] 3
#> 
#> attr(,"class")
#> [1] "multiScaleR_var"

## Surface metrics with an estimated hard radius, a fixed radius, and
## kernel weights (sigma estimated)
surface_var("elevation", metric = "sq")
#> $type
#> [1] "surface"
#> 
#> $source
#> [1] "elevation"
#> 
#> $metric
#> [1] "sq"
#> 
#> $radius
#> [1] NA
#> 
#> $optimize
#> [1] TRUE
#> 
#> $weighted
#> [1] FALSE
#> 
#> $base
#> [1] 2.718282
#> 
#> $classes_max
#> [1] NA
#> 
#> $class
#> [1] NA
#> 
#> attr(,"class")
#> [1] "multiScaleR_var"
surface_var("ndvi", metric = "sa", radius = 300)
#> $type
#> [1] "surface"
#> 
#> $source
#> [1] "ndvi"
#> 
#> $metric
#> [1] "sa"
#> 
#> $radius
#> [1] 300
#> 
#> $optimize
#> [1] FALSE
#> 
#> $weighted
#> [1] FALSE
#> 
#> $base
#> [1] 2.718282
#> 
#> $classes_max
#> [1] NA
#> 
#> $class
#> [1] NA
#> 
#> attr(,"class")
#> [1] "multiScaleR_var"
surface_var("elevation", metric = "sq", weighted = TRUE)
#> $type
#> [1] "surface"
#> 
#> $source
#> [1] "elevation"
#> 
#> $metric
#> [1] "sq"
#> 
#> $radius
#> [1] NA
#> 
#> $optimize
#> [1] TRUE
#> 
#> $weighted
#> [1] TRUE
#> 
#> $base
#> [1] 2.718282
#> 
#> $classes_max
#> [1] NA
#> 
#> $class
#> [1] NA
#> 
#> attr(,"class")
#> [1] "multiScaleR_var"

## Compute covariates at survey points from the bundled landscape
land <- terra::rast(system.file("extdata", "landscape.tif",
                                package = "multiScaleR"))
# land3 is continuous. Build a three-class layer from it
land$cover <- terra::classify(land$land3,
                              rbind(c(-Inf, 0.45, 1),
                                    c(0.45, 0.55, 2),
                                    c(0.55, Inf, 3)))
data("surv_pts")
vars <- msr_vars(
  land1_prop   = kernel_var("land1"),
  cover_shdi   = landscape_var("cover", metric = "shdi", radius = 400),
  land3_rough  = surface_var("land3", metric = "sq", radius = 300)
)
inputs <- kernel_prep(pts = terra::vect(surv_pts)[1:20, ],
                      raster_stack = land, max_D = 500,
                      scale_vars = vars, verbose = FALSE)
# Standardized covariates, one column per declared covariate
head(inputs$kernel_dat)
#>   land1_prop cover_shdi land3_rough
#> 1  1.1809859  0.5977263  -0.9685262
#> 2 -0.2849056 -0.8066326   0.2836450
#> 3 -1.0591167 -0.7715155  -1.0619046
#> 4 -0.2894405  0.8180410  -1.3727348
#> 5  1.0626868  0.7701467   0.4892986
#> 6 -0.4052110  0.1007936   1.0414276
```
