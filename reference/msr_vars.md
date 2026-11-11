# Define multiScaleR covariate transformations

These helpers define named model covariates as transformations of one or
more source raster layers. They are optional; if omitted,
[`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
preserves the default behavior where each raster layer becomes one
optimized kernel-weighted covariate with the same name.

Use `msr_vars()` to collect one or more `kernel_var()`,
`landscape_var()`, or `surface_var()` specifications into a single
object that is passed to the `scale_vars` argument of
[`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
and
[`kernel_scale.raster`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md).

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
```

## Arguments

- ...:

  Named `kernel_var()`, `landscape_var()`, or `surface_var()`
  specifications. Each argument must be named; the name becomes the
  covariate column name in the model data frame. All names must be
  unique and cannot be `"x"`, which is reserved for internal coordinate
  data. At least one specification is required.

- source:

  Character scalar. Name of the source raster layer in `raster_stack`
  from which the covariate is derived. Must exactly match a layer name
  in the raster stack provided to
  [`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md).

- metric:

  Character scalar. The metric to compute within the circular buffer.
  For `landscape_var()` this is a landscape ecology metric computed on a
  categorical raster; for `surface_var()` this is a continuous surface
  texture metric (`"sa"` or `"sq"`) computed on a continuous raster.
  Must be one of the supported metrics (see Details). Matched
  case-insensitively.

- radius:

  Optional positive numeric. Fixed buffer radius in the same units as
  the projection. When `NULL` (default), the radius is treated as a free
  parameter and optimized alongside the model. When supplied, the metric
  is computed at this fixed radius and no scale optimization is
  performed for this covariate.

- base:

  Positive numeric (not equal to 1). Logarithm base used by every
  entropy-based metric: the diversity metrics (`"shdi"`, `"shei"`,
  `"msidi"`, `"msiei"`) and the information-theory metrics (`"ent"`,
  `"condent"`, `"joinent"`, `"mutinf"`). Default: `exp(1)` (natural log,
  so entropies are in nats). Use `2` for bits or `10` for Hartley units.
  To reproduce the information-theory values reported by FRAGSTATS and
  the `landscapemetrics` package, set `base = 2`. The base cancels in
  the ratio metrics (`"relmutinf"`, `"iji"`), so they are unaffected by
  it.

- classes_max:

  Optional positive numeric. Maximum number of possible patch types in
  the landscape, used only for relative patch richness (`"rpr"`). If
  `NULL` (default), the observed number of classes within each buffer is
  used as the denominator.

- class:

  Optional integer-like scalar. The focal class for the class-level
  metrics (`"clumpy"`, `"pland"`, `"ca"`), which describe one category
  rather than the whole landscape. It must match one of the integer
  class codes in the source raster. It is required for those three
  metrics and is not accepted for any other metric (which describe the
  whole buffer). Where the focal class is absent from a buffer,
  `"pland"` and `"ca"` return 0 (zero cover) and `"clumpy"` returns
  `NA`.

- weighted:

  Logical, for `surface_var()` only. When `FALSE` (default), the surface
  metric is computed over a hard-edged circular neighborhood and the
  optimized scale parameter is that neighborhood radius. When `TRUE`,
  each cell is weighted by the model's distance kernel (the `kernel`
  passed to
  [`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md))
  and the optimized parameter is the kernel scale (sigma), exactly as
  for a kernel-weighted mean. The weighted form estimates the scale of
  effect of heterogeneity with a likelihood that is smooth in sigma. It
  is available only for the distribution metrics (`"sa"`, `"sq"`,
  `"ssk"`, `"sku"`), not the slope metric (`"sdq"`), and the scale is
  always optimized, so `radius` must be `NULL`.

## Value

- `msr_vars()`:

  A data frame of class `"multiScaleR_vars"` with one row per covariate.
  Columns are `covariate` (the name assigned in `...`), `source` (source
  raster layer name), `type` (`"kernel"`, `"landscape"`, or
  `"surface"`), `metric` (landscape or surface metric, or `NA`),
  `radius` (fixed radius or `NA` when optimized), `optimize` (logical
  indicating whether the scale is estimated), `weighted` (logical, the
  kernel-weighted form of a surface metric), `base` (logarithm base for
  entropy metrics), `classes_max` (richness denominator for `"rpr"`, or
  `NA`), and `class` (focal class for the class-level metrics, or `NA`).

- `kernel_var()`:

  An internal list of class `"multiScaleR_var"` representing a
  kernel-weighted mean covariate. Pass one or more of these inside
  `msr_vars()`.

- `landscape_var()`:

  An internal list of class `"multiScaleR_var"` representing a landscape
  metric covariate. Pass one or more of these inside `msr_vars()`.

- `surface_var()`:

  An internal list of class `"multiScaleR_var"` representing a
  continuous surface texture metric covariate. Pass one or more of these
  inside `msr_vars()`.

## Details

**Kernel vs. landscape covariates**

`kernel_var(source)` defines a kernel-weighted mean of the continuous
raster values within a circular neighborhood. The neighborhood radius
(sigma) is optimized by
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

`landscape_var(source, metric)` derives a landscape ecology metric from
a categorical (or thresholded) raster layer within a circular buffer.
Raster classes must be encoded as finite integer-like values. The buffer
radius can be fixed or optimized.

`surface_var(source, metric)` derives a continuous surface texture
metric from a continuous raster layer (elevation, NDVI, canopy height,
temperature, moisture, and similar surfaces) within a circular buffer.
These metrics summarize within-neighborhood heterogeneity rather than an
average value, so they complement `kernel_var()` (which summarizes the
mean). The buffer radius can be fixed or optimized. Note that these
covariates measure variability, not central tendency: changing the
radius changes both the ecological scale and the statistical quantity
being summarized.

By default a surface metric uses a hard-edged neighborhood. Setting
`weighted = TRUE` instead weights each cell by the model's distance
kernel and optimizes the kernel scale (sigma), which estimates the scale
of effect of heterogeneity with a smooth likelihood (see the `weighted`
argument). This applies to the distribution metrics (`"sa"`, `"sq"`,
`"ssk"`, `"sku"`).

**Supported landscape metrics**

Composition metrics (require a categorical raster):

- `"shdi"`:

  Shannon diversity index.

- `"shei"`:

  Shannon evenness index.

- `"sidi"`:

  Simpson diversity index.

- `"siei"`:

  Simpson evenness index. Can be uninformative for optimized scales when
  local evenness is uniform across sample points.

- `"msidi"`:

  Modified Simpson diversity index.

- `"msiei"`:

  Modified Simpson evenness index. Can be uninformative for optimized
  scales when local evenness is uniform across sample points.

- `"pr"`:

  Patch richness (number of distinct classes).

- `"prd"`:

  Patch richness density (pr per 100 ha).

- `"rpr"`:

  Relative patch richness (pr / classes_max \* 100).

- `"ta"`:

  Total area of the buffer (ha).

Edge metrics (require adjacency information; cell IDs are cached
internally):

- `"ed"`:

  Edge density (m/ha).

- `"te"`:

  Total edge length (m).

- `"lsi"`:

  Landscape shape index.

Adjacency metrics (configuration; built from cell-pair counts):

- `"ai"`:

  Aggregation index (%). Area-weighted percentage of like adjacencies
  relative to the maximum possible for each class. Higher values mean
  classes are more spatially aggregated.

- `"pladj"`:

  Proportion of like adjacencies (%). Percentage of all cell-pair
  boundaries that fall between cells of the same class.

- `"contag"`:

  Contagion index (%). Ranges from near 0 (maximum interspersion of
  classes) to 100 (a single class). Requires at least 2 classes within
  the buffer.

- `"iji"`:

  Interspersion and juxtaposition index (%). The evenness of the
  adjacencies among the different classes, ignoring like adjacencies.
  Ranges from near 0 (one class pair dominates the between-class
  boundaries) to 100 (every class is equally adjacent to every other).
  Undefined, and returned as `NA`, for fewer than 3 classes. For
  optimized scales, prescreen for usable variation because local windows
  with too few classes or between-class adjacencies can produce all `NA`
  values.

Information theory metrics (configuration and composition; built from
the class-pair co-occurrence matrix following Nowosad and Stepinski,
2019). The `base` argument sets the logarithm base; set `base = 2` for
bits, as in `landscapemetrics`:

- `"ent"`:

  Marginal entropy, H(x). The diversity of the cell values entering
  adjacencies; it increases with the number and evenness of classes.

- `"joinent"`:

  Joint entropy, H(x, y). The diversity of the class-pair adjacencies
  themselves; the total complexity of the pattern.

- `"condent"`:

  Conditional entropy, H(y \| x). The configurational complexity that
  remains once composition is accounted for: the uncertainty about a
  neighbor's class given the focal cell's class.

- `"mutinf"`:

  Mutual information, I(y, x). How much knowing a cell's class tells you
  about its neighbor; it rises as like classes cluster together. This is
  the configuration signal separated from composition.

- `"relmutinf"`:

  Relative mutual information, I(y, x) / H(x). Mutual information
  rescaled by marginal entropy to the range 0 to 1, so that
  configuration can be compared across landscapes of differing
  composition. It is 1 when mutual information is 0 (the
  `landscapemetrics` convention) and is unaffected by `base`.

Class-level metrics (describe a single focal class; require the `class`
argument):

- `"clumpy"`:

  Clumpiness index for the focal class. Ranges from -1 (maximally
  disaggregated) through 0 (randomly distributed) to 1 (maximally
  clumped), correcting the proportion of like adjacencies for the
  class's own abundance. Returned as `NA` where the class fills or is
  absent from the buffer.

- `"pland"`:

  Percentage of landscape (%). The focal class's share of the valid
  buffer cells.

- `"ca"`:

  Class area (hectares) of the focal class within the buffer.

Adjacency and information-theory projections are class-count limited to
avoid accidental use of continuous rasters and to keep FFT projection
times bounded. The current ceilings are 200 classes for the
like-adjacency metrics (`"ai"`, `"pladj"`, `"clumpy"`) and 50 classes
for the metrics built from the full class-pair matrix (`"contag"`,
`"iji"`, and the information-theory family), whose cost scales with the
square of the class count.

Before optimizing the scale of `"iji"`, `"siei"`, or `"msiei"`, inspect
a pilot
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
result or a few fixed radii. These metrics can be undefined or nearly
constant when buffers contain too few classes, too few between-class
adjacencies, or nearly identical evenness at all sample points. In that
case, use a larger radius, increase `max_D`, or choose a more stable
diversity, entropy, or edge metric.

**Supported surface texture metrics**

Surface metrics (require a continuous raster):

- `"sa"`:

  Average roughness: the mean absolute deviation of the neighborhood
  values from their mean.

- `"sq"`:

  Root mean square (RMS) roughness: the sample standard deviation of the
  neighborhood values. More sensitive to large deviations than `"sa"`.

- `"ssk"`:

  Skewness: the asymmetry of the neighborhood value distribution.
  Positive when occasional high values dominate, negative when
  occasional low values do.

- `"sku"`:

  Kurtosis (excess): the peakedness and tail weight of the neighborhood
  value distribution relative to a normal distribution (which has excess
  kurtosis 0).

- `"sdq"`:

  Root mean square slope: the RMS of the local gradient magnitude, a
  measure of how steeply the surface changes (terrain or structural
  ruggedness).

- `"sdr"`:

  Surface area ratio: the percentage by which the true 3D surface area
  exceeds the flat projected area. Behaves as a continuous edge-density
  analog (more relief gives a larger ratio).

The dispersion and shape metrics (`"sa"`, `"sq"`, `"ssk"`, `"sku"`)
match the `geodiv` package (`sq` uses an N - 1 denominator as
[`stats::sd()`](https://rdrr.io/r/stats/sd.html) does; `ssk` uses the
bias-adjusted skewness and `sku` the excess kurtosis, matching
`geodiv`'s defaults). `"sdq"` is reported as a true slope (the gradient
is divided by the cell resolution);
[`geodiv::sdq()`](https://rdrr.io/pkg/geodiv/man/sdq.html) instead uses
a cell spacing of one, so the two agree after multiplying by the
resolution. `"sdr"` is computed in physical units (real elevation values
and the real cell size), so a tilted plane of slope `s` returns
`(sqrt(1 + s^2) - 1) * 100`;
[`geodiv::sdr()`](https://rdrr.io/pkg/geodiv/man/sdr.html) instead
rescales the surface to a unit range and so is dimensionless. The
physical definition is the interpretable one for terrain and, unlike the
rescaled version, can be projected by FFT.

Raster projection of `"sq"`, `"ssk"`, `"sku"`, `"sdq"`, and `"sdr"` uses
fast Fourier transform (FFT) convolution; projection of `"sa"` uses an
exact compiled masked focal pass, because average roughness has no
closed-form FFT decomposition. `"sdq"` and `"sdr"` use each cell's
neighbors, so (like the landscape edge metrics) they cache raster cell
IDs internally, and their projection agrees with point extraction to
within a small boundary tolerance rather than bit-for-bit.

**Mixing covariate types**

Multiple covariate types can be combined in one `msr_vars()` call. For
example, a kernel-weighted mean of forest cover and a fixed-radius edge
density metric can both be defined and passed together to
[`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md):


    vars <- msr_vars(
      forest_mean  = kernel_var("forest"),
      forest_ed500 = landscape_var("forest", metric = "ed", radius = 500)
    )
    kernel_inputs <- kernel_prep(pts, rasters, max_D = 1000,
                                 scale_vars = vars)

Covariates with a fixed `radius` are computed once and not re-evaluated
during optimization, which can meaningfully reduce computation time.

## References

Nowosad, J., & Stepinski, T. F. (2019). Information theory as a
consistent framework for quantification and classification of landscape
patterns. *Landscape Ecology*, 34(9), 2091-2101.
[doi:10.1007/s10980-019-00830-x](https://doi.org/10.1007/s10980-019-00830-x)

## See also

[`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md),
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
[`kernel_scale.raster`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md).
For a worked walkthrough of the landscape metrics see
[`vignette("landscape_metric_covariates", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/landscape_metric_covariates.md),
and for the continuous-surface metrics see
[`vignette("surface_metric_covariates", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/surface_metric_covariates.md).

## Examples

``` r
## Kernel-weighted mean only (equivalent to default behavior)
vars <- msr_vars(
  forest_prop = kernel_var("forest")
)

## Combining kernel, landscape, and surface covariates
vars <- msr_vars(
  forest_prop = kernel_var("forest"),
  forest_ed   = landscape_var("forest", metric = "ed"),
  cover_shdi_500 = landscape_var("landcover", metric = "shdi", radius = 500),
  elev_sd     = surface_var("elevation", metric = "sq")
)

## landscape_var with natural-log Shannon diversity (the default)
landscape_var("landcover", metric = "shdi")
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

## landscape_var with log2 Shannon diversity and a fixed 250 m radius
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

## Configuration: interspersion (landscape-level), no class argument
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

## Information theory in bits (set base = 2 to match landscapemetrics)
landscape_var("landcover", metric = "mutinf", radius = 500, base = 2)
#> $type
#> [1] "landscape"
#> 
#> $source
#> [1] "landcover"
#> 
#> $metric
#> [1] "mutinf"
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

## Class-level metrics require `class`: clumpiness and cover of class 3
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

## surface_var: optimized RMS roughness, and fixed-radius average roughness
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

## Kernel-weighted roughness: optimize the scale of effect of heterogeneity
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

## Define a single optimized kernel-weighted mean covariate
kv <- kernel_var("forest")
```
