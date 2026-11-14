#' Define Model Covariates from Raster Layers
#'
#' @description
#' These helpers declare how each model covariate is computed from a source
#' raster layer, and whether its scale is estimated or fixed. Use
#' `kernel_var()` for a kernel-weighted mean, `landscape_var()` for a landscape
#' metric of a categorical raster, and `surface_var()` for a texture metric of
#' a continuous raster. Collect the specifications with `msr_vars()` and pass
#' the result to the `scale_vars` argument of [kernel_prep()] and
#' [kernel_scale.raster()].
#'
#' The helpers are optional. Without `scale_vars`, [kernel_prep()] turns each
#' raster layer into one kernel-weighted mean covariate with the layer's name.
#'
#' @param ... For `msr_vars()`, one or more named specifications created with
#'   `kernel_var()`, `landscape_var()`, or `surface_var()`. Each name becomes
#'   the covariate's column name in the model data. Names must be unique, and
#'   `"x"` is reserved for internal coordinate data. At least one
#'   specification is required. For `print()`, further arguments passed to
#'   [base::print.data.frame()].
#' @param source Character scalar. Name of the source layer in the
#'   `raster_stack` given to [kernel_prep()]. It must match a layer name
#'   exactly. Several covariates can share one source.
#' @param metric Character scalar. Code of the metric to compute within the
#'   buffer, matched case-insensitively. `landscape_var()` accepts the 25
#'   landscape metrics and `surface_var()` the six surface metrics defined in
#'   Details.
#' @param radius Optional positive numeric scalar. Buffer radius in map units.
#'   With the default `NULL`, the radius is estimated with the model. It is a
#'   hard radius, so it has no Hessian standard error (see
#'   [multiScaleR-kernels]). When supplied, the metric is computed once at
#'   this radius and enters the model as an ordinary covariate, with no scale
#'   to estimate. Keep it within the `max_D` of [kernel_prep()], because cells
#'   beyond `max_D` are never extracted.
#' @param weighted Logical, for `surface_var()` only. With `FALSE` (default),
#'   the metric summarizes every cell within a hard-edged circle, and the
#'   estimated scale is that radius. With `TRUE`, each cell is weighted by the
#'   model's distance kernel (the `kernel` given to [kernel_prep()]) and the
#'   estimated scale is the kernel's sigma, as for a kernel-weighted mean.
#'   Gaussian, negative exponential, and exponential power weights change
#'   smoothly with sigma, so sigma gets a standard error. With
#'   `kernel = "fixed"` the weights are 1 or 0, so the covariate is again a
#'   hard radius without a standard error. Weighting applies only to the
#'   distribution metrics (`"sa"`, `"sq"`, `"ssk"`, `"sku"`), not to `"sdq"`
#'   or `"sdr"`, and requires `radius = NULL`.
#' @param base Positive numeric scalar, not equal to 1. Base of the logarithm
#'   used by Shannon diversity (`"shdi"`) and the information-theory metrics
#'   (`"ent"`, `"condent"`, `"joinent"`, `"mutinf"`). Default: `exp(1)`, so
#'   entropies are in nats. Use `2` for bits, which reproduces FRAGSTATS and
#'   the `landscapemetrics` package for the information-theory metrics, or
#'   `10` for Hartleys. The base cancels in the ratio metrics (`"shei"`,
#'   `"contag"`, `"iji"`, `"relmutinf"`). The modified Simpson metrics
#'   (`"msidi"`, `"msiei"`) always use the natural log, as FRAGSTATS and
#'   `landscapemetrics` define them. For these and for metrics that take no
#'   logarithm, `landscape_var()` warns and ignores a non-default `base`.
#' @param classes_max Positive numeric scalar. Number of classes that could
#'   occur in the landscape. Required for relative patch richness (`"rpr"`),
#'   where it is the denominator, and unused by every other metric. The
#'   observed class count is never substituted, so `landscape_var()` stops if
#'   `"rpr"` is requested without it.
#' @param class Optional integer-like scalar. The focal class code for the
#'   class-level metrics (`"clumpy"`, `"pland"`, `"ca"`), which describe one
#'   class rather than the whole buffer. It must match one of the class codes
#'   in the source raster. It is required for those three metrics and an
#'   error for every other metric.
#' @param x A `multiScaleR_vars` object returned by `msr_vars()`.
#'
#' @return
#' `msr_vars()` returns a data frame of class
#' `c("multiScaleR_vars", "data.frame")` with one row per covariate and these
#' columns:
#' * `covariate`: the name given in `...`.
#' * `source`: the source raster layer.
#' * `type`: `"kernel"`, `"landscape"`, or `"surface"`.
#' * `metric`: the metric code, or `NA` for a kernel mean.
#' * `radius`: the fixed radius in map units, or `NA` when the scale is
#'   estimated.
#' * `optimize`: `TRUE` when the scale is estimated with the model.
#' * `weighted`: `TRUE` for a kernel-weighted surface metric.
#' * `base`: the logarithm base actually used.
#' * `classes_max`: the denominator of `"rpr"`, or `NA`.
#' * `class`: the focal class of a class-level metric, or `NA`.
#'
#' `kernel_var()`, `landscape_var()`, and `surface_var()` each return a list of
#' class `"multiScaleR_var"` holding one specification. Pass it to
#' `msr_vars()`. It has no use on its own.
#'
#' `print()` displays the specification table and returns `x` invisibly.
#'
#' @details
#' # Choose a covariate type
#'
#' | Helper | Source raster | What it summarizes | Scale |
#' |:-------|:--------------|:-------------------|:------|
#' | `kernel_var()` | Any numeric layer | Kernel-weighted mean of the cell values | Sigma, estimated, with a standard error |
#' | `landscape_var()` | Categorical, with integer class codes | Composition, edges, aggregation, information theory, or one focal class | Hard radius, estimated without a standard error, or fixed |
#' | `surface_var()` | Continuous, such as elevation or canopy height | Variation, distribution shape, slope, or relief | Hard radius (estimated or fixed), or sigma with `weighted = TRUE` |
#'
#' A kernel mean describes how much of something surrounds a point. Landscape
#' and surface metrics describe how it is arranged or how much it varies, so
#' they complement the mean rather than replace it. Note that for a metric,
#' changing the radius changes both the spatial scale and the statistical
#' quantity summarized. Landscape metrics require finite, integer-like class
#' codes, and both metric types require square raster cells.
#'
#' Whether a scale is a smooth kernel width or a hard radius decides how it
#' can be estimated and reported. [multiScaleR-kernels] explains the
#' difference, and why a hard radius has a stepped likelihood and no Hessian
#' standard error.
#'
#' # Combine covariate types
#'
#' One `msr_vars()` call can mix all three types, including several covariates
#' from the same source layer. For example, a kernel-weighted mean of forest
#' cover and a fixed-radius edge density of the same layer:
#'
#' ```
#' vars <- msr_vars(
#'   forest_mean  = kernel_var("forest"),
#'   forest_ed500 = landscape_var("forest", metric = "ed", radius = 500)
#' )
#' kernel_inputs <- kernel_prep(pts, rasters, max_D = 1000,
#'                              scale_vars = vars)
#' ```
#'
#' A covariate with a fixed `radius` is computed once and is not re-evaluated
#' during optimization, which can save considerable time. Under
#' `kernel = "expow"`, a hard-radius covariate does not use the kernel shape.
#' Its reported shape and shape standard error are `NA`, and its unused shape
#' position is left out of model-selection parameter counts.
#'
#' Not every fitting route accepts every covariate type. Kernel means,
#' weighted `"sq"`, and fixed-radius metrics work with every route, including
#' the analytic backends of [multiScale_optim()] and [multiScale_bayes()],
#' within each route's kernel restrictions. Those routes need a Gaussian or
#' negative exponential kernel for weighted `"sq"`. An estimated hard radius,
#' and any weighted metric other than `"sq"`, needs the default refit
#' backend. See [multiScaleR-fitting-methods] and
#' `vignette("v07-fitting-methods", package = "multiScaleR")` for the full
#' feature table.
#'
#' # Landscape metric definitions
#'
#' Each landscape metric summarizes the cells whose centers lie within the
#' radius. In the formulas, `n` is the number of valid (non-missing) cells,
#' `n_i` the number in class `i`, `p_i = n_i / n`, `m` the number of classes
#' present, `res` the cell size in map units, and `A = n * res^2 / 10000` the
#' buffer area (hectares when map units are meters). `g_ij` counts pairs of
#' cells that share a side, have classes `i` and `j`, and both lie in the
#' buffer. Each pair is counted in both directions, and `q_ij = g_ij / sum(g)`.
#' Logarithms use `base` where marked, and the natural log otherwise.
#'
#' **Composition** (how many classes, and how evenly they share the buffer):
#'
#' | Code | Name | Definition | Units or range |
#' |:-----|:-----|:-----------|:---------------|
#' | `shdi` | Shannon diversity | `-sum(p_i * log(p_i))`, with `base` | Log units, 0 or more |
#' | `shei` | Shannon evenness | `shdi / log(m)`, with `base`. 0 when `m = 1` | 0 to 1 |
#' | `sidi` | Simpson diversity | `1 - sum(p_i^2)` | 0 to 1 |
#' | `siei` | Simpson evenness | `sidi / (1 - 1 / m)`. `NA` when `m = 1` | 0 to 1 |
#' | `msidi` | Modified Simpson diversity | `-log(sum(p_i^2))` | Natural log units, 0 or more |
#' | `msiei` | Modified Simpson evenness | `msidi / log(m)`. Undefined when `m = 1` | 0 to 1 |
#' | `pr` | Patch richness | `m` | Number of classes |
#' | `prd` | Patch richness density | `m / A * 100` | Classes per 100 ha |
#' | `rpr` | Relative patch richness | `m / classes_max * 100` | Percent |
#' | `ta` | Total area | `A` | Hectares |
#'
#' **Edge** (how much boundary separates different classes). Here `E` is the
#' number of cell sides shared by two buffer cells of different classes, and
#' `B` the number of cell sides on the buffer's outer boundary or next to a
#' missing cell:
#'
#' | Code | Name | Definition | Units or range |
#' |:-----|:-----|:-----------|:---------------|
#' | `te` | Total edge | `E * res` | Map units (meters) |
#' | `ed` | Edge density | `E * res / A` | Meters per hectare |
#' | `lsi` | Landscape shape index | `(E + B) / Pmin(n)`, where `Pmin(n)` is the smallest possible perimeter, in cell sides, of `n` cells | Ratio, 1 or more |
#'
#' **Aggregation** (whether cells of a class sit together or are mixed with
#' other classes):
#'
#' | Code | Name | Definition | Units or range |
#' |:-----|:-----|:-----------|:---------------|
#' | `pladj` | Proportion of like adjacencies | `sum(g_ii) / sum(g) * 100` | Percent |
#' | `ai` | Aggregation index | `sum(p_i * (g_ii / 2) / gmax_i) * 100`, where `gmax_i` is the most like adjacencies `n_i` cells can have | Percent |
#' | `contag` | Contagion | `(1 + sum(q_ij * log(q_ij)) / (2 * log(m))) * 100`. `NA` when `m < 2` | 0 (interspersed) to 100 (one class) |
#' | `iji` | Interspersion and juxtaposition | `-sum(e * log(e)) / log(m * (m - 1) / 2) * 100`, where `e` is each class pair's share of the adjacencies between different classes. `NA` when `m < 3` or no such adjacencies exist | Near 0 (one pair dominates) to 100 (all pairs equally adjacent) |
#'
#' **Information theory** (pattern complexity, split into composition and
#' configuration, following Nowosad and Stepinski 2019). With
#' `q_i = sum_j(q_ij)`:
#'
#' | Code | Name | Definition | Units or range |
#' |:-----|:-----|:-----------|:---------------|
#' | `ent` | Marginal entropy, H(x) | `-sum(q_i * log(q_i))`, with `base` | Log units. Rises with the number and evenness of classes |
#' | `joinent` | Joint entropy, H(x, y) | `-sum(q_ij * log(q_ij))`, with `base` | Log units. Total complexity of the pattern |
#' | `condent` | Conditional entropy, H(y given x) | `joinent - ent` | Log units. Uncertainty about a neighbor's class given a cell's class |
#' | `mutinf` | Mutual information, I(y, x) | `ent - condent` | Log units. How much a cell's class tells you about its neighbor's |
#' | `relmutinf` | Relative mutual information | `mutinf / ent`, and 1 when `mutinf = 0` | 0 to 1, unaffected by `base` |
#'
#' Both clumping and regular alternation make neighbors predictable and raise
#' `mutinf`. Its size also depends on `ent`, so it is not a
#' composition-free aggregation index. `relmutinf` rescales it so that
#' configuration can be compared across buffers with different composition.
#'
#' **Class level** (one focal class `k`, set with `class`):
#'
#' | Code | Name | Definition | Units or range |
#' |:-----|:-----|:-----------|:---------------|
#' | `pland` | Percentage of landscape | `n_k / n * 100`. 0 when the class is absent | Percent |
#' | `ca` | Class area | `n_k * res^2 / 10000`. 0 when the class is absent | Hectares |
#' | `clumpy` | Clumpiness index | With `G = g_kk / (4 * n_k - Pmin(n_k))` and `P = n_k / n`: `(G - P) / (1 - P)` when `G >= P` or `P >= 0.5`, otherwise `(G - P) / (-P)`, as in `landscapemetrics`. `NA` when the class is absent or fills the buffer | Unitless. 0 when like adjacencies match the share `P` |
#'
#' # Surface metric definitions
#'
#' Each surface metric summarizes the values `z` of the `n` valid cells whose
#' centers lie within the radius. `s` is the sample standard deviation
#' (`"sq"`) and `res` the cell size in map units.
#'
#' | Code | Name | Definition | Units |
#' |:-----|:-----|:-----------|:------|
#' | `sa` | Average roughness | `sum(abs(z - mean(z))) / n` | Raster units |
#' | `sq` | Root mean square roughness | `sqrt(sum((z - mean(z))^2) / (n - 1))`. `NA` when `n < 2` | Raster units |
#' | `ssk` | Skewness | `sqrt(n * (n - 1)) / (n - 2) * (sum((z - mean(z))^3) / n) / s^3`. `NA` when `n < 3` or `s = 0` | Unitless. Positive when a few high values stand out |
#' | `sku` | Excess kurtosis | `(sum((z - mean(z))^4) / n) / s^4 - 3`. `NA` when `s = 0` | Unitless. 0 for normally distributed values |
#' | `sdq` | Root mean square slope | Square root of the mean of `((z - z_right)^2 + (z - z_below)^2) / res^2` over cells whose right and lower neighbors are both in the buffer | Raster units per map unit |
#' | `sdr` | Surface area ratio | `(sum(a) / (k * res^2) - 1) * 100`, where `a` is the three-dimensional area of each of the `k` blocks of 2 x 2 cell centers in the buffer | Percent |
#'
#' `"sa"`, `"sq"`, `"ssk"`, and `"sku"` match the `geodiv` package. `"sdq"` is
#' a true slope (rise over run). `geodiv::sdq()` assumes a cell spacing of 1,
#' so it equals `"sdq"` multiplied by `res`. `"sdr"` uses the real values and
#' cell size, so a plane tilted at slope `t` returns
#' `(sqrt(1 + t^2) - 1) * 100`. `geodiv::sdr()` rescales the surface to a
#' unit range first, so its values differ. Read `"sdq"` and `"sdr"` as a slope
#' and an area ratio only when the raster values and the map share units, as
#' elevation in meters does.
#'
#' # Kernel-weighted surface metrics
#'
#' With `weighted = TRUE`, each cell gets the kernel weight `w = w(d)` of its
#' distance `d` from the point (see [multiScaleR-kernels]), using the cells
#' extracted within `max_D`. With the weighted mean
#' `mu = sum(w * z) / sum(w)`:
#'
#' | Code | Weighted definition |
#' |:-----|:--------------------|
#' | `sa` | `sum(w * abs(z - mu)) / sum(w)` |
#' | `sq` | `sqrt(sum(w * (z - mu)^2) / sum(w))`, and 0 for constant values |
#' | `ssk` | `(sum(w * (z - mu)^3) / sum(w)) / sq^3`. `NA` when `sq = 0` |
#' | `sku` | `(sum(w * (z - mu)^4) / sum(w)) / sq^4 - 3`. `NA` when `sq = 0` |
#'
#' Weighted `"sq"` divides by the sum of weights, so it is a population-style
#' spread and differs from the `n - 1` form of the unweighted metric. The
#' weighted skewness and kurtosis use plain standardized moments, without the
#' small-sample adjustment of the unweighted `"ssk"`. Weighted `"sq"` is
#' smooth in sigma except where the weighted spread is zero. Absolute
#' deviations (`"sa"`) and the standardized higher moments have additional
#' points where they are not smooth or not defined. Within one evaluation,
#' weighted metrics reuse weights and centered moments only when their
#' source, scale, shape, and valid cells match exactly. Their scales are
#' still estimated separately. Under `kernel = "fixed"`, the fitted scale is
#' a hard cutoff radius. It has no standard error or Wald limits, although
#' the weighted definitions above still apply.
#'
#' # Check a metric before you optimize it
#'
#' Before estimating the radius of `"iji"`, `"siei"`, or `"msiei"`, inspect a
#' pilot [kernel_prep()] result or a few fixed radii. These metrics can be
#' undefined or nearly constant when buffers hold too few classes, too few
#' adjacencies between different classes, or nearly identical evenness at
#' every point. `"iji"` in particular returns `NA` with fewer than three
#' classes, and small buffers can make it `NA` everywhere. If so, use a larger
#' radius, increase `max_D`, or choose a more stable diversity, entropy, or
#' edge metric. A candidate scale that leaves a covariate undefined at a
#' modeled observation is rejected during optimization, so every candidate is
#' compared on the same observations.
#'
#' The adjacency and information-theory metrics accept a limited number of
#' classes, which guards against passing a continuous raster by mistake and
#' keeps raster projection times bounded. The limits are 200 classes for the
#' like-adjacency metrics (`"ai"`, `"pladj"`, `"clumpy"`) and 50 for the
#' metrics built from the full class-pair table (`"contag"`, `"iji"`, and the
#' information-theory metrics), whose cost grows with the square of the class
#' count.
#'
#' During optimization, the adjacency metrics reuse sorted cell and
#' neighbor-pair distances instead of rescanning every buffered cell at each
#' radius, and metrics on the same source and radius share their counts. This
#' saves time at the cost of retained memory, which
#' [estimate_multiscale_ram()] includes. It does not change radius boundaries
#' or metric definitions, and it does not make a hard radius smooth or remove
#' concerns about local optima.
#'
#' # Map a metric
#'
#' [kernel_scale.raster()] computes the same covariates for every raster cell.
#' Projection of `"sq"`, `"ssk"`, `"sku"`, `"sdq"`, `"sdr"`, and the landscape
#' metrics uses fast Fourier transform (FFT) convolution by default.
#' Unweighted `"sa"` uses an exact compiled moving-window pass instead,
#' because average roughness cannot be written as a convolution. `"sdq"` and
#' `"sdr"` use each cell's neighbors, so their projected values agree with
#' point extraction to within a small tolerance at the buffer boundary rather
#' than exactly. Set `fft = FALSE` in [kernel_scale.raster()] to use direct
#' moving-window calculations instead.
#'
#' @references
#' Nowosad, J., & Stepinski, T. F. (2019). Information theory as a consistent
#' framework for quantification and classification of landscape patterns.
#' *Landscape Ecology*, 34(9), 2091-2101.
#' \doi{10.1007/s10980-019-00830-x}
#'
#' @seealso [kernel_prep()], [multiScale_optim()], [kernel_scale.raster()],
#'   and [multiScaleR-kernels]. For worked examples, see
#'   `vignette("v04-landscape-metrics", package = "multiScaleR")` for the
#'   landscape metrics and
#'   `vignette("v05-surface-metrics", package = "multiScaleR")` for the
#'   surface metrics.
#'
#' @examples
#' ## Kernel-weighted mean only (the same as omitting scale_vars)
#' vars <- msr_vars(forest_prop = kernel_var("forest"))
#' vars
#'
#' ## Combine kernel, landscape, and surface covariates. The printed table
#' ## shows which scales are estimated (optimize = TRUE) and which are fixed.
#' vars <- msr_vars(
#'   forest_prop    = kernel_var("forest"),
#'   forest_ed      = landscape_var("forest", metric = "ed"),
#'   cover_shdi_500 = landscape_var("landcover", metric = "shdi", radius = 500),
#'   elev_sd        = surface_var("elevation", metric = "sq")
#' )
#' print(vars)
#'
#' ## Shannon diversity in bits within a fixed 250 m radius
#' landscape_var("landcover", metric = "shdi", radius = 250, base = 2)
#'
#' ## Interspersion of the whole buffer needs no focal class
#' landscape_var("landcover", metric = "iji", radius = 500)
#'
#' ## Class-level metrics need `class`, here clumpiness and cover of class 3
#' landscape_var("landcover", metric = "clumpy", radius = 500, class = 3)
#' landscape_var("landcover", metric = "pland", radius = 500, class = 3)
#'
#' ## Surface metrics with an estimated hard radius, a fixed radius, and
#' ## kernel weights (sigma estimated)
#' surface_var("elevation", metric = "sq")
#' surface_var("ndvi", metric = "sa", radius = 300)
#' surface_var("elevation", metric = "sq", weighted = TRUE)
#'
#' ## Compute covariates at survey points from the bundled landscape
#' land <- terra::rast(system.file("extdata", "landscape.tif",
#'                                 package = "multiScaleR"))
#' # land3 is continuous. Build a three-class layer from it
#' land$cover <- terra::classify(land$land3,
#'                               rbind(c(-Inf, 0.45, 1),
#'                                     c(0.45, 0.55, 2),
#'                                     c(0.55, Inf, 3)))
#' data("surv_pts")
#' vars <- msr_vars(
#'   land1_prop   = kernel_var("land1"),
#'   cover_shdi   = landscape_var("cover", metric = "shdi", radius = 400),
#'   land3_rough  = surface_var("land3", metric = "sq", radius = 300)
#' )
#' inputs <- kernel_prep(pts = terra::vect(surv_pts)[1:20, ],
#'                       raster_stack = land, max_D = 500,
#'                       scale_vars = vars, verbose = FALSE)
#' # Standardized covariates, one column per declared covariate
#' head(inputs$kernel_dat)
#'
#' @rdname msr_vars
#' @export
msr_vars <- function(...) {
  specs <- list(...)
  if (length(specs) == 0) {
    stop("At least one variable specification must be provided.", call. = FALSE)
  }

  spec_names <- names(specs)
  if (is.null(spec_names) || any(!nzchar(spec_names))) {
    stop("All `msr_vars()` specifications must be named.", call. = FALSE)
  }
  if (anyDuplicated(spec_names)) {
    stop("All `msr_vars()` names must be unique.", call. = FALSE)
  }
  validate_reserved_variable_names(spec_names, "covariate name in `msr_vars()`")

  for (i in seq_along(specs)) {
    if (!inherits(specs[[i]], "multiScaleR_var")) {
      stop(
        "All `msr_vars()` inputs must be created with `kernel_var()`, `landscape_var()`, or `surface_var()`.",
        call. = FALSE
      )
    }
  }

  out <- data.frame(
    covariate = spec_names,
    source = vapply(specs, `[[`, character(1), "source"),
    type = vapply(specs, `[[`, character(1), "type"),
    metric = vapply(specs, function(x) .msr_chr_or_na(x$metric), character(1)),
    radius = vapply(specs, function(x) .msr_num_or_na(x$radius), numeric(1)),
    optimize = vapply(specs, `[[`, logical(1), "optimize"),
    weighted = vapply(specs, function(x) isTRUE(x$weighted), logical(1)),
    base = vapply(specs, `[[`, numeric(1), "base"),
    classes_max = vapply(specs, function(x) .msr_num_or_na(x$classes_max), numeric(1)),
    class = vapply(specs, function(x) .msr_num_or_na(x$class), numeric(1)),
    stringsAsFactors = FALSE
  )

  class(out) <- c("multiScaleR_vars", "data.frame")
  out
}


.msr_hard_radius_covariates <- function(scale_vars, kernel = NULL,
                                         covariates = scale_vars$covariate) {
  # A uniform kernel changes support in cell steps for every fitted scale.
  # Explicit names also cover saved fits without scale specifications.
  if (identical(kernel, "fixed")) return(as.character(covariates))
  if (is.null(scale_vars) || is.null(scale_vars$type) ||
      is.null(scale_vars$covariate)) return(character())
  weighted <- if ("weighted" %in% names(scale_vars)) {
    as.logical(scale_vars$weighted)
  } else {
    rep(FALSE, nrow(scale_vars))
  }
  weighted[is.na(weighted)] <- FALSE
  scale_vars$covariate[
    scale_vars$type == "landscape" |
      (scale_vars$type == "surface" & !weighted)
  ]
}


#' @rdname msr_vars
#' @export
kernel_var <- function(source) {
  validate_character_scalar(source, "source")

  structure(
    list(
      type = "kernel",
      source = source,
      metric = NA_character_,
      radius = NA_real_,
      optimize = TRUE,
      weighted = FALSE,
      base = exp(1),
      classes_max = NA_real_,
      class = NA_real_
    ),
    class = "multiScaleR_var"
  )
}


#' @rdname msr_vars
#' @export
landscape_var <- function(source,
                          metric,
                          radius = NULL,
                          base = exp(1),
                          classes_max = NULL,
                          class = NULL) {
  validate_character_scalar(source, "source")
  metric <- match.arg(tolower(metric), .msr_landscape_metrics())
  validate_scalar_numeric(base, "base", positive = TRUE)
  if (base == 1) {
    stop("`base` must not equal 1.", call. = FALSE)
  }

  if (!is.null(radius)) {
    validate_scalar_numeric(radius, "radius", positive = TRUE)
  }
  if (!is.null(classes_max)) {
    validate_scalar_numeric(classes_max, "classes_max", positive = TRUE)
  }

  # Modified Simpson diversity is -ln(sum(p^2)) by definition (FRAGSTATS and
  # landscapemetrics offer no base choice), so the metric code never applies
  # `base` to it. Warn instead of silently ignoring a requested base, and store
  # the natural log so the specification records what is actually computed.
  if (metric %in% c("msidi", "msiei") && !isTRUE(all.equal(base, exp(1)))) {
    warning(
      sprintf("`base` does not apply to `%s`, which always uses the natural log; ignoring `base = %s`.",
              metric, format(base)),
      call. = FALSE
    )
    base <- exp(1)
  }

  # Only the entropy metrics use `base`. Shannon evenness, contagion, IJI and
  # relative mutual information are ratios of logarithms, so every base gives
  # the same value. The remaining metrics take no logarithm, so a requested
  # base would have no effect. Say so, as for the modified Simpson metrics.
  base_free <- setdiff(.msr_landscape_metrics(),
                       c("shdi", "shei", "ent", "condent", "joinent", "mutinf",
                         "relmutinf", "iji", "contag", "msidi", "msiei"))
  if (metric %in% base_free && !isTRUE(all.equal(base, exp(1)))) {
    warning(
      sprintf("`base` does not apply to `%s`, so `base = %s` is ignored.",
              metric, format(base)),
      call. = FALSE
    )
    base <- exp(1)
  }

  # Relative patch richness has no defined value without the maximum class
  # count: the metric code returns NA, which previously surfaced later as an
  # unrelated "all NA" error in kernel_prep(). Fail here, where the covariate
  # is defined and the missing argument can be named.
  if (metric == "rpr" && is.null(classes_max)) {
    stop(
      "`classes_max` is required for relative patch richness (`rpr`): supply the maximum number of classes possible in the landscape.",
      call. = FALSE
    )
  }

  # The class-level metrics (clumpy, pland, ca) describe a single focal class and
  # therefore require `class`; the landscape-level metrics do not accept it.
  is_class_metric <- metric %in% .msr_class_metrics()
  if (is_class_metric) {
    if (is.null(class)) {
      stop(
        sprintf("`class` is required for the class-level metric `%s`.", metric),
        call. = FALSE
      )
    }
    validate_scalar_numeric(class, "class", integerish = TRUE)
  } else if (!is.null(class)) {
    stop(
      sprintf("`class` applies only to the class-level metrics (%s); `%s` is a landscape-level metric.",
              paste(.msr_class_metrics(), collapse = ", "), metric),
      call. = FALSE
    )
  }

  structure(
    list(
      type = "landscape",
      source = source,
      metric = metric,
      radius = if (is.null(radius)) NA_real_ else radius,
      optimize = is.null(radius),
      weighted = FALSE,
      base = base,
      classes_max = if (is.null(classes_max)) NA_real_ else classes_max,
      class = if (is.null(class)) NA_real_ else as.numeric(class)
    ),
    class = "multiScaleR_var"
  )
}


#' @rdname msr_vars
#' @export
surface_var <- function(source,
                        metric,
                        radius = NULL,
                        weighted = FALSE) {
  validate_character_scalar(source, "source")
  metric <- match.arg(tolower(metric), .msr_surface_metrics())
  validate_scalar_logical(weighted, "weighted")

  if (isTRUE(weighted)) {
    if (metric %in% .msr_surface_neighbor_metrics()) {
      stop(
        sprintf(
          "`weighted = TRUE` is not supported for the slope metric `%s`; kernel weighting applies to the distribution metrics (%s).",
          metric,
          paste(setdiff(.msr_surface_metrics(), .msr_surface_neighbor_metrics()),
                collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (!is.null(radius)) {
      stop(
        "`weighted = TRUE` surface metrics optimize the kernel scale, so `radius` must be NULL.",
        call. = FALSE
      )
    }
  }

  if (!is.null(radius)) {
    validate_scalar_numeric(radius, "radius", positive = TRUE)
  }

  structure(
    list(
      type = "surface",
      source = source,
      metric = metric,
      radius = if (is.null(radius)) NA_real_ else radius,
      optimize = is.null(radius),
      weighted = weighted,
      base = exp(1),
      classes_max = NA_real_,
      class = NA_real_
    ),
    class = "multiScaleR_var"
  )
}


#' @rdname msr_vars
#' @export
#' @method print multiScaleR_vars
print.multiScaleR_vars <- function(x, ...) {
  cat("multiScaleR variable specifications:\n")
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}


.msr_chr_or_na <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    NA_character_
  } else {
    as.character(x[[1]])
  }
}


.msr_num_or_na <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    NA_real_
  } else {
    as.numeric(x[[1]])
  }
}


.msr_landscape_metrics <- function() {
  c("shdi", "shei", "sidi", "siei", "msidi", "msiei", "pr", "prd", "rpr",
    "ta", "ed", "te", "lsi", "ai", "pladj", "contag", "iji", "ent", "condent",
    "joinent", "mutinf", "relmutinf", "clumpy", "pland", "ca")
}


.msr_composition_metrics <- function() {
  c("shdi", "shei", "sidi", "siei", "msidi", "msiei", "pr", "prd", "rpr",
    "ta")
}


# Class-level metrics target a single focal class, supplied via the `class`
# argument of `landscape_var()`: clumpiness (CLUMPY, an adjacency metric) and the
# basic composition metrics percentage of landscape (PLAND) and class area (CA).
.msr_class_metrics <- function() {
  c("clumpy", "pland", "ca")
}


.msr_surface_metrics <- function() {
  c("sa", "sq", "ssk", "sku", "sdq", "sdr")
}


# Surface metrics that require per-cell neighbor information (cell IDs), as
# opposed to the pointwise distribution metrics: `sdq` (RMS slope, forward
# differences to each cell's neighbors) and `sdr` (surface area ratio, which
# triangulates each cell quad).
.msr_surface_neighbor_metrics <- function() {
  c("sdq", "sdr")
}


.msr_edge_metrics <- function() {
  c("ed", "te", "lsi")
}


# Landscape-level metrics derived from the per-window class-pair co-occurrence
# matrix: aggregation/contagion (ai, pladj, contag), interspersion (iji), and the
# Nowosad & Stepinski (2019) information-theory family (ent, condent, joinent,
# mutinf, relmutinf). The class-level CLUMPY is built from the same matrix but is
# routed via `.msr_class_metrics()` because it additionally needs a focal class.
.msr_adjacency_metrics <- function() {
  c("ai", "pladj", "contag", "iji", "ent", "condent", "joinent", "mutinf",
    "relmutinf")
}


.msr_default_scale_vars <- function(raster_stack) {
  vars <- data.frame(
    covariate = names(raster_stack),
    source = names(raster_stack),
    type = "kernel",
    metric = NA_character_,
    radius = NA_real_,
    optimize = TRUE,
    weighted = FALSE,
    base = exp(1),
    classes_max = NA_real_,
    class = NA_real_,
    stringsAsFactors = FALSE
  )
  class(vars) <- c("multiScaleR_vars", "data.frame")
  vars
}


.msr_validate_scale_vars <- function(scale_vars, raster_stack, kernel = "gaussian") {
  if (is.null(scale_vars)) {
    scale_vars <- .msr_default_scale_vars(raster_stack)
  }

  if (!inherits(scale_vars, "multiScaleR_vars")) {
    stop("`scale_vars` must be created with `msr_vars()`.", call. = FALSE)
  }

  scale_vars <- as.data.frame(scale_vars, stringsAsFactors = FALSE)
  missing_sources <- setdiff(scale_vars$source, names(raster_stack))
  if (length(missing_sources) > 0) {
    stop(
      paste0(
        "The following `scale_vars` sources or optimized covariate layers are not raster layers: ",
        paste(missing_sources, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (anyDuplicated(scale_vars$covariate)) {
    stop("`scale_vars` covariate names must be unique.", call. = FALSE)
  }

  # Hard-neighborhood metrics use their radius alone. An expow mixture may
  # retain legacy shape slots, but those slots never enter these transforms.

  class(scale_vars) <- c("multiScaleR_vars", "data.frame")
  scale_vars
}


.msr_optimized_scale_vars <- function(scale_vars) {
  if (is.null(scale_vars)) {
    return(NULL)
  }
  scale_vars[scale_vars$optimize, , drop = FALSE]
}


.msr_scale_vars_in_model <- function(scale_vars, multiScaleR) {
  if (is.null(scale_vars) ||
      is.null(multiScaleR) ||
      !inherits(multiScaleR, "multiScaleR")) {
    return(scale_vars)
  }

  used_covariates <- tryCatch(
    intersect(scale_vars$covariate,
      if (methods::is(multiScaleR$opt_mod, "unmarkedFit"))
        .unmarked_model_predictors(multiScaleR$opt_mod) else
        .model_predictors(.analysis_model(multiScaleR$opt_mod))),
    error = function(e) character(0)
  )

  if (length(used_covariates) == 0 && !is.null(multiScaleR$scale_est)) {
    used_covariates <- intersect(scale_vars$covariate,
                                 row.names(multiScaleR$scale_est))
  }

  if (length(used_covariates) == 0) {
    return(scale_vars)
  }

  scale_vars[scale_vars$covariate %in% used_covariates, , drop = FALSE]
}


.msr_optimized_covariates <- function(scale_vars, cov_df = NULL) {
  if (is.null(scale_vars)) {
    return(colnames(cov_df[[1]]))
  }
  .msr_optimized_scale_vars(scale_vars)$covariate
}


.msr_needs_cells <- function(scale_vars) {
  any(scale_vars$type == "landscape" &
        scale_vars$metric %in% c(.msr_edge_metrics(), .msr_adjacency_metrics(),
                                 "clumpy")) ||
    any(scale_vars$type == "surface" &
          scale_vars$metric %in% .msr_surface_neighbor_metrics())
}


.msr_extract_cells <- function(cov_df) {
  cells <- attr(cov_df, "cell", exact = TRUE)
  if (is.null(cells)) {
    stop(
      "Raster cell IDs are required for this landscape metric but were not cached.",
      call. = FALSE
    )
  }
  cells
}


.msr_eval_scale_vars <- function(d,
                                 cov_df,
                                 scale_vars,
                                 sigma,
                                 shape,
                                 kernel,
                                 unit_conv,
                                 resolution,
                                 n_cols,
                                 covariates = scale_vars$covariate,
                                 validate = TRUE,
                                 prepared_neighbors = NULL,
                                 prepared_adjacency = NULL,
                                 prepared_edges = NULL) {
  out <- numeric(length(covariates))
  names(out) <- covariates

  opt_vars <- .msr_optimized_scale_vars(scale_vars)
  param_covariates <- covariates[covariates %in% opt_vars$covariate]

  weighted_values <- .msr_shared_weighted_metrics(d, cov_df, scale_vars,
    covariates, sigma, shape, kernel, param_covariates, validate)

  adjacency_values <- if (isTRUE(validate)) numeric(0) else
    .msr_eval_prepared_adjacency(prepared_adjacency, scale_vars, covariates,
                               sigma, param_covariates, unit_conv)
  edge_values <- if (isTRUE(validate)) numeric(0) else
    .msr_eval_edge_neighbors(prepared_edges, scale_vars, covariates,
                             sigma, param_covariates, unit_conv)

  for (j in seq_along(covariates)) {
    covariate <- covariates[[j]]
    spec <- scale_vars[scale_vars$covariate == covariate, , drop = FALSE]
    if (nrow(spec) != 1) {
      stop("Could not find a unique `scale_vars` specification for `", covariate, "`.",
           call. = FALSE)
    }

    if (covariate %in% names(weighted_values)) {
      out[[j]] <- weighted_values[[covariate]]
      next
    }

    if (covariate %in% names(edge_values)) {
      out[[j]] <- edge_values[[covariate]]
      next
    }
    if (covariate %in% names(adjacency_values)) {
      out[[j]] <- adjacency_values[[covariate]]
      next
    }
    prepared <- prepared_neighbors[[covariate]]
    if (!isTRUE(validate) && !is.null(prepared)) {
      radius <- if (isTRUE(spec$optimize))
        sigma[[match(covariate, param_covariates)]] else spec$radius / unit_conv
      # Preserve the legacy empty-neighborhood error path.
      if (is.finite(radius) && radius >= prepared$min_distance) {
        out[[j]] <- .msr_eval_surface_neighbors(prepared, radius)
        next
      }
    }
    values <- cov_df[, spec$source, drop = FALSE]

    if (identical(spec$type, "kernel")) {
      param_idx <- match(covariate, param_covariates)
      out[[j]] <- scale_type(
        d = d,
        kernel = kernel,
        sigma = sigma[[param_idx]],
        shape = if (is.null(shape)) NULL else shape[[param_idx]],
        r_stack.df = values
      )[[1]]
      next
    }

    radius <- if (isTRUE(spec$optimize)) {
      sigma[[match(covariate, param_covariates)]]
    } else {
      spec$radius / unit_conv
    }

    metric <- spec$metric
    if (identical(spec$type, "surface")) {
      if (isTRUE(spec$weighted)) {
        param_idx <- match(covariate, param_covariates)
        out[[j]] <- .surface_weighted_by_buffer(
          d = d,
          r_stack.df = values,
          kernel = kernel,
          sigma = sigma[[param_idx]],
          shape = if (is.null(shape)) NULL else shape[[param_idx]],
          metric = metric,
          validate = validate
        )[[1]]
      } else if (metric %in% .msr_surface_neighbor_metrics()) {
        out[[j]] <- .surface_neighbor_by_buffer(
          d = d,
          r_stack.df = values,
          cells = .msr_extract_cells(cov_df),
          radius = radius,
          resolution = resolution,
          n_cols = n_cols,
          metric = metric,
          validate = validate
        )[[1]]
      } else {
        out[[j]] <- .surface_metric_by_buffer(
          d = d,
          r_stack.df = values,
          radius = radius,
          metric = metric,
          validate = validate
        )[[1]]
      }
    } else if (metric %in% c(.msr_composition_metrics(), "pland", "ca")) {
      out[[j]] <- .landscape_composition_by_buffer(
        d = d,
        r_stack.df = values,
        radius = radius,
        metric = metric,
        base = spec$base,
        resolution = resolution,
        classes_max = if (is.na(spec$classes_max)) NULL else spec$classes_max,
        focal_class = if (is.na(spec$class)) NULL else spec$class,
        validate = validate
      )[[1]]
    } else if (metric %in% .msr_edge_metrics()) {
      out[[j]] <- .landscape_edge_by_buffer(
        d = d,
        r_stack.df = values,
        cells = .msr_extract_cells(cov_df),
        radius = radius,
        resolution = resolution,
        n_cols = n_cols,
        metric = metric,
        validate = validate
      )[[1]]
    } else if (metric %in% c(.msr_adjacency_metrics(), "clumpy")) {
      out[[j]] <- .landscape_adjacency_by_buffer(
        d = d,
        r_stack.df = values,
        cells = .msr_extract_cells(cov_df),
        radius = radius,
        n_cols = n_cols,
        metric = metric,
        base = spec$base,
        focal_class = if (is.na(spec$class)) NULL else spec$class,
        validate = validate
      )[[1]]
    } else {
      stop("Unsupported landscape metric '", metric, "'.", call. = FALSE)
    }
  }

  out
}


.msr_merge_scl_params <- function(primary, fallback) {
  if (is.null(fallback)) {
    return(primary)
  }
  if (is.null(primary)) {
    return(fallback)
  }

  out <- fallback
  out$mean[names(primary$mean)] <- primary$mean
  out$sd[names(primary$sd)] <- primary$sd
  out
}

# Full-length computational shape vector for the legacy expow parameter layout.
# A hard-radius transform ignores this placeholder; its reported shape is NA.
.msr_shape_values <- function(object) {
  if (is.null(object$shape_est)) return(NULL)
  value <- stats::setNames(object$shape_est$Mean, rownames(object$shape_est))
  hard <- intersect(names(value), .msr_hard_radius_covariates(object$kernel_inputs$scale_vars,
    object$kernel_inputs$kernel, names(value)))
  value[hard] <- 1
  value
}


.msr_scale_vars_raster <- function(raster_stack,
                                   scale_vars,
                                   sigma,
                                   shape,
                                   kernel,
                                   pct_wt,
                                   fft,
                                   na.rm,
                                   verbose,
                                   buffer_radius = NULL) {
  opt_vars <- .msr_optimized_scale_vars(scale_vars)
  out <- vector("list", nrow(scale_vars))
  names(out) <- scale_vars$covariate
  convolve <- if (isTRUE(fft)) fft_convolution else .msr_focal_convolution

  for (i in seq_len(nrow(scale_vars))) {
    spec <- scale_vars[i, , drop = FALSE]
    source_raster <- raster_stack[[spec$source]]

    if (identical(spec$type, "kernel")) {
      param_idx <- match(spec$covariate, opt_vars$covariate)
      if (isTRUE(verbose)) {
        cat(paste0(
          "\nSmoothing spatRaster variable '", spec$covariate,
          "' from layer '", spec$source,
          "' at sigma = ", floor(sigma[[param_idx]]), "\n"
        ))
      }
      out[[i]] <- .msr_kernel_raster_one(
        raster = source_raster,
        sigma = sigma[[param_idx]],
        shape = if (is.null(shape)) NULL else shape[[param_idx]],
        kernel = kernel,
        pct_wt = pct_wt,
        fft = fft,
        na.rm = na.rm,
        buffer_radius = buffer_radius
      )
      names(out[[i]]) <- spec$covariate
      next
    }

    radius <- if (isTRUE(spec$optimize)) {
      sigma[[match(spec$covariate, opt_vars$covariate)]]
    } else {
      spec$radius
    }
    if (!is.null(buffer_radius) && !isTRUE(spec$weighted) && radius > buffer_radius)
      stop("A metric radius exceeds the prepared buffer; its truncated footprint cannot be projected by this matching mode.", call. = FALSE)

    if (isTRUE(verbose)) {
      if (identical(spec$type, "surface") && isTRUE(spec$weighted)) {
        cat(paste0(
          "\nCalculating weighted surface texture metric '", spec$metric,
          "' for variable '", spec$covariate,
          "' at sigma = ", floor(radius), "\n"
        ))
      } else {
        metric_label <- if (identical(spec$type, "surface")) {
          "surface texture metric"
        } else {
          "landscape metric"
        }
        cat(paste0(
          "\nCalculating ", metric_label, " '", spec$metric,
          "' for variable '", spec$covariate,
          "' at radius = ", floor(radius), "\n"
        ))
      }
    }

    if (identical(spec$type, "surface")) {
      if (isTRUE(spec$weighted)) {
        param_idx <- match(spec$covariate, opt_vars$covariate)
        out[[i]] <- .surface_weighted_raster_fft(
          raster = source_raster,
          kernel = kernel,
          sigma = sigma[[param_idx]],
          shape = if (is.null(shape)) NULL else shape[[param_idx]],
          metric = spec$metric,
          pct_wt = pct_wt,
          na.rm = na.rm,
          buffer_radius = buffer_radius, convolve = convolve
        )
      } else {
        out[[i]] <- .surface_metric_raster_fft(
          raster = source_raster,
          radius = radius,
          metric = spec$metric,
          na.rm = na.rm, convolve = convolve
        )
      }
    } else if (spec$metric %in% c(.msr_composition_metrics(), "pland", "ca")) {
      out[[i]] <- .landscape_composition_raster_fft(
        raster = source_raster,
        radius = radius,
        metric = spec$metric,
        base = spec$base,
        classes_max = if (is.na(spec$classes_max)) NULL else spec$classes_max,
        focal_class = if (is.na(spec$class)) NULL else spec$class,
        na.rm = na.rm, convolve = convolve
      )
    } else if (spec$metric %in% .msr_edge_metrics()) {
      out[[i]] <- .landscape_edge_raster_fft(
        raster = source_raster,
        radius = radius,
        metric = spec$metric,
        na.rm = na.rm, convolve = convolve
      )
    } else {
      out[[i]] <- .landscape_adjacency_raster_fft(
        raster = source_raster,
        radius = radius,
        metric = spec$metric,
        base = spec$base,
        focal_class = if (is.na(spec$class)) NULL else spec$class,
        na.rm = na.rm, convolve = convolve
      )
    }
    if (!isTRUE(na.rm)) {
      support <- if (identical(spec$type, "surface") && isTRUE(spec$weighted))
        .surface_kernel_weight_matrix(source_raster, radius,
          if (is.null(shape)) NULL else shape[[param_idx]], kernel, pct_wt) != 0 else
        .landscape_circle_kernel(radius, terra::res(source_raster)[1L]) != 0
      finite <- !is.na(terra::as.matrix(source_raster, wide = TRUE))
      count <- convolve(matrix(as.numeric(finite), nrow(finite), ncol(finite)),
        support * 1, fun = "sum", na.rm = TRUE)
      value <- terra::values(out[[i]])
      value[as.vector(round(count) != sum(support)), 1L] <- NA_real_
      terra::values(out[[i]]) <- value
    }
    names(out[[i]]) <- spec$covariate
  }

  terra::rast(out)
}


.msr_kernel_raster_one <- function(raster,
                                   sigma,
                                   shape,
                                   kernel,
                                   pct_wt,
                                   fft,
                                   na.rm,
                                   buffer_radius = NULL) {
  if (!is.null(buffer_radius)) {
    if (!isTRUE(na.rm))
      stop("Prepared-buffer projection requires finite-cell normalization.", call. = FALSE)
    convolve <- if (isTRUE(fft)) fft_convolution else .msr_focal_convolution
    weights <- .surface_kernel_weight_matrix(raster, sigma, shape, kernel,
                                              pct_wt, buffer_radius)
    values <- terra::as.matrix(raster, wide = TRUE)
    valid <- is.finite(values)
    zeroed <- values
    zeroed[!valid] <- 0
    # Sum convolutions normalize against actual mapped finite support, including
    # raster edges. Keep the original central-cell mask as the projection domain.
    s0 <- convolve(matrix(as.numeric(valid), nrow(values), ncol(values)),
                           weights, fun = "sum", na.rm = TRUE)
    s1 <- convolve(zeroed, weights, fun = "sum", na.rm = TRUE)
    result <- s1 / s0
    result[!is.finite(s0) | s0 <= 0] <- NA_real_
    result[t(!valid)] <- NA_real_
    out <- terra::rast(raster)
    terra::values(out) <- as.vector(result)
    return(out)
  }
  mx <- kernel_dist(kernel = kernel,
                    dimension = "1d",
                    sigma = sigma,
                    shape = shape,
                    prob = pct_wt)

  r_res <- terra::res(raster)[1]
  focal_d <- ceiling(mx / r_res) * 2
  if ((focal_d %% 2) == 0) {
    focal_d <- focal_d + 1
  }

  mat <- matrix(0, focal_d, focal_d)
  r_wt <- terra::rast(mat)
  terra::crs(r_wt) <- terra::crs(raster)
  cntr_crd <- terra::xyFromCell(r_wt, ceiling(length(mat) / 2))
  cell_crds <- terra::crds(r_wt)
  r_wt[] <- fields::rdist(cntr_crd, cell_crds)[1, ] * r_res
  r_wt[] <- scale_type_r(d = as.vector(r_wt),
                         kernel = kernel,
                         sigma = sigma,
                         shape = shape,
                         output = "wts")

  wt_mat <- terra::as.matrix(r_wt, wide = TRUE)
  if (isTRUE(fft)) {
    mat <- terra::as.matrix(raster, wide = TRUE)
    out_mat <- fft_convolution(mat,
                               wt_mat,
                               fun = "mean",
                               na.rm = na.rm)
    out <- terra::rast(raster)
    terra::values(out) <- as.vector(out_mat)
    out
  } else {
    terra::focal(raster,
                 w = wt_mat,
                 fun = "mean",
                 na.rm = na.rm,
                 na.policy = "omit",
                 expand = FALSE)
  }
}
