#' Simulated counts for the hab example
#'
#' Simulated counts at the 75 points in [pts], for learning to estimate a
#' scale of effect with the `hab.tif` raster.
#'
#' @name count_data
#' @usage data(count_data)
#' @format A data frame with 75 rows, one per point in [pts], and 2 columns:
#' * `y`: integer count at each point.
#' * `hab`: kernel-weighted mean of the binary `hab.tif` raster around each
#'   point, centered and scaled to unit variance.
#'
#' @details Counts were simulated from a Poisson distribution with an
#'   intercept of 0.5, a `hab` coefficient of 0.75 (log scale), and a Gaussian
#'   scale of effect (sigma) of 75 map units. A fit that recovers a sigma near
#'   75 and a coefficient near 0.75 is working as intended.
#' @source Simulated for this package.
#' @seealso [pts], [hab], [multiScale_optim()].
#' @examples
#' data(count_data)
#' head(count_data)
#' # Starting model: hab will be replaced by values at each candidate scale
#' glm(y ~ hab, family = poisson, data = count_data)
NULL

#' Simulated counts for the vignette example
#'
#' Simulated counts at the 100 points in [surv_pts], used as the running
#' example in the vignettes with the `landscape.tif` raster.
#'
#' @name landscape_counts
#' @usage data(landscape_counts)
#' @format A data frame with 100 rows, one per point in [surv_pts], and 2
#'   columns:
#' * `counts`: integer count at each point.
#' * `site`: a habitat variable measured at the point itself, which has no
#'   scale of effect.
#'
#' @details Counts were simulated from a Poisson distribution with an
#'   intercept of 0.25 and log-scale coefficients of -0.5 for `land1`, 0.3 for
#'   `site`, and 0.7 for `land2`. The true Gaussian scales of effect (sigma) are
#'   250 map units for `land1` and 500 for `land2`. The `land1` and `land2`
#'   covariates are not stored here; build them from [landscape] with
#'   [kernel_prep()].
#' @source Simulated for this package.
#' @seealso [surv_pts], [landscape],
#'   `vignette("v01-quickstart", package = "multiScaleR")`.
#' @examples
#' data(landscape_counts)
#' summary(landscape_counts)
NULL

#' Survey points for the hab example
#'
#' Locations of the 75 survey points that go with [count_data] and the
#' `hab.tif` raster.
#'
#' @name pts
#' @usage data(pts)
#' @format An `sf` point object with 75 features, no coordinate reference
#'   system, and coordinates in the map units of `hab.tif` (from 400 to 3600).
#'   Its one attribute column, `obs`, repeats the counts in `count_data$y`.
#' @source Simulated for this package.
#' @seealso [count_data], [hab], [kernel_prep()].
#' @examples
#' data(pts)
#' pts
#' hab <- terra::rast(system.file("extdata", "hab.tif",
#'                                package = "multiScaleR"))
#' terra::plot(hab)
#' terra::points(terra::vect(pts), pch = 19)
NULL

#' Survey points for the vignette example
#'
#' Locations of the 100 survey points that go with [landscape_counts] and the
#' `landscape.tif` raster.
#'
#' @name surv_pts
#' @usage data(surv_pts)
#' @format An `sf` point object with 100 features, no coordinate reference
#'   system, and coordinates in the map units of `landscape.tif`. Its one
#'   attribute column, `obs`, is a numeric column carried over from the
#'   simulation; the vignettes use `landscape_counts$counts` as the response.
#' @source Simulated for this package.
#' @seealso [landscape_counts], [landscape],
#'   `vignette("v01-quickstart", package = "multiScaleR")`.
#' @examples
#' data(surv_pts)
#' land <- terra::rast(system.file("extdata", "landscape.tif",
#'                                 package = "multiScaleR"))
#' terra::plot(land$land1)
#' terra::points(terra::vect(surv_pts), pch = 19)
NULL

#' Example habitat raster
#'
#' A binary habitat raster stored as `hab.tif` in the package's `extdata`
#' folder, used with [pts] and [count_data].
#'
#' @name hab
#' @format A `SpatRaster` with one layer, `hab`, of 400 by 400 cells of 10
#'   map units (an extent of 0 to 4000 in each direction) and no coordinate
#'   reference system. Cells are 1 for habitat and 0 otherwise.
#' @source Simulated for this package.
#' @seealso [pts], [count_data].
#' @examples
#' hab <- terra::rast(system.file("extdata", "hab.tif",
#'                                package = "multiScaleR"))
#' hab
#' terra::plot(hab)
NULL

#' Example landscape raster
#'
#' Three simulated landscape layers stored as `landscape.tif` in the
#' package's `extdata` folder, used with [surv_pts] and [landscape_counts] in
#' the vignettes.
#'
#' @name landscape
#' @format A `SpatRaster` with three layers of 200 by 200 cells of 20 map
#'   units (an extent of 0 to 4000 in each direction) and no coordinate
#'   reference system. Values range from 0 to 1:
#' * `land1`: binary surface with low spatial autocorrelation.
#' * `land2`: continuous surface with low spatial autocorrelation.
#' * `land3`: continuous surface with high spatial autocorrelation.
#'
#' @source Simulated for this package.
#' @seealso [surv_pts], [landscape_counts].
#' @examples
#' land <- terra::rast(system.file("extdata", "landscape.tif",
#'                                 package = "multiScaleR"))
#' land
#' terra::plot(land)
NULL
