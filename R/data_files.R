#' Example data frame
#'
#' Example count data to be used for optimizing scales of effect
#'
#' @format ## `count_data`
#' @name count_data
#' @usage data(count_data)
#' A data frame with 100 rows and 2 columns. Data were simulated from a Poisson distribution with an intercept of 0.25; land1 effect = -0.5; site effect = 0.3; land2 effect = 0.7. True simulated Gaussian scale effects (sigma): land1 = 250; land2 = 500.
#' \describe{
#'   \item{counts}{Simulated counts at spatial locations}
#'   \item{site}{Scaled values from the local point level}
#'
#'   ...
#' }
NULL

#' Example data frame
#'
#' Example count data to be used vignette document example
#'
#' @format ## `landscape_counts`
#' @name landscape_counts
#' @usage data(landscape_counts)
#' A data frame with 100 rows and 2 columns. Data were simulated from a Poisson distribution with an intercept of 0.5, a slope of 0.75, and scale of effect (sigma) of 75
#' \describe{
#'   \item{y}{Simulated counts at spatial locations}
#'   \item{hab}{A scaled, weighted mean estimate of the habitat surrounding a point}
#'
#'   ...
#' }
NULL

#' Spatial sample points
#'
#' Example point file for optimizing scales of effect
#'
#' @format ## `pts`
#' @name pts
#' @usage data(pts)
#' An sf class point object:
#' \describe{
#'   \item{pts}{spatial location of points}
#'   ...
#' }
NULL

#' Spatial sample points
#'
#' Example point file for use with vignette document example
#'
#' @format  `surv_pts`
#' @name surv_pts
#' @usage data(surv_pts)
#' An sf class point object:
#' \describe{
#'   \item{pts}{100 spatial point locations}
#'   ...
#' }
NULL

#' Example raster
#'
#' Example habitat raster for optimizing scales of effect
#'
#' @format ## `hab`
#' @name hab
#' @usage
#' hab <- terra::rast(system.file("data/hab.tif", package = 'multiScaleR'))
#' A binary SpatRaster object
#' \describe{
#'   \item{hab}{A binary raster}
#'   ...
#' }
NULL

#' Simulated raster
#'
#' Raster data for use with vignette example
#'
#' @format `landscape_rast`
#' @name landscape
#' @usage
#' land_rast <- terra::rast(system.file("data/landscape.tif", package = 'multiScaleR'))
#'
#' \describe{
#'   A spatRaster object with three surfaces:
#'   \item{land1}{A binary landscape surface with low correlation}
#'   \item{land2}{A continuous landscape surface with low autocorrelation}
#'   \item{land3}{A continuous landscape surface with high autocorrelation}
#'   ...
#' }
NULL
