#' Example data frame
#'
#' Example count data to be used for optimizing scales of effect
#'
#' @format ## `count_data`
#' @name count_data
#' @usage data(count_data)
#' A data frame with 75 rows and 2 columns. Data were simulated from a Poisson distribution with an intercept of 0.5, a slope of 0.75, and scale of effect (sigma) of 75
#' \describe{
#'   \item{y}{Simulated counts at spatial locations}
#'   \item{had}{A scaled, weighted mean estimate of the habitat surrounding a point}#'
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
