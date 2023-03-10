% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_prep.R
\name{kernel_prep}
\alias{kernel_prep}
\title{Kernel Scale Preparation}
\usage{
kernel_prep(pts,
            raster_stack,
            max_D,
            kernel = 'gaussian',
            sigma = NULL,
            shape = NULL,
            projected = TRUE)
}
\arguments{
\item{pts}{Point locations provided as `SpatVector` or `sf` objects}

\item{raster_stack}{Raster layer(s) of class `SpatRaster`}

\item{max_D}{The maximum distance to consider during the scale optimization}

\item{kernel}{Kernel function to be used ('gaussian', 'exp', 'fixed', 'expow'; Default: 'gaussian')}

\item{sigma}{Initial values for optimizing the scale parameter. Default: NULL, initial values will be automatically generated. This is recommended.}

\item{shape}{Initial values for optimizing the shape parameter if using exponential power kernel. Default: NULL, starting values will be automatically generated. This is recommended.}

\item{projected}{Logical. Are `pts` and `raster_stack` projected. Function currently requires that both are projected, Default: TRUE}
}
\value{
A list of class `multiscaleR` with necessary elements to conduct scale optimization using the `multiScale_optim` function
}
\description{
Function to prepare data inputs for kernel scale analysis
}
\details{
Spatial point locations and raster layers should have a defined projection and be the same CRS. If providing starting values for `sigma` or `shape`, it must be a vector of length equal to the number of raster layers for which scale is being assessed and should be provided in the unit of the used projection. When specifying `max_D`, ensure that your raster layers adequately extend beyond the points provided so that the surrounding landscape can be meaningfully sampled during scale optimization.
}
\examples{
library(terra)
pts <- vect(cbind(c(3,5,7),
                  c(7,5,3)))

mat_list <- list(r1 = rast(matrix(rnorm(100),
                                  nrow = 10)),
                 r2 = rast(matrix(rnorm(100),
                                  nrow = 10)))
rast_stack <- rast(mat_list)
kernel_inputs <- kernel_prep(pts = pts,
                             raster_stack = rast_stack,
                             max_D = 2,
                             kernel = 'gaussian',
                             sigma = NULL)
}
