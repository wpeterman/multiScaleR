#' Retrieve diagnostics from multiScaleR objects
#'
#' Returns structured warning/diagnostic information stored on fitted
#' \code{multiScaleR} objects. Diagnostics are populated automatically during
#' \code{\link{multiScale_optim}} and flag potential issues with the
#' optimization result.
#'
#' @param object An object to inspect. Must be of class \code{multiScaleR}.
#' @param ... Additional arguments passed to methods.
#'
#' @return A named list with up to three elements:
#' \describe{
#'   \item{\code{max_distance}}{A list describing whether the estimated scale of
#'     effect approaches or exceeds \code{max_D}. Fields include:
#'     \code{triggered} (logical), \code{variables} (names of affected
#'     covariates), \code{effective_distance} (estimated 90\% kernel distance
#'     per covariate), \code{max_D} (the limit used), \code{ratio}
#'     (\code{max_D / effective_distance}), and \code{suggested_max_D} (a
#'     recommended minimum \code{max_D} value). Triggered when
#'     \code{max_D / effective_distance < 2}.}
#'   \item{\code{sigma_precision}}{A list describing whether sigma estimates are
#'     imprecise. Fields include: \code{triggered} (logical), \code{variables}
#'     (names of affected covariates), \code{estimate} (sigma means),
#'     \code{se} (standard errors), and \code{se_to_mean} (ratio of SE to
#'     mean). Triggered when \code{SE / Mean >= 0.5} for any covariate. \code{NULL}
#'     when no precision concern was flagged.}
#'   \item{\code{shape_precision}}{Identical structure to \code{sigma_precision}
#'     but for the shape parameter of the exponential power kernel.  \code{NULL}
#'     unless \code{kernel = "expow"} was used and a precision concern arose.}
#' }
#'
#' @details
#' All three diagnostics are evaluated automatically at the end of
#' \code{\link{multiScale_optim}}. Console warnings are printed when any
#' diagnostic is triggered. \code{diagnostics()} provides programmatic access
#' to the same information without re-running the model.
#'
#' When \code{max_distance} is triggered, consider re-running
#' \code{\link{kernel_prep}} with a larger \code{max_D} value. The suggested
#' minimum is stored in \code{suggested_max_D}.
#'
#' When \code{sigma_precision} is triggered, the variable may not have a
#' meaningful scale of effect, or the data may be insufficient to estimate it
#' precisely. Interpret affected covariates with caution.
#'
#' @seealso \code{\link{multiScale_optim}}, \code{\link{kernel_prep}}
#'
#' @examples
#' \donttest{
#' data('pts')
#' data('count_data')
#' hab <- terra::rast(system.file('extdata', 'hab.tif', package = 'multiScaleR'))
#'
#' kernel_inputs <- kernel_prep(pts = pts,
#'                              raster_stack = hab,
#'                              max_D = 250,
#'                              kernel = 'gaussian')
#'
#' mod <- glm(y ~ hab, family = poisson, data = count_data)
#' opt <- multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs)
#'
#' ## Access diagnostics
#' diag <- diagnostics(opt)
#' diag$max_distance$triggered
#' diag$sigma_precision
#' }
#'
#' @export
diagnostics <- function(object, ...) {
  UseMethod("diagnostics")
}


#' @rdname diagnostics
#' @export
diagnostics.multiScaleR <- function(object, ...) {
  if (!is.null(object$diagnostics)) {
    return(object$diagnostics)
  }

  list(
    max_distance = NULL,
    sigma_precision = NULL,
    shape_precision = NULL
  )
}
