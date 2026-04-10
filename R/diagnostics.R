#' Retrieve diagnostics from multiScaleR objects
#'
#' Returns structured warning/diagnostic information stored on fitted
#' \code{multiScaleR} objects.
#'
#' @param object An object to inspect.
#' @param ... Additional arguments passed to methods.
#'
#' @return For \code{multiScaleR} objects, a named list of diagnostics.
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
