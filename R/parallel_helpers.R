#' Identify the originating package of an S3 or S4 model object
#' and load it on a PSOCK cluster
#'
#' @param model An R model object (S3 or S4)
#' @param cl A PSOCK cluster (created with parallel::makeCluster)
#'
#' @return Invisibly returns the package name loaded
#' @details For internal use
#' @rdname cluster_prep
#' @keywords internal
#' @importFrom utils getS3method
#' @importFrom methods getClassDef
cluster_prep <- function(model, cl) {

  pkg <- extract_namespace(model)

  if(is.null(pkg)){
    model_classes <- class(model)

    # Try S3 generic method (e.g., predict.lm)
    pkg <- NULL
    for (cls in model_classes) {
      method <- tryCatch(
        getS3method("predict", cls, optional = TRUE),
        error = function(e) NULL
      )
      if (!is.null(method)) {
        env <- environment(method)
        if (!is.null(env) && !is.null(pkgname <- environmentName(env)) &&
            pkgname %in% loadedNamespaces()) {
          pkg <- pkgname
          break
        }
      }
    }

    # Try S4 class definition
    if (is.null(pkg) && isS4(model)) {
      for (cls in model_classes) {
        s4class <- tryCatch(getClassDef(cls), error = function(e) NULL)
        if (!is.null(s4class) && !is.null(pkgname <- s4class@package)) {
          pkg <- pkgname
          break
        }
      }
    }

    if (is.null(pkg)) {
      stop("Could not determine the package for the model object.")
    }
  }

  # Export pkg to workers
  clusterExport(cl, varlist = "pkg", envir = environment())

  # Load required packages on each worker
  clusterEvalQ(cl, {
    library(pkg, character.only = TRUE)
    library("multiScaleR")
  })

  invisible(pkg)
}




#' Extract the Namespace from a Model Call
#'
#' This function attempts to extract the package namespace from the call
#' used to fit a model object, assuming the function was called using the
#' `pkg::fun()` syntax.
#'
#' @param x A fitted model object.
#'
#' @return A character string with the namespace (package name), or `NULL`
#'   if the namespace cannot be determined.
#'
#' @details For internal use
#' @rdname extract_namespace
#' @keywords internal
#' @importFrom insight get_call
extract_namespace <- function(x) {
  fc <- tryCatch(
    get_call(x),
    error = function(e) NULL
  )
  if (!is.null(fc)) {
    fc <- paste(fc, collapse = "")  # ensure single string
    if (grepl("::", fc)) sub("::.*", "", fc) else NULL
  }
}
