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
  analysis_mod <- .analysis_model(model)
  if (is.null(pkg) && !identical(analysis_mod, model)) {
    pkg <- extract_namespace(analysis_mod)
  }
  if(is.null(pkg)){
    pkg <- extract_call_function_package(model)
  }
  if (is.null(pkg) && !identical(analysis_mod, model)) {
    pkg <- extract_call_function_package(analysis_mod)
  }

  msr_path <- normalizePath(getNamespaceInfo("multiScaleR", "path"),
                            winslash = "/", mustWork = FALSE)
  msr_system_path <- normalizePath(system.file(package = "multiScaleR"),
                                   winslash = "/", mustWork = FALSE)
  msr_version <- as.character(utils::packageVersion("multiScaleR"))
  msr_load_all <- nzchar(msr_system_path) && !identical(msr_path, msr_system_path)
  lib_paths <- .libPaths()
  r_libs_user <- Sys.getenv("R_LIBS_USER", unset = "")
  r_libs <- Sys.getenv("R_LIBS", unset = "")

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

  # Export package metadata to workers
  clusterExport(cl,
                varlist = c("pkg", "msr_path", "msr_version", "msr_load_all",
                            "lib_paths", "r_libs_user", "r_libs"),
                envir = environment())

  # Load required packages on each worker. In development sessions, the master
  # process may be using pkgload/devtools from a source tree; PSOCK workers must
  # use that same source tree rather than silently loading an older installed
  # multiScaleR from .libPaths().
  worker_info <- clusterEvalQ(cl, {
    if (nzchar(r_libs_user)) {
      Sys.setenv(R_LIBS_USER = r_libs_user)
    }
    if (nzchar(r_libs)) {
      Sys.setenv(R_LIBS = r_libs)
    }
    .libPaths(unique(c(lib_paths, .libPaths())))

    if (!is.null(pkg) && !identical(pkg, "multiScaleR")) {
      library(pkg, character.only = TRUE)
    }

    if (isTRUE(msr_load_all)) {
      if (!requireNamespace("pkgload", quietly = TRUE)) {
        stop(
          paste0(
            "PSOCK workers need the `pkgload` package to load the same ",
            "development copy of multiScaleR as the main R session. ",
            "Install `pkgload`, reinstall multiScaleR, or use `n_cores = 1`."
          ),
          call. = FALSE
        )
      }
      pkgload::load_all(msr_path, quiet = TRUE)
    } else {
      library("multiScaleR")
    }

    list(
      version = as.character(utils::packageVersion("multiScaleR")),
      path = normalizePath(getNamespaceInfo("multiScaleR", "path"),
                           winslash = "/", mustWork = FALSE)
    )
  })

  worker_versions <- vapply(worker_info, `[[`, character(1), "version")
  if (any(worker_versions != msr_version)) {
    stop(
      paste0(
        "PSOCK workers loaded multiScaleR version(s) ",
        paste(unique(worker_versions), collapse = ", "),
        " but the main R session is using version ", msr_version, ". ",
        "Reinstall multiScaleR or use `n_cores = 1`."
      ),
      call. = FALSE
    )
  }

  if (isTRUE(msr_load_all)) {
    worker_paths <- vapply(worker_info, `[[`, character(1), "path")
    if (any(worker_paths != msr_path)) {
      stop(
        paste0(
          "PSOCK workers did not load the same development copy of multiScaleR ",
          "as the main R session. Use `n_cores = 1` or reinstall/load the ",
          "package consistently before running in parallel."
        ),
        call. = FALSE
      )
    }
  }

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


extract_call_function_package <- function(x) {
  fc <- tryCatch(
    get_call(x),
    error = function(e) NULL
  )

  if (is.null(fc) || !is.call(fc)) {
    return(NULL)
  }

  fit_fun <- fc[[1]]
  if (is.call(fit_fun) && identical(fit_fun[[1]], as.name("::"))) {
    return(as.character(fit_fun[[2]]))
  }

  if (!is.symbol(fit_fun)) {
    return(NULL)
  }

  fit_fun <- as.character(fit_fun)

  attached <- tryCatch(
    utils::find(fit_fun, mode = "function"),
    error = function(e) character()
  )
  attached <- attached[grepl("^package:", attached)]
  attached <- sub("^package:", "", attached)
  attached <- attached[attached %in% loadedNamespaces()]
  if (length(attached) > 0) {
    return(attached[[1]])
  }

  loaded <- loadedNamespaces()
  has_fun <- vapply(
    loaded,
    function(ns) {
      exists(fit_fun, envir = asNamespace(ns), mode = "function",
             inherits = FALSE)
    },
    logical(1)
  )
  loaded <- loaded[has_fun]
  if (length(loaded) > 0) {
    return(loaded[[1]])
  }

  NULL
}
