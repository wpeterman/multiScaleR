#' Check Optional Software for Bayesian Scale Models
#'
#' Check whether the R packages needed by a Bayesian inference method can be
#' loaded and meet the package's minimum versions. Nothing is installed, no
#' model is fitted, and compiler or sampling performance is not tested.
#'
#' @param method Inference method to check: \code{"sample"} (default) needs
#'   RTMB, tmbstan, rstan, and posterior; \code{"laplace"} and
#'   \code{"marginal"} need RTMB; \code{"prepare"} needs none of these.
#' @param correction \code{"none"} (default) or \code{"importance"} for
#'   \code{method = "marginal"}. Importance correction also needs posterior.
#' @return A \code{multiScaleR_bayes_setup} list with \code{method},
#'   \code{correction}, \code{ready} (whether package checks passed),
#'   \code{packages} (a data frame with \code{package}, \code{minimum},
#'   installed \code{version}, \code{loadable}, \code{meets_version},
#'   \code{passed}, and a diagnostic \code{message}), and \code{notes}.
#'   A passed check does not certify a working compiler, successful sampling,
#'   convergence, or adequate memory. Fit a small model next.
#' @details Preparation can validate your data and formulas before installing
#'   sampling software. For installation instructions, including Windows
#'   toolchain guidance, see
#'   \code{vignette("bayesian-scales", package = "multiScaleR")}.
#' @seealso \code{\link{multiScale_bayes}}
#' @examples
#' check_bayes_setup(method = "prepare")
#' @export
check_bayes_setup <- function(method = c("sample", "laplace", "marginal", "prepare"),
                              correction = c("none", "importance")) {
  method <- match.arg(method); correction <- match.arg(correction)
  if (correction != "none" && method != "marginal")
    stop("Importance correction is available only with method = 'marginal'.", call. = FALSE)
  required <- switch(method, sample = c(RTMB = "2.0", tmbstan = "1.2.1",
    rstan = NA_character_, posterior = "1.6.0"), laplace = c(RTMB = "2.0"),
    marginal = c(RTMB = "2.0"), prepare = character())
  if (correction == "importance") required <- c(required, posterior = "1.6.0")
  packages <- if (length(required)) do.call(rbind, lapply(seq_along(required), function(i)
    .msr_bayes_package_status(names(required)[i], required[i]))) else
      data.frame(package = character(), minimum = character(), version = character(),
        loadable = logical(), meets_version = logical(), passed = logical(), message = character())
  structure(list(method = method, correction = correction, ready = all(packages$passed),
    packages = packages, notes = c("Checks R packages only; does not install software or run a model.",
      "Compiler configuration, memory and sampling diagnostics require a separate small fit.")),
    class = "multiScaleR_bayes_setup")
}

.msr_bayes_package_status <- function(package, minimum) {
  message <- ""
  loadable <- tryCatch(requireNamespace(package, quietly = TRUE), error = function(e) {
    message <<- conditionMessage(e); FALSE
  })
  version <- suppressWarnings(tryCatch(as.character(if (loadable) getNamespaceVersion(package) else
    utils::packageVersion(package)), error = function(e) NA_character_))
  meets <- !is.na(version) && (is.na(minimum) || utils::compareVersion(version, minimum) >= 0L)
  if (!loadable && !nzchar(message)) message <- "Package is missing or cannot be loaded."
  if (loadable && !meets) message <- paste("Update to version", minimum, "or later.")
  data.frame(package = package, minimum = unname(minimum), version = version,
    loadable = loadable, meets_version = meets, passed = loadable && meets, message = message)
}

#' @rdname check_bayes_setup
#' @param x A setup-check object.
#' @param ... Reserved for extensions; unsupported arguments cause an error.
#' @export
print.multiScaleR_bayes_setup <- function(x, ...) {
  if (length(list(...))) stop("Unused print arguments.", call. = FALSE)
  cat("Bayesian package setup:", x$method, "\n")
  cat(if (x$ready) "Package checks passed.\n" else "Packages need attention.\n")
  if (nrow(x$packages)) print(x$packages, row.names = FALSE)
  else cat("Preparation needs no optional Bayesian packages.\n")
  cat(paste(x$notes, collapse = "\n"), "\n")
  invisible(x)
}
