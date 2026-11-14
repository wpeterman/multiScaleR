#' Check Optional Software for Bayesian Scale Models
#'
#' Check whether the R packages that a [multiScale_bayes()] computation
#' method needs can be loaded and meet the package's minimum versions.
#' Nothing is installed, no model is fitted, and neither the compiler nor
#' sampling performance is tested.
#'
#' @param method Computation method to check, as in [multiScale_bayes()].
#'   `"sample"` (default, NUTS) needs RTMB, tmbstan, rstan, and posterior.
#'   `"laplace"` and `"marginal"` need RTMB. `"slice"` and `"prepare"` need
#'   no optional packages to run in one process. Slice chain diagnostics need
#'   posterior (version 1.6.0 or later).
#' @param correction `"none"` (default) or `"importance"`. Importance
#'   correction applies only to `method = "marginal"` and also needs
#'   posterior.
#' @param cores Requested number of worker processes; default 1. A value
#'   greater than one also checks processx, which parallel fits need. It has
#'   no effect for `"prepare"`, which always uses one process.
#'
#' @return An object of class `multiScaleR_bayes_setup`, a list with these
#'   elements.
#'   * `method`, `correction`: the settings checked.
#'   * `ready`: `TRUE` when every listed package loads and meets its minimum
#'     version.
#'   * `packages`: a data frame with one row per package and columns
#'     `package`, `minimum` (the minimum version, or `NA` when any version
#'     works), installed `version`, `loadable`, `meets_version`, `passed`, and
#'     a diagnostic `message` that says what to install or update.
#'   * `notes`: what the check does and does not cover.
#'
#'   Printing the object lists the packages and the notes, and returns it
#'   invisibly. A passed check does not certify a working compiler,
#'   successful sampling, convergence, or adequate memory. Fit a small model
#'   next.
#'
#' @details Preparation (`method = "prepare"`) validates your data and
#'   formulas before you install any sampling software. Installing the NUTS
#'   packages from source can require a C++ compiler, which on Windows comes
#'   with Rtools. `vignette("v13-bayesian-scales", package = "multiScaleR")`
#'   covers installation, and
#'   `vignette("v15-bayesian-computation", package = "multiScaleR")` compares
#'   the methods.
#'
#' @seealso [multiScale_bayes()]
#' @examples
#' check_bayes_setup(method = "prepare")  # always ready: nothing to install
#'
#' # Packages for slice sampling on two workers. `ready` is FALSE when
#' # processx is missing. The message column says what to install.
#' setup <- check_bayes_setup(method = "slice", cores = 2)
#' setup$ready
#' setup$packages[c("package", "loadable", "passed")]
#' @export
check_bayes_setup <- function(method = c("sample", "laplace", "marginal", "prepare", "slice"),
                              correction = c("none", "importance"), cores = 1L) {
  method <- match.arg(method); correction <- match.arg(correction)
  validate_scalar_numeric(cores, "cores", integerish = TRUE, lower = 1,
    upper = .Machine$integer.max)
  if (correction != "none" && method != "marginal")
    stop("Importance correction is available only with method = 'marginal'.", call. = FALSE)
  required <- switch(method, sample = c(RTMB = "2.0", tmbstan = "1.2.1",
    rstan = NA_character_, posterior = "1.6.0"), laplace = c(RTMB = "2.0"),
    marginal = c(RTMB = "2.0"), prepare = character(), slice = character())
  if (correction == "importance") required <- c(required, posterior = "1.6.0")
  if (method != "prepare" && cores > 1L) required <- c(required, processx = NA_character_)
  packages <- if (length(required)) do.call(rbind, lapply(seq_along(required), function(i)
    .msr_bayes_package_status(names(required)[i], required[i]))) else
      data.frame(package = character(), minimum = character(), version = character(),
        loadable = logical(), meets_version = logical(), passed = logical(), message = character())
  notes <- if (method == "slice")
    "Install posterior (>= 1.6.0) to compute chain diagnostics." else
      c("Checks R packages only; does not install software or run a model.",
        "Compiler configuration, memory and sampling diagnostics require a separate small fit.")
  structure(list(method = method, correction = correction, ready = all(packages$passed),
    packages = packages, notes = notes),
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
#' @param x A `multiScaleR_bayes_setup` object from `check_bayes_setup()`.
#' @param ... Not used. Supplying any further argument is an error.
#' @export
print.multiScaleR_bayes_setup <- function(x, ...) {
  if (length(list(...))) stop("Unused print arguments.", call. = FALSE)
  cat("Bayesian package setup:", x$method, "\n")
  cat(if (x$ready) "Package checks passed.\n" else "Packages need attention.\n")
  if (nrow(x$packages)) print(x$packages, row.names = FALSE)
  else if (x$method == "slice") cat("Slice sampling needs no optional fitting packages.\n")
  else cat("Preparation needs no optional Bayesian packages.\n")
  cat(paste(x$notes, collapse = "\n"), "\n")
  invisible(x)
}
