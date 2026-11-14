#' List the Starting-Value Order of a Distance-Removal Model
#'
#' Read the fitted coefficient and log SD values of an `unmarked`
#' distance-removal (GDR) model in the order that `unmarked::gdistremoval()`
#' expects for its `starts` argument, with a map that names the model
#' component and position of each value. Use it to build extra starting
#' values for [make_gdr_refitter()] and [fit_gdr_response_field()].
#'
#' @param model An `unmarkedFitGDR` fitted by `unmarked::gdistremoval()` with
#'   `unmarked` 1.5.2, using the C or TMB engine. Its stored estimates and
#'   optimizer values must be complete and agree with each other. Models with
#'   `offset()` terms and zero-inflated TMB models are not supported. Finite
#'   values do not show that the model converged.
#' @return A plain list:
#'   * `starts`: an unnamed numeric vector of the fitted values in starting
#'     order, on the link scales (log, logit, and log for SDs). Random-effect
#'     values are not included.
#'   * `map`: a data frame with one row per element of `starts`:
#'     `constructor_position`, the TMB parameter `block` (such as
#'     `"beta_lambda"` or `"lsigma_lambda"`), `within_block` (position in the
#'     block), `fitted_position` (position in the fitted optimizer vector),
#'     `component` (`"lambda"`, `"dist"`, `"rem"` and so on), `kind`
#'     (`"regression"` or `"log_sd"`) and `label`, the coefficient or variance
#'     name stored by `unmarked`. For a TMB-engine model, paste `block` and
#'     `within_block` as `"beta_lambda[1]"` to name a start for
#'     [make_gdr_refitter()].
#'   * `fitted_values`: the fitted optimizer vector with its own names and
#'     order. `coordinate_count` is its length, not a count for AIC.
#'     `random_modes_omitted` gives the number of random-effect values left
#'     out in each component.
#'   * `optimizer`: the stored objective `value` and `convergence` code, as
#'     reported by the fit.
#'   * `provider`: the `unmarked` version, the engine and, for TMB, the TMB
#'     version. `format_version` identifies the layout of this list.
#'   * `fingerprint`: a SHA-256 of this list. It does not identify the data,
#'     settings or likelihood.
#'   * `fit_eligibility_established`, `likelihood_evaluated` and
#'     `uncertainty_evaluated`: always `FALSE`. The function only reads stored
#'     values.
#' @details
#' # Why the order matters
#'
#' With the TMB engine, the order of the `starts` argument can differ from
#' the order of the fitted optimizer vector. For example, `unmarked` may
#' place the negative binomial dispersion after the regression and variance
#' blocks in the starts, but next to abundance in the fitted vector. Repeated
#' coefficient labels cannot resolve this, so the map identifies each value
#' by its block and position and checks every value against the stored
#' estimates.
#'
#' # What the function does not do
#'
#' It reads stored values only. It does not evaluate the likelihood, its
#' gradient or Hessian, and it leaves the model unchanged. A failed or
#' nonconverged fit can still have finite values, and its stored diagnostics
#' keep their meaning. Nonfinite or inconsistent values stop with an error
#' rather than being dropped or repaired.
#'
#' The list can be saved with `saveRDS()` without saving a fitted model. It
#' cannot restore a fit. Refitting still needs the original data,
#' transformations, formulas and optimizer settings.
#' @seealso [make_gdr_refitter()], [fit_gdr_response_field()],
#'   [expand_gdr_intervals()]
#' @examples
#' frame <- unmarked::unmarkedFrameGDR(matrix(c(0, 0, 1, 0), 2, byrow = TRUE),
#'   matrix(c(0, 0, 0, 1, 0, 0), 2, byrow = TRUE),
#'   dist.breaks = c(0, 25, 75), unitsIn = "m", period.lengths = c(1, 2, 1))
#' # A fixed point for inspecting coordinates, not a converged fit.
#' point <- unmarked::gdistremoval(data = frame, engine = "C", K = 30,
#'   starts = c(log(8), log(40), qlogis(.3)), se = FALSE,
#'   control = list(maxit = 0))
#' snapshot <- gdr_start_values(point)
#' # Log density, log distance scale and logit removal, in starting order
#' snapshot$starts
#' # Component, block position and stored label of each value
#' snapshot$map[, c("block", "within_block", "component", "kind", "label")]
#' @export
gdr_start_values <- function(model) {
  fail <- function(message) stop(message, call. = FALSE)
  if (!methods::is(model, "unmarkedFitGDR"))
    fail("model must contain the qualified native GDR coordinate metadata.")
  if (!methods::is(model@data, "unmarkedFrameGDR") ||
      length(model@mixture) != 1L || !model@mixture %in% c("P", "NB", "ZIP") ||
      length(model@keyfun) != 1L || !model@keyfun %in% c("halfnorm", "exp", "hazard", "uniform"))
    fail("Stored native GDR mixture, key and frame metadata must remain qualified.")
  .msr_gdr_native_guard(model)
  version <- as.character(utils::packageVersion("unmarked"))
  if (version != "1.5.2") fail("This native GDR coordinate layout requires separate provider qualification.")
  expected <- c("lambda", if (model@mixture != "P") "alpha",
    if (model@data@numPrimary > 1L) "phi", if (model@keyfun != "uniform") "dist",
    if (model@keyfun == "hazard") "scale", "rem")
  components <- model@estimates@estimates
  if (!identical(names(components), expected))
    fail("Stored native components do not match the complete active GDR layout.")
  beta <- sigma <- labels_beta <- labels_sigma <- list()
  omitted <- stats::setNames(integer(length(expected)), expected)
  valid_labels <- function(labels, values) is.character(labels) &&
    length(labels) == length(values) && !anyNA(labels) && all(nzchar(labels))
  for (name in expected) {
    estimate <- components[[name]]
    fixed <- estimate@fixed
    values <- estimate@estimates
    if (!is.numeric(fixed) || !length(fixed) || anyNA(fixed) ||
        !identical(as.numeric(fixed), as.numeric(seq_along(fixed))) ||
        length(values) < length(fixed))
      fail("Stored native fixed coordinates are incomplete or noncanonical.")
    beta[[name]] <- unname(values[fixed])
    labels_beta[[name]] <- names(values)[fixed]
    random <- estimate@randomVarInfo
    sigma[[name]] <- if (length(random)) unname(random$estimates) else numeric()
    labels_sigma[[name]] <- if (length(random)) random$names else character()
    if ((!is.numeric(sigma[[name]])) ||
        (length(random) && (is.null(random$estimates) || is.null(random$names))))
      fail("Stored native variance coordinates are incomplete.")
    if (!valid_labels(labels_beta[[name]], beta[[name]]) ||
        !valid_labels(labels_sigma[[name]], sigma[[name]]))
      fail("Original native coefficient and variance labels must be retained.")
    omitted[name] <- length(values) - length(fixed)
  }
  blocks <- c(rep(paste0("beta_", expected), lengths(beta)),
    rep(paste0("lsigma_", expected), lengths(sigma)))
  reference <- c(unlist(beta, use.names = FALSE), unlist(sigma, use.names = FALSE))
  labels <- c(unlist(labels_beta, use.names = FALSE), unlist(labels_sigma, use.names = FALSE))
  values <- model@opt$par
  tmb <- !is.null(model@TMB)
  if (!tmb && (any(lengths(sigma)) || any(omitted)))
    fail("The native C coordinate layout cannot contain random or variance blocks.")
  fitted_blocks <- if (tmb) names(values) else blocks
  incoming <- if (tmb) model@TMB$starts_order else blocks
  checked <- .msr_gdr_exact_start_map(blocks, fitted_blocks, values)
  if (!isTRUE(all.equal(checked$starts, reference, tolerance = 0, check.attributes = FALSE)))
    fail("Stored native component values and optimizer coordinates disagree.")
  result <- .msr_gdr_exact_start_map(incoming, fitted_blocks, values)
  native_to_input <- .msr_gdr_exact_start_map(incoming, blocks, reference)$map$fitted_position
  result$map$component <- sub("^(beta|lsigma)_", "", result$map$block)
  result$map$kind <- ifelse(startsWith(result$map$block, "lsigma_"), "log_sd", "regression")
  result$map$label <- labels[native_to_input]
  optimizer <- model@opt[intersect(c("value", "convergence"), names(model@opt))]
  if (any(!vapply(optimizer, function(x) is.numeric(x) && length(x) == 1L && is.null(dim(x)), logical(1))))
    fail("Stored optimizer diagnostics must retain their original scalar numeric form.")
  result$fitted_values <- values
  result$coordinate_count <- length(values)
  result$random_modes_omitted <- omitted
  result$optimizer <- optimizer
  result$provider <- list(unmarked = version, engine = if (tmb) "TMB" else "C",
    TMB = if (tmb) as.character(utils::packageVersion("TMB")) else NULL)
  result$format_version <- 1L
  result$fit_eligibility_established <- FALSE
  result$likelihood_evaluated <- FALSE
  result$uncertainty_evaluated <- FALSE
  result$fingerprint <- .msr_gdr_hash(result)
  result
}
