#' Inspect Native GDR Constructor Start Coordinates
#'
#' Extract stored numeric optimizer values in the original constructor order,
#' retaining a reversible component and within-component position map.
#' @param model An `unmarkedFitGDR` from the qualified unmarked 1.5.2 layout.
#'   Stored fixed coefficients, log-SDs, optimizer values and constructor block
#'   metadata must be complete and consistent. Existing offset and TMB ZIP
#'   restrictions remain in force. A finite coordinate snapshot does not
#'   establish that the model converged or is eligible for inference.
#' @return A plain list with:
#'   \itemize{
#'     \item `starts`: unnamed numeric constructor starts on the native link
#'       scales. Log-SDs remain log-SDs; latent random modes are not included.
#'     \item `map`: `constructor_position`, native `block`, `within_block`,
#'       `fitted_position`, `component`, `kind` (`"regression"` or `"log_sd"`),
#'       and the original stored coefficient or variance `label`.
#'     \item `fitted_values`: the original fixed optimizer vector, including
#'       its names and order; `coordinate_count`: its length, not a new AIC
#'       parameter-count estimate; `random_modes_omitted`: the stored number
#'       of latent coordinates in each native component.
#'     \item `optimizer`: stored scalar `value` and `convergence` when present.
#'       Nonfinite diagnostic values are retained; they do not validate a fit.
#'     \item `provider`: unmarked version, native engine and, for TMB, the
#'       installed TMB version. `format_version` identifies the snapshot schema.
#'     \item `fingerprint`: SHA256 of this coordinate snapshot. It does not
#'       identify the complete data, likelihood, controls or prediction context.
#'     \item `qualification`: coordinate-only scope.
#'       `fit_eligibility_established`, `likelihood_evaluated` and
#'       `uncertainty_evaluated` are all `FALSE`.
#'   }
#' @details
#' Native TMB incoming starts and fitted optimizer coordinates can have different
#'   orders. For example, the constructor may place NB dispersion after the
#'   regression and variance blocks, while the fitted vector places it beside
#'   abundance. Repeated coefficient labels do not resolve this difference.
#'   The map uses the native block and its within-block position, checking every
#'   coordinate against the stored component values before returning starts.
#'
#' This function reads stored metadata only. It does not call a saved likelihood,
#'   update, score, Hessian, covariance or TMB tape, and it leaves the fitted model
#'   unchanged. A failed or nonconverged fit can still have finite coordinates;
#'   its stored diagnostics retain their original interpretation. Nonfinite or
#'   inconsistent coordinates stop rather than being dropped or repaired.
#'
#' The list can be saved with `saveRDS()` without serializing a fitted model or
#'   executable objective. Restoring an objective still requires the original
#'   data, transformations, formulas and evaluated controls, plus all original
#'   numerical and scientific checks. The snapshot is not a likelihood cache or
#'   a portable fitted model. It does not restore latent modes or qualify native
#'   TMB likelihood-only fitting; top-level `se = FALSE` and optimizer controls
#'   do not establish that behavior in the inspected provider.
#'
#' Priors, normalized targets, parameter penalties, diagnostics, native
#'   coefficient covariance and uncertainty distinctions are unchanged.
#' @seealso [make_gdr_refitter()], [expand_gdr_intervals()]
#' @examples
#' frame <- unmarked::unmarkedFrameGDR(matrix(c(0, 0, 1, 0), 2, byrow = TRUE),
#'   matrix(c(0, 0, 0, 1, 0, 0), 2, byrow = TRUE),
#'   dist.breaks = c(0, 25, 75), unitsIn = "m", period.lengths = c(1, 2, 1))
#' # A fixed point for inspecting coordinates, not a converged fit.
#' point <- unmarked::gdistremoval(data = frame, engine = "C", K = 30,
#'   starts = c(log(8), log(40), qlogis(.3)), se = FALSE,
#'   control = list(maxit = 0))
#' snapshot <- gdr_start_values(point)
#' snapshot$starts
#' snapshot$map
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
