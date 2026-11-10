#' Covariance for Native GDR Response Fields
#'
#' Return retained joint covariance on the same physical coefficient scales as
#' `coef()`, without rerunning a TMB tape. Cross-component and fixed/field
#' covariance are retained. Native log-SD parameters remain in the diagnostic
#' packet and parameter penalty, rather than being added to regression output.
#' @param object A fitted [unmarkedFitGDRField-class].
#' @param type Optional active component. If omitted, return all regression
#'   components in native `coef()` order.
#' @param altNames Logical; default `TRUE` uses native component-prefixed names.
#' @param method Native interface control. Only `"hessian"` is qualified;
#'   bootstrap covariance requires a separately validated field bootstrap.
#' @param fixedOnly Logical; default `TRUE` returns regression covariance.
#'   `FALSE` also includes physical random-field/observer coefficients.
#' @param ... Reserved controls. Unrecognized arguments stop.
#' @return A named covariance matrix aligned with `coef(object, ...)`. Physical
#'   abundance states use their original Gaussian SD and fixed correlation.
#'   Full fixed/latent cross-terms are retained. This covariance holds ecological
#'   scales and correlation range fixed and adds no observation noise.
#' @importMethodsFrom unmarked vcov
#' @md
#' @export
methods::setMethod("vcov", "unmarkedFitGDRField", function(object, type,
    altNames = TRUE, method = "hessian", fixedOnly = TRUE, ...) {
  field <- .msr_gdr_field_guard(object)
  validate_scalar_logical(altNames, "altNames")
  validate_scalar_logical(fixedOnly, "fixedOnly")
  if (!identical(method, "hessian") || length(list(...)))
    stop("Only retained joint GDR field covariance is qualified; field bootstrap needs separate validation.", call. = FALSE)
  components <- if (missing(type)) names(field$bases) else type
  if ((!missing(type) && length(type) != 1L) || !is.character(components) || anyNA(components) ||
      any(!components %in% names(field$bases)) || anyDuplicated(components))
    stop("Supply one active GDR regression component or omit type.", call. = FALSE)
  point <- field$point
  indices <- unlist(lapply(components, function(name)
    c(which(names(point) == paste0("beta_", name)),
      if (!fixedOnly) which(names(point) == paste0("b_", name)))), use.names = FALSE)
  labels <- if (missing(type)) names(unmarked::coef(object, altNames = altNames, fixedOnly = fixedOnly)) else
    names(unmarked::coef(object, type = type, altNames = altNames, fixedOnly = fixedOnly))
  if (length(indices) != length(labels)) stop("Native coefficient covariance labels changed.", call. = FALSE)
  G <- matrix(0, length(indices), length(point))
  G[cbind(seq_along(indices), indices)] <- 1
  lambda <- which(names(point) == "b_lambda")
  rows <- match(lambda, indices)
  if (any(!is.na(rows))) G[rows[!is.na(rows)], lambda] <- field$lower[!is.na(rows), , drop = FALSE]
  covariance <- G %*% field$joint_covariance %*% t(G)
  dimnames(covariance) <- list(labels, labels)
  covariance
})

#' @rdname unmarkedFitGDRField-class
#' @param parm,type,level,method Native interval controls. Only normal intervals
#'   using the retained component covariance are qualified for this field.
#' @importMethodsFrom unmarked confint
#' @export
methods::setMethod("confint", "unmarkedFitGDRField", function(object, parm,
    level = .95, type, method = c("normal", "profile")) {
  .msr_gdr_field_guard(object)
  method <- match.arg(method)
  if (method != "normal")
    stop("Spatial GDR profile intervals require separate field qualification.", call. = FALSE)
  methods::callNextMethod()
})
