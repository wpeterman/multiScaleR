#' Covariance Matrix for a Distance-Removal Model with a Fixed-Range Field
#'
#' Return the covariance matrix of the coefficients of a model fitted by
#' [fit_gdr_response_field()]. The matrix is retained from the fit, so this
#' method does not rerun TMB. It matches `coef()` on the link scale and keeps
#' the covariances between components and between the fixed coefficients and
#' the field values.
#'
#' @param object An [unmarkedFitGDRField][unmarkedFitGDRField-class] returned
#'   by [fit_gdr_response_field()].
#' @param type Optional character naming one component: `"lambda"`
#'   (abundance), `"dist"` (distance) or `"rem"` (removal). If omitted, all
#'   components are returned together, in the order of `coef(object)`.
#' @param altNames Logical, default `TRUE`. `TRUE` labels rows and columns
#'   with component prefixes, such as `lambda(Int)`, as `coef()` does.
#' @param method Character. Only `"hessian"`, the default, is available. It
#'   returns the covariance computed from the curvature of the likelihood
#'   during fitting. Bootstrap covariance is not available for this model.
#' @param fixedOnly Logical, default `TRUE`. `TRUE` returns the fixed
#'   regression coefficients only. `FALSE` adds the field value of every block
#'   (labeled `lambda(b_lambda)`) and any observer intercepts in the distance
#'   component.
#' @param ... Not used. Unrecognized arguments stop the call.
#' @return A square numeric matrix with row and column names matching
#'   `coef(object, type, altNames = altNames, fixedOnly = fixedOnly)`. Entries
#'   are on the link scale of each coefficient: log for abundance and
#'   distance, logit for removal. Block field values are physical deviations
#'   on the log density scale, not TMB's internal whitened values. The square
#'   roots of the diagonal are the standard errors.
#' @details
#' The log SDs of the random effects are estimated parameters, but they are
#' not rows of this matrix. Use `unmarked::sigma(object)` to see the field SD
#' with its confidence limits.
#'
#' The matrix treats the ecological scales of effect and the correlation
#' range as known, and it contains no observation noise. `confint()` uses it
#' for normal (Wald) intervals, and
#' [predict()][predict,unmarkedFitGDRField-method] uses it for standard
#' errors.
#' @seealso [fit_gdr_response_field()], [unmarkedFitGDRField-class],
#'   [predict,unmarkedFitGDRField-method]
#' @examples
#' \donttest{
#' # The field needs unmarked 1.5.2 and TMB 1.9.25.
#' if (requireNamespace("TMB", quietly = TRUE) &&
#'     utils::packageVersion("unmarked") == "1.5.2" &&
#'     utils::packageVersion("TMB") == "1.9.25") {
#'   # Eight survey blocks 1 km apart, six point counts in each
#'   block <- factor(rep(paste0("b", 1:8), each = 6))
#'   hab <- rep(seq(-1, 1, length.out = 6), 8)
#'   shift <- c(-0.4, -0.3, -0.1, 0, 0.1, 0.2, 0.35, 0.4)[block]
#'   sim <- sim_dat_gdr(density = exp(log(400) + 0.5 * hab + shift),
#'     distance_scale = 40, removal = 0.2, dist_breaks = c(0, 25, 50, 75),
#'     interval_lengths = c(2, 3, 5), interval_units = "minutes",
#'     site_covs = data.frame(hab = hab, block = block), seed = 1)
#'   template <- unmarked::gdistremoval(lambdaformula = ~ hab + (1 | block),
#'     removalformula = ~ 1, distanceformula = ~ 1, data = sim$observed,
#'     keyfun = "halfnorm", mixture = "P", output = "density",
#'     unitsOut = "kmsq", K = 60, engine = "TMB",
#'     starts = c(log(400), 0, log(40), qlogis(0.2), log(0.3)))
#'   block_xy <- cbind(x = 500000 + 1000 * (as.integer(block) - 1),
#'                     y = 4400000)
#'   rownames(block_xy) <- rownames(unmarked::siteCovs(sim$observed))
#'   cold <- c("beta_lambda[1]" = log(300), "beta_lambda[2]" = 0,
#'             "beta_dist[1]" = log(30), "beta_rem[1]" = qlogis(0.4),
#'             "lsigma_lambda[1]" = log(0.5))
#'   settings <- list(method = "BFGS", se = TRUE, threads = 1L,
#'                    control = list(maxit = 2500L, reltol = 1e-10))
#'   field_fit <- fit_gdr_response_field(template, block_xy, crs = 32617,
#'     range = 2000, cold_starts = list(cold = cold),
#'     audit_dir = tempfile("gdr-field-"), tmb_settings = settings)
#'
#'   # Abundance coefficients: covariance and standard errors (log scale)
#'   V <- vcov(field_fit, type = "lambda")
#'   print(V)
#'   print(sqrt(diag(V)))
#'   # All components together: abundance, distance, then removal
#'   print(round(vcov(field_fit), 4))
#'   # Adding the eight block field values gives a 10 x 10 matrix
#'   print(dim(vcov(field_fit, type = "lambda", fixedOnly = FALSE)))
#' }
#' }
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
#' @param parm,type,level,method Arguments of `confint()`. `type` names the
#'   component (`"lambda"`, `"dist"` or `"rem"`), `parm` optionally selects
#'   coefficients within it, and `level` is the confidence level (default
#'   0.95). `method = "normal"` (the default) gives Wald intervals from the
#'   retained covariance. `method = "profile"` stops for this class.
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
