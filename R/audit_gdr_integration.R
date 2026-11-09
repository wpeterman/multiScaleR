# Reproduce the qualified provider's distance arithmetic for diagnosis, not
# fitting. Keep its cancellation/quadrature errors visible beside the stable law.
.msr_gdr_audit_native_bins <- function(scale, breaks, key, area, weights) {
  if (!is.matrix(area) || !is.matrix(weights) || !identical(dim(area), dim(weights)) ||
      nrow(area) != length(scale) || ncol(area) != length(breaks) - 1L ||
      any(!is.finite(area)) || any(area <= 0) || any(!is.finite(weights)) || any(weights < 0))
    stop("Native distance geometry is not aligned with the audit rows.")
  out <- area
  rates <- unique(scale)
  rate_row <- match(scale, rates)
  for (j in seq_len(ncol(area))) {
    if (key == "uniform") {
      integral <- rep(1, length(rates))
    } else if (key == "halfnorm") {
      s2 <- rates^2
      p1 <- 1 - exp(-breaks[j + 1L]^2 / (2 * s2))
      p2 <- 1 - exp(-breaks[j]^2 / (2 * s2))
      integral <- (s2 * p1 - s2 * p2) * 2 * pi
    } else {
      # unmarked1.5.2 uses100 trapezoids per point-transect exponential bin.
      # Evaluate each distinct scale once, retaining the native summation order.
      h <- (breaks[j + 1L] - breaks[j]) / 100
      f <- function(x) exp(-x / rates) * x
      total <- numeric(length(rates))
      for (k in 1:99) total <- total + f(breaks[j] + k * h)
      integral <- h / 2 * (f(breaks[j]) + 2 * total + f(breaks[j + 1L])) * 2 * pi
    }
    detection <- if (key == "uniform") rep(1, nrow(area)) else integral[rate_row] / area[, j]
    out[, j] <- detection * weights[, j]
  }
  out
}

#' Audit native GDR integration at the stored coefficient estimates
#'
#' Separate finite latent-count truncation from distance-integration differences
#' without refitting the model or changing its likelihood.
#' @param model A native C `unmarkedFitGDR` from unmarked 1.5.2, with one primary
#'   period, complete paired margins, fixed effects, Poisson or NB abundance,
#'   and a half-normal, exponential or uniform point-transect distance key.
#'   Active offsets and other native domains retain their separate guards.
#' @param cutoffs Optional finite positive integer vector of additional native
#'   latent-count cutoffs. Values must be unique, at least the largest observed
#'   total and at most 100000. The stored cutoff is always included first.
#'   At most 32 distinct cutoffs and ten million event-by-latent-count positions
#'   per evaluation are admitted. These bounds limit diagnostic work, not model
#'   parameters.
#' @param design Logical. If `TRUE`, add compact rank and conditioning summaries
#'   of the stored native component matrices. Reconstructed state spline recipes
#'   are reported only after exact agreement with the stored native state matrix.
#'   Default: `FALSE`, preserving the original audit output.
#' @param design_data Optional data frame of state-design scenarios, in the
#'   model's retained predictor units. Requires `design = TRUE`. This evaluates
#'   the verified training basis, not ecological kernels, response predictions
#'   or new observation/random-effect contexts.
#' @return A list with the `cutoffs` table, `references` for the infinite-sum
#'   native-bin and stable analytic laws, `distance` comparison, stored-point
#'   agreement, coefficient and provider fingerprints, and explicit numerical
#'   qualification flags. The additive `finite_poisson` reference contains its
#'   status, compact evaluations and a cutoff comparison with native nll values,
#'   differences and agreement flags. It does not supply fitting admission.
#'   Errors and warnings from individual numerical
#'   evaluations remain in the corresponding records. No model, optimizer,
#'   covariance or likelihood environment is stored in the returned audit.
#'   With `design = TRUE`, `design` contains native matrix summaries, requested
#'   spline calls, retained knots/boundaries and optional scenario support counts.
#'   Reconstruction errors and native warnings remain explicit. No design
#'   matrices, formula environments or individual scenario rows are returned.
#' @details
#' Original evaluated designs, area conversions, units and durations come from
#'   the stored C likelihood. Additional evaluations use those same coefficients
#'   and one thread, with only the audit cutoff changed. The original model and
#'   its cutoff, parameter count, coefficient uncertainty and ecological scales
#'   are unchanged. No optimization, Hessian or `sdreport()` is performed.
#'
#' `truncation_loss` is native finite-cutoff nll minus the infinite-sum nll
#'   using the same native distance bins. Qualification requires a finite loss
#'   between -1e-11 and 1e-9. This is evaluated at each supplied cutoff; it is
#'   not an adaptive change to an ecological-scale optimization target.
#'
#' For Poisson abundance, the additional finite-cutoff reference uses the exact
#' identity `Pois(n; lambda*q) * P(Pois(lambda*(1-q)) <= K-n)` at each event,
#' followed by the same native-bin distance/removal multinomial terms. It avoids
#' a latent-abundance loop for this reference and preserves the original finite
#' cutoff. The cutoff is not renormalized into a truncated abundance model.
#' The original native evaluations and truncation-loss gate still run unchanged.
#' `finite_poisson$comparison$agrees` requires finite values, no evaluation
#' warnings and an absolute nll difference at most 1e-9. Errors stay in its
#' evaluation records. NB retains its original audit and is explicitly marked
#' unqualified for this additional Poisson identity. Expected means retain
#' their underlying native meaning; no uncertainty or TMB integration is added.
#'
#' The separate distance comparison uses the infinite-sum native-bin law and
#'   stable analytic point-annulus integrals. Native exponential bins use 100
#'   trapezoids per bin, so the analytic exponential target generally differs.
#'   Native half-normal arithmetic can lose precision near the uniform limit.
#'   Invalid native probabilities are retained as failed reference records;
#'   probabilities are never clamped or renormalized to pass qualification.
#'
#' Numerical agreement does not establish convergence, global optimality,
#'   biological fit, density identification, or calibrated uncertainty. The
#'   observed count law, expected circle abundance, density and realized latent
#'   population are distinct quantities. TMB/random-effect/ZIP/missing-history
#'   and multiple-primary extensions require separate qualification.
#'
#' Optional design diagnostics use native columns without additional rescaling.
#' Rank uses QR tolerance `1e-7`; this descriptive tolerance does not change any
#' fitting or numerical-admission gate. The condition number is infinite when
#' column rank is deficient.
#'
#' At most 256 columns and two million matrix positions per component are
#' summarized. Larger designs retain an explicit work-limit status; scenario
#' work is bounded before its model matrix is allocated. These are read-only
#' audit limits and do not constrain or qualify a fitted model.
#'
#' Basis reconstruction must reproduce the original
#' matrix exactly, including factor/interacted columns, before native training
#' prediction calls can be used. No knot is repaired or moved by the audit.
#' Standard scale/poly/ns/bs preprocessing uses the existing training-only
#' machinery shared with CV. Unqualified custom transformations fail the
#' scenario audit explicitly rather than being learned from scenario rows.
#' Requested calls and actual retained knots are distinct because native spline
#' constructors may adjust a knot that coincides with a boundary.
#'
#' Support counts are provided for bare numeric spline inputs. Other expressions
#' keep their native matrix/recipe audit with an explicit support status.
#' Support calculations do not reevaluate those source expressions. Counts
#' outside observed training ranges and outside spline boundaries describe
#' different extrapolation questions. A full-rank basis does not establish a
#' converged fit, ecological threshold, response shape or calibrated uncertainty.
#' @seealso [make_gdr_refitter()], [validate_unmarked_data()]
#' @export
audit_gdr_integration <- function(model, cutoffs = NULL, design = FALSE, design_data = NULL) {
  validate_scalar_logical(design, "design")
  if (!is.null(design_data) && (!isTRUE(design) || !is.data.frame(design_data) ||
      !nrow(design_data) || anyDuplicated(names(design_data))))
    stop("design_data requires design = TRUE and a nonempty data frame with unique columns.", call. = FALSE)
  settings <- .msr_gdr_C_settings(model)
  frame <- model@data
  if (frame@numPrimary != 1L || !model@mixture %in% c("P", "NB") ||
      !model@keyfun %in% c("halfnorm", "exp", "uniform") ||
      !identical(frame@survey, "point") || !frame@unitsIn %in% c("m", "km"))
    stop("This audit requires single-primary complete fixed-effect C P/NB point-transect models with qualified distance keys.")
  env <- environment(model@nllFun)
  binding <- function(name) get(name, envir = env, inherits = FALSE)
  gd <- binding("gd")
  yd <- frame@yDistance
  yr <- frame@yRemoval
  complete <- function(y) is.matrix(y) && nrow(y) > 0L && ncol(y) > 0L &&
    all(is.finite(y)) && all(y >= 0) && all(y == floor(y))
  if (!complete(yd) || !complete(yr) || nrow(yd) != nrow(yr) ||
      any(rowSums(yd) != rowSums(yr)) || !identical(gd$yDist, as.vector(t(yd))) || !identical(gd$yRem, as.vector(t(yr))))
    stop("The audit requires complete paired native response margins without dropped or changed rows.")
  db <- frame@dist.breaks
  pl <- frame@period.lengths
  if (length(db) != ncol(yd) + 1L || any(!is.finite(db)) || db[1L] != 0 || any(diff(db) <= 0) ||
      length(pl) != ncol(yr) || any(!is.finite(pl)) || any(pl <= 0) || any(pl != floor(pl)))
    stop("The audit requires finite increasing distance breaks from zero and positive integer removal durations.")
  beta <- unmarked::coef(model)
  if (!is.numeric(beta) || !length(beta) || any(!is.finite(beta)) ||
      is.null(names(beta)) || anyNA(names(beta)) || anyDuplicated(names(beta)))
    stop("Stored coefficients must be finite and uniquely named.")
  blocks <- c("lambda", if (model@mixture == "NB") "alpha",
    if (model@keyfun != "uniform") "dist", "rem")
  indices <- stats::setNames(lapply(blocks, function(component) {
    selected <- names(unmarked::coef(model, type = component))
    index <- match(selected, names(beta))
    if (!length(index) || anyNA(index) || anyDuplicated(index))
      stop("Native coefficient component identities are incomplete.")
    index
  }), blocks)
  if (!identical(unname(unlist(indices)), seq_along(beta)))
    stop("Native coefficient order does not match the complete C constructor blocks.")
  n <- nrow(yd)
  active <- setdiff(blocks, "alpha")
  for (component in active) {
    x <- gd[[paste0("X_", component)]]
    rows <- if (component == "rem") n * ncol(yr) else n
    if (!is.matrix(x) || nrow(x) != rows || ncol(x) != length(indices[[component]]) || any(!is.finite(x)))
      stop("Native component design dimensions or rows are not qualified: ", component, ".")
  }
  if (model@mixture == "NB" && length(indices$alpha) != 1L)
    stop("The native NB dispersion coordinate must be scalar.")
  if (!is.null(cutoffs) && (!is.numeric(cutoffs) || !is.null(dim(cutoffs))))
    stop("Additional cutoffs must be a numeric vector of positive integers.")
  all_cutoffs <- c(settings$K, cutoffs)
  if (!is.numeric(all_cutoffs) || !is.null(dim(cutoffs)) || !length(all_cutoffs) ||
      any(!is.finite(all_cutoffs)) || any(all_cutoffs != floor(all_cutoffs)) ||
      any(all_cutoffs < max(1, rowSums(yd))) || any(all_cutoffs > 100000) || anyDuplicated(cutoffs))
    stop("Supply unique positive integer cutoffs covering the observed totals, each at most100000.")
  all_cutoffs <- unique(all_cutoffs)
  if (length(all_cutoffs) > 32L || as.double(n) * (max(all_cutoffs) + 1) > 1e7)
    stop("The requested cutoff audit exceeds its declared diagnostic work bound.")
  dll <- getLoadedDLLs()[["unmarked"]][["path"]]
  if (!is.character(dll) || length(dll) != 1L || !file.exists(dll))
    stop("The actual native C provider DLL identity is unavailable.")
  dll_hash <- digest::digest(file = dll, algo = "sha256")
  capture <- function(expr) {
    warnings <- character()
    error <- NULL
    value <- withCallingHandlers(tryCatch(force(expr), error = function(e) {
      error <<- conditionMessage(e)
      NULL
    }), warning = function(w) warnings <<- c(warnings, conditionMessage(w)))
    list(value = value, warnings = warnings, error = error)
  }
  linear <- function(component) as.numeric(gd[[paste0("X_", component)]] %*% beta[indices[[component]]])
  abundance <- exp(linear("lambda")) * binding("A")
  scale <- if (model@keyfun == "uniform") rep(Inf, n) else exp(linear("dist"))
  removal <- matrix(stats::plogis(linear("rem")), n, ncol(yr), byrow = TRUE)
  size <- if (model@mixture == "P") Inf else exp(beta[indices$alpha])
  stable <- capture(.msr_gdr_exact_observed(abundance, scale, removal, yd, yr,
    frame@dist.breaks, frame@period.lengths, model@keyfun, size))
  native <- capture({
    pd <- .msr_gdr_audit_native_bins(scale, frame@dist.breaks, model@keyfun, binding("a"), binding("u"))
    pr <- .msr_gdr_exact_removal_bins(removal, frame@period.lengths)
    .msr_gdr_exact_score_bins(abundance, pd, pr, yd, yr, size,
      probability_method = paste0("unmarked1.5.2_native_", model@keyfun,
        if (model@keyfun == "exp") "_100_trapezoids" else "_arithmetic"))
  })
  compact <- function(reference) list(nll = if (is.null(reference$value)) NA_real_ else -reference$value$logLik,
    method = if (is.null(reference$value)) NA_character_ else reference$value$probability_method,
    warnings = reference$warnings, error = reference$error)
  references <- list(native_bins = compact(native), stable_analytic = compact(stable))
  evaluations <- lapply(all_cutoffs, function(k) capture({
    # A child binding shadows K/threads without changing the original closure.
    audit_env <- new.env(parent = env)
    audit_env$K <- as.integer(k)
    audit_env$threads <- 1L
    fn <- model@nllFun
    environment(fn) <- audit_env
    value <- fn(unname(beta))
    if (length(value) != 1L || !is.finite(value)) stop("Native cutoff evaluation is nonfinite.")
    as.numeric(value)
  }))
  table <- data.frame(cutoff = all_cutoffs, original_cutoff = all_cutoffs == settings$K,
    nll = vapply(evaluations, function(x) if (is.null(x$value)) NA_real_ else x$value, numeric(1)),
    warning_count = vapply(evaluations, function(x) length(x$warnings), integer(1)))
  table$truncation_loss <- table$nll - references$native_bins$nll
  table$truncation_qualified <- !length(native$warnings) & is.finite(table$truncation_loss) &
    table$truncation_loss >= -1e-11 & table$truncation_loss <= 1e-9 & table$warning_count == 0L
  finite_poisson <- list(status = "not_qualified_for_NB", comparison = NULL,
    evaluations = NULL, fitting_admission_established = FALSE)
  if (model@mixture == "P") {
    # Keep one invariant point and only scalar records across the ladder.
    detection <- if (is.null(native$value)) NULL else
      rowSums(native$value$distance_probability) * rowSums(native$value$removal_probability)
    observed <- rowSums(yd)
    records <- lapply(all_cutoffs, function(k) compact(capture({
      if (is.null(native$value)) stop("Native-bin reference failed: ", native$error)
      correction <- .msr_gdr_poisson_cutoff_correction(abundance, detection, observed, k)
      list(logLik = sum(native$value$per_event$logLik + correction),
        probability_method = paste0(native$value$probability_method, "_finite_K_P_CDF"))
    })))
    comparison <- data.frame(cutoff = all_cutoffs,
      nll = vapply(records, function(x) x$nll, numeric(1)),
      warning_count = vapply(records, function(x) length(x$warnings), integer(1)))
    comparison$native_minus_identity_nll <- table$nll - comparison$nll
    comparison$agrees <- is.finite(comparison$nll) & is.finite(table$nll) &
      comparison$warning_count == 0L & table$warning_count == 0L &
      !length(native$warnings) & abs(comparison$native_minus_identity_nll) <= 1e-9
    finite_poisson <- list(status = "evaluated_closed_form_reference", comparison = comparison,
      evaluations = records, fitting_admission_established = FALSE)
  }
  stored_nll <- -as.numeric(unmarked::logLik(model))
  difference <- references$native_bins$nll - references$stable_analytic$nll
  result <- list(cutoffs = table, evaluations = evaluations, references = references,
    finite_poisson = finite_poisson,
    distance = list(native_minus_analytic_nll = difference,
      analytic_equivalence_qualified = is.finite(difference) && abs(difference) <= 1e-9,
      native_probability_error = native$error, stable_probability_error = stable$error),
    stored_point = list(nll = stored_nll, reevaluated_nll = table$nll[1L],
      agrees = is.finite(table$nll[1L]) && is.finite(stored_nll) && abs(table$nll[1L] - stored_nll) <= 1e-8),
    coefficients = beta, actual_cutoff = settings$K, mixture = model@mixture, keyfun = model@keyfun,
    units = list(input = frame@unitsIn, output = model@output, density = model@unitsOut),
    fingerprint = .msr_gdr_hash(list(data = frame, design = settings$design, beta = beta,
      actual_cutoff = settings$K, formulas = lapply(model@formlist, deparse))),
    provider = list(unmarked_version = as.character(utils::packageVersion("unmarked")),
      native_DLL_sha256 = dll_hash,
      native_likelihood_sha256 = .msr_gdr_hash(list(formals = deparse(formals(model@nllFun)),
        body = deparse(body(model@nllFun)))),
      audit_adapter_sha256 = .msr_gdr_hash(lapply(list(audit_gdr_integration,
        .msr_gdr_audit_native_bins, .msr_gdr_exact_observed, .msr_gdr_exact_score_bins,
        .msr_gdr_exact_distance_bins, .msr_gdr_exact_removal_bins,
        .msr_gdr_finite_poisson_score_bins, .msr_gdr_poisson_cutoff_correction,
        .msr_gdr_design_audit, .msr_gdr_design_summary, .msr_cv_prediction_terms,
        .msr_fixed_function_name, .msr_standard_prediction_function,
        .msr_prediction_pointwise_functions), function(fn)
          list(formals = deparse(formals(fn)), body = deparse(body(fn)))))),
    parameter_count = .msr_parameter_count(model), uncertainty = "none_read_only_numerical_audit",
    fitting_admission_established = FALSE, density_identification_established = FALSE)
  if (isTRUE(design)) result$design <- .msr_gdr_design_audit(gd, active,
    model@formlist$lambda, frame@siteCovs, design_data)
  result
}
