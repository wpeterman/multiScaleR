# Methods for multiScaleR_joint objects (experimental).

#' Methods for joint scale fits
#'
#' @description
#' Print, summarize, and extract results from a joint fit made with
#' [multiScale_joint()]. `logLik()` makes [stats::AIC()] and [stats::BIC()]
#' work, so several structures fitted to the same data sets can be compared
#' directly. `anova()` gives likelihood-ratio tests between nested structures.
#' A test of a `"common"` fit against a `"footprint_scaled"` fit tests
#' `kappa = 0`, on the boundary of the parameter space, and uses the matching
#' mixture of chi-square distributions (half the chi-square p-value for one
#' `kappa`).
#'
#' @param x,object A `multiScaleR_joint` object.
#' @param digits Number of significant digits to print.
#' @param ... Further `multiScaleR_joint` objects for `anova()`; ignored by the
#'   other methods.
#'
#' @return `print()` returns `x` invisibly. `summary()` returns an object of
#'   class `summary_multiScaleR_joint` holding the parameter, scale,
#'   comparison, and per-data-set coefficient tables. `logLik()` returns a
#'   `logLik` object whose `df` counts every native parameter plus the joint
#'   scale parameters, `kappa`, and shared coefficients. `coef()` returns the
#'   joint parameters on their natural scales (scales in map units). `vcov()`
#'   returns their covariance on the optimizer scale (log scales, `kappa`,
#'   coefficients). `nobs()` returns the total number of observations. `anova()`
#'   returns a data frame of likelihood-ratio tests.
#'
#' @seealso [multiScale_joint()]
#' @name multiScaleR_joint-methods
NULL


#' @rdname multiScaleR_joint-methods
#' @export
print.multiScaleR_joint <- function(x, digits = 4, ...) {
  cat("\nJoint scale-of-effect fit (experimental)\n")
  sets <- names(x$fits)
  n <- vapply(x$fits, function(f) .msr_model_nobs(f$opt_mod), numeric(1))
  cls <- vapply(x$fits, function(f) class(.analysis_model(f$opt_mod))[1L], character(1))
  cat("Data sets: ", paste(sprintf("%s (%s, n = %d)", sets, cls, as.integer(n)),
                           collapse = "; "), "\n", sep = "")
  cat("Kernel: ", x$kernel, "\n", sep = "")

  cat("\nScale structure\n")
  print(x$structure, row.names = FALSE)
  if (!is.null(x$footprint)) {
    cat("Footprint SD (map units): ",
        paste(sprintf("%s = %s", names(x$footprint), format(x$footprint, digits = digits)),
              collapse = ", "), "\n", sep = "")
  }

  level <- format(100 * x$level)
  cat("\nJoint parameters (scales in map units)\n")
  print(.msr_joint_print_table(x$parameters, c("parameter", "estimate", "SE"), level,
                               digits), row.names = FALSE)

  cat("\nApparent scale in each data set\n")
  print(.msr_joint_print_table(x$scales, c("data_set", "covariate", "sigma", "SE"), level,
                               digits), row.names = FALSE)

  aic <- -2 * x$logLik + 2 * x$df
  cat(sprintf("\nlogLik = %s (df = %d), AIC = %s, n = %d\n",
              format(x$logLik, digits = digits + 2), as.integer(x$df),
              format(aic, digits = digits + 2), as.integer(x$nobs)))
  .msr_joint_print_comparison(x$comparison, digits)

  if (!is.null(x$residual_correlation)) {
    cat("\nResidual correlation at shared clusters (Spearman)\n")
    print(x$residual_correlation, row.names = FALSE, digits = 3)
  }
  notes <- .msr_joint_notes(x)
  if (length(notes)) cat("\n", paste0("Note: ", notes, collapse = "\n"), "\n", sep = "")
  invisible(x)
}


.msr_joint_print_table <- function(tab, keep, level, digits) {
  out <- tab[, intersect(keep, names(tab)), drop = FALSE]
  use_profile <- "profile_lower" %in% names(tab) && any(is.finite(tab$profile_lower))
  lo <- if (use_profile) tab$profile_lower else tab$wald_lower
  hi <- if (use_profile) tab$profile_upper else tab$wald_upper
  label <- if (use_profile) "profile" else "Wald"
  out[[sprintf("%s%% lower (%s)", level, label)]] <- lo
  out[[sprintf("%s%% upper (%s)", level, label)]] <- hi
  num <- vapply(out, is.numeric, logical(1))
  out[num] <- lapply(out[num], signif, digits = digits)
  out
}


.msr_joint_print_comparison <- function(cmp, digits) {
  if (cmp$df_difference <= 0) {
    cat("Same number of parameters as separate fits; no likelihood-ratio test.\n")
    return(invisible(NULL))
  }
  # Every structure lies inside the separate-scales model, away from its
  # boundary, so this is an ordinary chi-square test.
  cat(sprintf("Against separate scales: LR chi-square = %s, df = %d, p = %s\n",
              format(cmp$statistic, digits = digits), as.integer(cmp$df_difference),
              format.pval(cmp$p_value, digits = 3)))
  invisible(NULL)
}


.msr_joint_notes <- function(x) {
  notes <- character()
  if (length(x$uncertainty$at_bound))
    notes <- c(notes, paste0("At a search bound (SE held fixed): ",
                             paste(x$uncertainty$at_bound, collapse = ", "), "."))
  if ("bound_limited" %in% names(x$parameters) && any(x$parameters$bound_limited))
    notes <- c(notes, paste0("Profile interval reaches a search bound for: ",
                             paste(x$parameters$parameter[x$parameters$bound_limited],
                                   collapse = ", "), "."))
  if (identical(x$uncertainty$method, "hessian"))
    notes <- c(notes, "Standard errors assume the data sets are independent given the covariates.")
  if (identical(x$uncertainty$method, "jackknife"))
    notes <- c(notes, sprintf("Wald SEs from a delete-a-group jackknife (%d replicates); profile limits assume independence.",
                              as.integer(x$uncertainty$jackknife$replicates)))
  if (length(x$share))
    notes <- c(notes, "Shared coefficients use one standardization pooled across data sets (see $standardization).")
  notes
}


#' @rdname multiScaleR_joint-methods
#' @export
summary.multiScaleR_joint <- function(object, ...) {
  coefficients <- lapply(object$fits, function(f) {
    model <- .analysis_model(f$opt_mod)
    tab <- tryCatch(as.data.frame(insight::get_parameters(model)), error = function(e) NULL)
    tab
  })
  structure(list(call = object$call, structure = object$structure,
                 parameters = object$parameters, scales = object$scales,
                 comparison = object$comparison, coefficients = coefficients,
                 residual_correlation = object$residual_correlation,
                 logLik = object$logLik, df = object$df, nobs = object$nobs,
                 level = object$level, footprint = object$footprint,
                 standardization = object$standardization, joint = object),
            class = "summary_multiScaleR_joint")
}


#' @rdname multiScaleR_joint-methods
#' @export
print.summary_multiScaleR_joint <- function(x, digits = 4, ...) {
  cat("\nCall:\n")
  print(x$call)
  print(x$joint, digits = digits)
  cat("\nNative coefficients in each data set at the joint scales\n")
  for (d in names(x$coefficients)) {
    cat("\n", d, "\n", sep = "")
    tab <- x$coefficients[[d]]
    if (is.null(tab)) cat("  (not available for this model class)\n") else
      print(tab, row.names = FALSE, digits = digits)
  }
  if (!is.null(x$standardization)) {
    cat("\nPooled standardization of shared covariates\n")
    print(data.frame(covariate = names(x$standardization$center),
                     center = x$standardization$center,
                     spread = x$standardization$spread), row.names = FALSE,
          digits = digits)
  }
  invisible(x)
}


#' @rdname multiScaleR_joint-methods
#' @exportS3Method stats::logLik
logLik.multiScaleR_joint <- function(object, ...) {
  structure(object$logLik, df = object$df, nobs = object$nobs, class = "logLik")
}


#' @rdname multiScaleR_joint-methods
#' @exportS3Method stats::nobs
nobs.multiScaleR_joint <- function(object, ...) object$nobs


#' @rdname multiScaleR_joint-methods
#' @exportS3Method stats::coef
coef.multiScaleR_joint <- function(object, ...) {
  stats::setNames(object$parameters$estimate, object$parameters$parameter)
}


#' @rdname multiScaleR_joint-methods
#' @exportS3Method stats::vcov
vcov.multiScaleR_joint <- function(object, ...) object$vcov


#' @rdname multiScaleR_joint-methods
#' @exportS3Method stats::anova
anova.multiScaleR_joint <- function(object, ...) {
  fits <- c(list(object), list(...))
  if (!all(vapply(fits, inherits, logical(1), "multiScaleR_joint")))
    stop("anova() compares multiScaleR_joint fits only.", call. = FALSE)
  if (length(fits) < 2L) stop("Supply two or more joint fits to compare.", call. = FALSE)
  sets <- lapply(fits, function(f) sort(names(f$fits)))
  n <- vapply(fits, `[[`, numeric(1), "nobs")
  if (length(unique(sets)) != 1L || length(unique(n)) != 1L)
    stop("Joint fits must combine the same data sets with the same observations.",
         call. = FALSE)
  label <- vapply(fits, function(f) paste(sprintf("%s=%s", f$structure$covariate,
                                                  f$structure$structure),
                                          collapse = "; "), character(1))
  shared <- vapply(fits, function(f) if (length(f$share))
    paste0(" share=", paste(f$share, collapse = ",")) else "", character(1))
  tab <- data.frame(df = vapply(fits, `[[`, numeric(1), "df"),
                    logLik = vapply(fits, `[[`, numeric(1), "logLik"),
                    row.names = make.unique(paste0(label, shared)))
  tab$AIC <- -2 * tab$logLik + 2 * tab$df
  ord <- order(tab$df)
  tab <- tab[ord, , drop = FALSE]
  fits <- fits[ord]
  tab$LR <- c(NA, 2 * diff(tab$logLik))
  tab$df_difference <- c(NA, diff(tab$df))
  tab$p_value <- ifelse(is.finite(tab$df_difference) & tab$df_difference > 0,
                        stats::pchisq(pmax(tab$LR, 0), tab$df_difference, lower.tail = FALSE),
                        NA_real_)
  # A common scale is "footprint_scaled" with kappa = 0, the edge of kappa >= 0.
  # With q such kappas the statistic follows a mixture of chi-square
  # distributions with 0 to q degrees of freedom and binomial(q, 1/2) weights,
  # which halves the p-value when q = 1.
  boundary <- FALSE
  for (i in seq_len(nrow(tab))[-1L]) {
    q <- .msr_joint_kappa_boundary(fits[[i - 1L]], fits[[i]])
    if (q > 0L && isTRUE(tab$df_difference[i] == q)) {
      lr <- max(tab$LR[i], 0)
      tab$p_value[i] <- if (lr <= 0) 1 else
        sum(choose(q, seq_len(q)) * 0.5^q * stats::pchisq(lr, seq_len(q), lower.tail = FALSE))
      boundary <- TRUE
    }
  }
  attr(tab, "heading") <- paste0(
    "Likelihood-ratio tests between joint fits (valid only for nested structures)\n",
    if (boundary) paste0("A test of a common scale against footprint_scaled tests kappa = 0, ",
                         "on the boundary of the parameter space;\nits p-value uses the ",
                         "matching chi-square mixture (half the chi-square p-value for one kappa).\n"))
  class(tab) <- c("anova", "data.frame")
  tab
}


# Number of covariates that are "footprint_scaled" in the larger fit and
# "common" in the smaller one, with every other structure and shared
# coefficient the same. Their kappas are 0 under the smaller fit.
.msr_joint_kappa_boundary <- function(small, big) {
  s <- stats::setNames(small$structure$structure, small$structure$covariate)
  b <- stats::setNames(big$structure$structure, big$structure$covariate)
  if (!setequal(names(s), names(b)) || !setequal(small$share, big$share)) return(0L)
  b <- b[names(s)]
  hit <- s == "common" & b == "footprint_scaled"
  if (!any(hit) || any(s[!hit] != b[!hit])) return(0L)
  sum(hit)
}
