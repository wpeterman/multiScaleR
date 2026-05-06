#' Plot Marginal Effects from a Fitted Model
#'
#' Generates marginal effect plots with 95% confidence intervals for each covariate
#' in a fitted model stored within a `multiScaleR` object.
#'
#' @param x A `multiScaleR` object containing at least the elements `opt_mod` (the fitted model)
#'   and `scl_params` (a list with `mean` and `sd` for each covariate used for scaling).
#' @param ylab Character. Y-axis label for the marginal effect plots. Default is `"Estimated response"`.
#' @param length.out Integer. Number of points at which to evaluate the marginal effect curve.
#'   Default is 100.
#' @param type For `unmarked` models, Default is `"state"`
#' @param link Logical. An optional switch to predict values on the response scale. Default = `FALSE`. If predicted values seem incorrect, try switching to `TRUE`
#' @return A named list of `ggplot` objects, one for each covariate, showing the predicted
#'   response and 95% confidence interval across the observed range of that covariate
#'   while holding other covariates at their mean values.
#'
#' @details For `unmarked` models, predictions are made using `type = "state"` and
#'   the `predict` method for state variables. For other models (e.g., `lm`, `glm`),
#'   predictions are made using the standard `predict(..., se.fit = TRUE)` call and
#'   transformed by the model's inverse link function. When `x$opt_mod` is a
#'   wrapper object around another fitted model, `plot_marginal_effects()`
#'   uses the nested analysis model for predictor discovery, data recovery,
#'   and prediction when possible.
#'
#' @export
#' @importFrom unmarked coef predict
#' @importFrom insight get_parameters link_inverse find_predictors get_predicted
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line xlab ylab labs
#' @importFrom cowplot theme_cowplot
#' @importFrom utils globalVariables

plot_marginal_effects <- function(x,
                                  ylab = "Estimated response",
                                  length.out = 100,
                                  type = "state",
                                  link = FALSE) {
  if(!inherits(x, "multiScaleR")){
    stop("`x` must be a fitted `multiScaleR` object.")
  }
  validate_character_scalar(ylab, "ylab")
  validate_scalar_numeric(length.out, "length.out", integerish = TRUE, positive = TRUE)
  validate_scalar_logical(link, "link")

  if(is.null(x$opt_mod)){
    stop("`x` must contain an `opt_mod` fitted model.")
  }
  if(!is.list(x$scl_params) ||
     !all(c("mean", "sd") %in% names(x$scl_params))){
    stop("`x` must contain `scl_params` with `mean` and `sd` elements.")
  }

  mod <- x$opt_mod
  analysis_mod <- .analysis_model(mod)
  scl <- x$scl_params
  namespace(analysis_mod)

  if(inherits(analysis_mod, "glm")){
    link <- TRUE
  }

  # These are populated for non-unmarked models; used when building plot subtitles
  poly_terms  <- character(0)
  poly_vars   <- character(0)
  inter_terms <- character(0)
  inter_vars  <- character(0)

  # Extract variables from fitted model
  if (any(grepl("^unmarked", class(mod)))) {
    vars <- names(coef(mod, altNames = FALSE))[-1]
    # Find the index of the second "Int"
    second_int_index <- which(vars == "Int")[1]

    # Subset the vector up to (but not including) the second "Int"
    vars <- vars[1:(second_int_index - 1)]
    # Strip I() terms: unmarked stores them as coefficient names; prediction handles
    # them automatically from the base variable column in newdata
    c_vars <- vars <- vars[!grepl("I\\(", vars)]
    dat <- mod@data@siteCovs
    dat_ns <- dat[which(names(scl$mean) %in% vars)]

    if(is.null(type) || !type %in% c("lambda", "state")){
      stop("`type` must be specified as either 'lambda' or 'state'.")
    }

    plot_list <- lapply(vars, function(v) {
      if(!v %in% names(scl$mean)){
        dat_ns <- dat[!vars %in% names(scl$mean)]
        x_unscaled <- x_seq <- seq(min(dat_ns[v]), max(dat_ns[v]), length.out = length.out)
      } else {
        x_seq <- seq(min(mod@data@siteCovs[[v]]), max(mod@data@siteCovs[[v]]), length.out = length.out)
        x_unscaled <- (x_seq * scl$sd[v]) + scl$mean[v]
      }

      # Construct newdata with v varying, others at zero (all unmarked covariates
      # are kernel-scaled with mean = 0, so this is equivalent to holding at mean)
      newdata <- as.data.frame(matrix(rep(0, each = length.out * length(vars)), ncol = length(vars)))
      names(newdata) <- vars
      newdata[[v]] <- x_seq

      preds <- tryCatch(
        predict(mod, newdata, type = type),
        error = function(e) {
          stop(
            sprintf("Failed to compute marginal effects for covariate '%s': %s", v, e$message),
            call. = FALSE
          )
        }
      )

      fit <- preds$Predicted
      se <- preds$SE
      lwr <- preds$lower
      upr <- preds$upper

      data.frame(x = x_unscaled, fit = fit, lwr = lwr, upr = upr, variable = v)
    })

  } else {
    c_vars <- find_predictors(analysis_mod)[[1]]
    vars <- unlist(find_predictors(analysis_mod))
    dat_all <- .model_data(analysis_mod)
    if(is.null(dat_all)){
      dat_all <- extract_model_data(analysis_mod)
    }
    if(is.null(dat_all)){
      stop("Could not recover the original model data needed to plot marginal effects.")
    }

    # Parse formula term labels to detect I() polynomial and interaction terms.
    # predict() handles I(x^2) automatically when newdata contains the base column,
    # so no special prediction logic is needed. Detection is used only for plot
    # annotations and the interaction warning below.
    fterms <- tryCatch(
      attr(stats::terms(formula(analysis_mod)), "term.labels"),
      error = function(e) character(0)
    )
    poly_terms  <- fterms[grepl("^I\\(", fterms)]
    inter_terms <- fterms[grepl(":", fterms, fixed = TRUE)]

    # Which base variables appear inside an I() term?
    poly_vars <- c_vars[sapply(c_vars, function(v) {
      length(poly_terms) > 0 &&
        any(grepl(paste0("\\b", v, "\\b"), poly_terms))
    })]

    # Which base variables appear in at least one interaction term?
    inter_vars <- if(length(inter_terms) > 0) {
      c_vars[sapply(c_vars, function(v) {
        any(grepl(paste0("\\b", v, "\\b"), inter_terms))
      })]
    } else character(0)

    if(length(inter_terms) > 0) {
      message(
        "Interaction term(s) detected: ", paste(inter_terms, collapse = ", "), ".\n",
        "Each panel shows the effect of that variable holding all others at their ",
        "sample mean. This represents one conditional slice of the interaction ",
        "surface and may not capture the full effect of each variable."
      )
    }

    # Hold non-focal predictors at representative observed values. Numeric
    # columns use their sample mean; non-numeric columns keep the first
    # non-missing observed value so special terms like strata() can still be
    # evaluated during prediction.
    pred_defaults <- stats::setNames(vector("list", length(vars)), vars)
    for (var in vars) {
      if (!var %in% names(dat_all)) {
        pred_defaults[[var]] <- 0
        next
      }

      values <- dat_all[[var]]
      if (is.numeric(values)) {
        pred_defaults[[var]] <- mean(values, na.rm = TRUE)
      } else {
        keep <- which(!is.na(values))[1]
        pred_defaults[[var]] <- if (is.na(keep)) values[1] else values[keep]
      }
    }

    plot_list <- lapply(c_vars, function(v) {
      if(!v %in% names(scl$mean)){
        dat <- dat_all[, -1, drop = FALSE]
        dat_ns <- dat[!vars %in% names(scl$mean)]
        x_unscaled <- x_seq <- seq(min(dat_ns[v]), max(dat_ns[v]), length.out = length.out)
      } else {
        rng_min <- suppressWarnings(apply(dat_all, 2, min)[-1])
        rng_max <- suppressWarnings(apply(dat_all, 2, max)[-1])
        x_seq <- seq(rng_min[v], rng_max[v], length.out = length.out)
        x_unscaled <- (x_seq * scl$sd[v]) + scl$mean[v]
      }

      # Construct newdata: focal variable sweeps its range; all others at sample mean.
      # R's predict() evaluates I(x^2), log(x), etc. from the base column, so
      # newdata only needs to contain the base variable names.
      newdata <- as.data.frame(
        lapply(pred_defaults, function(value) rep(value, length.out)),
        stringsAsFactors = FALSE
      )
      newdata[[v]] <- x_seq

      if(inherits(analysis_mod, "HLfit")){
        preds_ <- tryCatch(
          predict(analysis_mod, newdata = newdata, variances = list(respVar = TRUE), re.form = NA),
          error = function(e) {
            stop(
              sprintf("Failed to compute marginal effects for covariate '%s': %s", v, e$message),
              call. = FALSE
            )
          }
        )
        preds <- list(preds = as.vector(preds_),
                      se = sqrt(attr(preds_, "fixefVar")))
      } else if(inherits(analysis_mod, "zeroinfl")) {
        preds <- tryCatch(
          as.data.frame(get_predicted(analysis_mod, data = newdata)),
          error = function(e) {
            stop(
              sprintf("Failed to compute marginal effects for covariate '%s': %s", v, e$message),
              call. = FALSE
            )
          }
        )
      } else {
        preds <- tryCatch(
          predict(analysis_mod, newdata = newdata, se.fit = TRUE),
          error = function(e) {
            stop(
              sprintf("Failed to compute marginal effects for covariate '%s': %s", v, e$message),
              call. = FALSE
            )
          }
        )
      }

      if(inherits(analysis_mod, "zeroinfl")) {
        fit_ <- preds$Predicted
        lwr  <- preds$CI_low
        upr  <- preds$CI_high
        if(is.null(lwr) && is.null(upr)) lwr <- upr <- NA

      } else if(!is.null(link_inverse(analysis_mod)) && link){
        if(!inherits(preds, 'list') || (is.list(preds) && length(preds) == 1)){
          fit_ <- link_inverse(analysis_mod)(as.data.frame(preds)[, 1])
          lwr  <- upr <- NA
        } else {
          fit_ <- link_inverse(analysis_mod)(preds[[1]])
          lwr  <- link_inverse(analysis_mod)(preds[[1]] + qnorm(0.025) * preds[[2]])
          upr  <- link_inverse(analysis_mod)(preds[[1]] + qnorm(0.975) * preds[[2]])
        }

      } else {
        if(!inherits(preds, 'list') || (is.list(preds) && length(preds) == 1)){
          fit_ <- as.data.frame(preds)[[1]]
          lwr  <- upr <- NA
        } else {
          fit_ <- preds[[1]]
          lwr  <- preds[[1]] + qnorm(0.025) * preds[[2]]
          upr  <- preds[[1]] + qnorm(0.975) * preds[[2]]
        }
      }

      data.frame(x = x_unscaled, fit = fit_, lwr = lwr, upr = upr, variable = v)
    })
  }

  names(plot_list) <- c_vars

  # Build ggplots, adding subtitles for polynomial and interaction terms
  gg_list <- lapply(c_vars, function(v) {
    df <- plot_list[[v]]

    subtitle_parts <- character(0)

    # Note any I() terms involving this variable
    if(v %in% poly_vars) {
      rel_poly <- poly_terms[sapply(poly_terms, function(pt) {
        grepl(paste0("\\b", v, "\\b"), pt)
      })]
      subtitle_parts <- c(subtitle_parts,
                          paste0("Includes: ", paste(rel_poly, collapse = ", ")))
    }

    # Note the mean-conditioning for interaction terms
    if(v %in% inter_vars) {
      partner_terms <- inter_terms[sapply(inter_terms, function(t) {
        v %in% strsplit(t, ":")[[1]]
      })]
      others <- unique(unlist(lapply(
        strsplit(partner_terms, ":"),
        function(parts) parts[parts != v]
      )))
      subtitle_parts <- c(subtitle_parts,
                          paste0("At mean of: ", paste(others, collapse = ", ")))
    }

    subtitle <- if(length(subtitle_parts) > 0) paste(subtitle_parts, collapse = "; ") else NULL

    p <- ggplot(df, aes(x = x, y = fit)) +
      geom_line(linewidth = 1) +
      xlab(df$variable[1]) +
      ylab(ylab) +
      theme_cowplot()

    if(!is.null(subtitle))
      p <- p + labs(subtitle = subtitle)

    if(!all(is.na(df$lwr)) && !all(is.na(df$upr)))
      p <- p + geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25)

    p
  })

  names(gg_list) <- c_vars

  lapply(gg_list, print)
  invisible(gg_list)
}


safe_predict <- function(mod, newdata) {
  pkg <- extract_namespace(mod)

  if(is.null(pkg)){
    mod_classes <- class(mod)

    # Try S3 generic method (e.g., predict.lm)
    pkg <- NULL
    for (cls in mod_classes) {
      mthd <- tryCatch(
        getS3method("predict", cls, optional = TRUE),
        error = function(e) NULL
      )
      if (!is.null(mthd)) {
        env <- environment(mthd)
        if (!is.null(env) && !is.null(pkgname <- environmentName(env)) &&
            pkgname %in% loadedNamespaces()) {
          pkg <- pkgname
          break
        }
      }
    }

    # Try S4 class definition
    if (is.null(pkg) && isS4(mod)) {
      for (cls in mod_classes) {
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

  # Load namespace if needed
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
  if (!paste0("package:", pkg) %in% search()) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }

  # Dispatch predict
  predict(mod, newdata, se.fit = TRUE)
}

namespace <- function(x) {

  pkg <- extract_namespace(x)

  if (is.null(pkg)) {
    mod_classes <- class(x)

    # Try S3 generic method (e.g., predict.lm)
    pkg <- NULL
    for (cls in mod_classes) {
      mthd <- tryCatch(
        getS3method("predict", cls, optional = TRUE),
        error = function(e) NULL
      )
      if (!is.null(mthd)) {
        env <- environment(mthd)
        if (!is.null(env) && !is.null(pkgname <- environmentName(env)) &&
            pkgname %in% loadedNamespaces()) {
          pkg <- pkgname
          break
        }
      }
    }

    # Try S4 class definition
    if (is.null(pkg) && isS4(x)) {
      for (cls in mod_classes) {
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

  # Check if namespace is already loaded
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
  if (!paste0("package:", pkg) %in% search()) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }

  invisible(pkg)
}

utils::globalVariables(c("fit", "lwr", "upr"))
