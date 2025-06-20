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
#'   transformed by the model's inverse link function .
#'
#' @export
#' @importFrom unmarked coef predict
#' @importFrom insight get_parameters link_inverse find_predictors
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line xlab ylab
#' @importFrom cowplot theme_cowplot

plot_marginal_effects <- function(x,
                                  ylab = "Estimated response",
                                  length.out = 100,
                                  type = "state",
                                  link = FALSE) {
  mod <- x$opt_mod
  scl <- x$scl_params
  namespace(mod)

  if(inherits(mod, "glm")){
    link <- TRUE
  }

  # Extract variables from fitted model
  if (any(grepl("^unmarked", class(mod)))) {
    vars <- names(coef(mod, altNames = FALSE, type = 'state'))[-1]
    dat <- mod@data@siteCovs
    dat_ns <- dat[!vars %in% names(scl$mean)]

    if(is.null(type)){
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

      # Construct newdata with v varying, others at mean
      newdata <- as.data.frame(matrix(rep(0, each = length.out * length(vars)), ncol = length(vars)))
      names(newdata) <- vars
      newdata[[v]] <- x_seq

      preds <- predict(mod, newdata, type = type)

      fit <- preds$Predicted
      se <- preds$SE
      lwr <- preds$lower
      upr <- preds$upper

      data.frame(
        x = x_unscaled,
        fit = fit,
        lwr = lwr,
        upr = upr,
        variable = v
      )
    })
  } else {
    # vars <- names(coef(mod))[-1]
    vars <- find_predictors(mod)[[1]]

    plot_list <- lapply(vars, function(v) {
      if(!v %in% names(scl$mean)){
        dat <- get_data(mod)[,-1]
        dat_ns <- dat[!vars %in% names(scl$mean)]
        x_unscaled <- x_seq <- seq(min(dat_ns[v]), max(dat_ns[v]), length.out = length.out)
      } else {
        min <- suppressWarnings(apply(get_data(mod), 2, min)[-1])
        max <- suppressWarnings(apply(get_data(mod), 2, max)[-1])
        x_seq <- seq(min[v], max[v], length.out = length.out)
        x_unscaled <- (x_seq * scl$sd[v]) + scl$mean[v]
      }

      # Construct newdata with v varying, others at mean
      newdata <- as.data.frame(matrix(rep(0, each = length.out * length(vars)), ncol = length(vars)))
      names(newdata) <- vars
      newdata[[v]] <- x_seq


      # preds <- safe_predict(mod, newdata = newdata)
      if(inherits(mod, "HLfit")){
        preds_ <- predict(mod, newdata = newdata, variances = list(respVar = TRUE), re.form = NA)
        preds <- list(preds = as.vector(preds_),
                      se = sqrt(attr(preds_, "fixefVar")))
      } else {
        preds <- predict(mod, newdata = newdata, se.fit = T)
      }

      if(!is.null(link_inverse(mod)) && link){
        if(!inherits(preds,'list') || (is.list(preds) && length(preds) == 1)){
          fit_ <- link_inverse(mod)(as.data.frame(preds)[,1])
          lwr <- upr <- NA
        } else {
          fit_ <- link_inverse(mod)(preds[[1]])
          lwr <- link_inverse(mod)(preds[[1]] + qnorm(0.025)*preds[[2]])
          upr <- link_inverse(mod)(preds[[1]] + qnorm(0.975)*preds[[2]])
        }

      } else {
        if(!inherits(preds,'list') || (is.list(preds) && length(preds) == 1)){
          fit_ <- as.data.frame(preds)[[1]]
          lwr <- upr <- NA
        } else {
          fit_ <- (preds[[1]])
          lwr <- (preds[[1]] + qnorm(0.025)*preds[[2]])
          upr <- (preds[[1]] + qnorm(0.975)*preds[[2]])
        }
      }


      data.frame(
        x = x_unscaled,
        fit = fit_,
        lwr = lwr,
        upr = upr,
        variable = v
      )
    })
  }

  names(plot_list) <- vars

  # Build ggplots
  gg_list <- lapply(plot_list, function(df) {
    p <- ggplot(df, aes(x = x, y = fit)) +
      geom_line(linewidth = 1) +
      xlab(df$variable[1]) +
      ylab(ylab) +
      theme_cowplot()

    if (!all(is.na(df$lwr)) && !all(is.na(df$upr))) {
      p <- p + geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25)
    }

    p
  })

  names(gg_list) <- vars

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
  if (!isNamespaceLoaded(pkg)) {
    attachNamespace(pkg)
  }

  # Dispatch predict
  predict(mod, newdata, se.fit = TRUE)
}

namespace <- function(x) {
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
  if (!isNamespaceLoaded(pkg)) {
    # Load namespace if needed
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed.")
    }
    attachNamespace(pkg)
  }

  invisible(pkg)
}
