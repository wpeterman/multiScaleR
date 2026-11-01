#' Plot Marginal Effects from a Fitted Model
#'
#' Generates marginal effect plots with 95\% intervals for each
#' covariate in the optimized model stored within a \code{multiScaleR} object.
#' Each panel sweeps one covariate across its observed range while holding all
#' other numeric predictors at their sample mean. For categorical predictors,
#' the grid uses observed levels and other categorical inputs use an observed
#' reference level. For \code{unmarked} models,
#' other categorical predictors are held at their first observed level.
#'
#' @param x A fitted \code{multiScaleR} object returned by
#'   \code{\link{multiScale_optim}}. Must contain \code{opt_mod} (the optimized
#'   fitted model) and \code{scl_params} (a list with \code{mean} and \code{sd}
#'   elements for each covariate, used to back-transform scaled covariates to
#'   their original units for the x-axis). Also accepts \code{multiScaleR_bayes}
#'   fits for hypothetical fixed-reference covariate curves; see Details.
#' @param ylab Character scalar. Y-axis label for all marginal effect panels.
#'   Default: \code{"Estimated response"}.
#' @param length.out Positive integer. Number of equally spaced points at which
#'   to evaluate the marginal effect curve for each covariate. Default: \code{100}.
#' @param type Character. Prediction type for \code{unmarked} models. One of
#'   \code{"state"} (default) or \code{"lambda"}. Ignored for non-unmarked
#'   models.
#' @param link Logical. If \code{TRUE}, predictions are passed through the
#'   model's inverse link function before plotting (i.e., plotted on the
#'   response scale). For \code{glm} and spaMM \code{HLfit} objects this is applied automatically.
#'   Set to \code{TRUE} manually if predicted values appear to be on the link
#'   scale. Default: \code{FALSE}.
#' @param draws For Bayesian fits, a positive number of stored draws to select
#'   evenly across the retained sample, or explicit draw indices. The default
#'   \code{NULL} uses at most 1000 evenly spaced draws. Ignored for frequentist
#'   fits. Bayesian curves use the response scale regardless of \code{link}.
#'
#' @return A named list of \code{ggplot2} objects, one per covariate. Plots are
#'   printed as a side effect and the list is returned invisibly. Each plot
#'   shows:
#'   \itemize{
#'     \item The predicted response as a solid line for numeric predictors or
#'       points for categorical predictors, with x-axis values back-transformed
#'       to the original (unscaled) units when applicable.
#'     \item A shaded ribbon or categorical error bars for the 95\% confidence
#'       interval (when available from the model's \code{predict} method).
#'       Bayesian curves show posterior medians and central 95\% credible
#'       intervals, with the retained draw indices and prediction contract
#'       attached to the returned list.
#'     \item A subtitle noting any polynomial (\code{I(x^2)}) or interaction
#'       terms that involve the covariate.
#'   }
#'
#' @details
#' Marginal effects are computed by constructing a \code{newdata} data frame
#' where the focal covariate sweeps from its minimum to maximum observed value
#' and all other predictors are held at their sample mean. For kernel-scaled
#' covariates, whose scaled mean is 0, this is equivalent to holding them at
#' their average landscape value. Variables used only in model offsets are held
#' at their sample mean for prediction but do not receive marginal effect panels.
#' Quoted predictor names retain their original spelling in prediction grids.
#' Factor contrasts and fixed transformations are evaluated by the fitted model.
#'
#' For \code{glm} and most GLM-family models, predictions are transformed via
#' the inverse link function automatically (e.g., \code{exp()} for Poisson
#' log-link models). Confidence intervals are constructed from
#' \code{predict(..., se.fit = TRUE)} using \eqn{\pm 1.96 \times SE}.
#'
#' Native link-scale limits are transformed to the response scale and ordered
#' for increasing or decreasing links. Limits with nonfinite standard errors,
#' invalid native link/family endpoints, or an inverse-link pole inside the
#' interval are unavailable and produce a warning. Fitted values remain native;
#' the method does not clip or repair invalid confidence limits.
#'
#' For \code{unmarked} models, \code{predict(mod, newdata, type = type)} is
#' used, and lower/upper bounds come from the \code{lower} and \code{upper}
#' columns of the returned data frame.
#' Predictors are recovered from the selected component's formula, so factor
#' contrasts such as \code{Year2} are not mistaken for data columns. Factors
#' receive one point and confidence interval per observed level. Site keys
#' not used as fixed predictors and detection-only predictors do not receive
#' state-effect panels.
#' Numeric predictors not optimized by \code{multiScale_optim()} retain their
#' original units and are held at their sample mean in other panels.
#'
#' For \code{unmarked} random-effects models, prediction uses
#' \code{re.form = NA}, which sets random effects to zero. These curves are
#' conditional predictions at the stated covariate values, not averages over
#' years, sites, or the random-effect distribution. Confidence intervals are
#' conditional on the fitted spatial scale.
#' Explicit nondefault factor contrasts and learned fixed bases in ordinary
#' \code{pcount} (Poisson or negative-binomial) and \code{occu} fits use
#' \code{unmarked::linearComb()} and \code{unmarked::backTransform()} to preserve
#' the fitted coefficient basis and link. Learned prediction calls are trained
#' once on retained site covariates and held fixed across plotting grids.
#' Other unmarked classes requiring that fallback stop rather than return a
#' different basis or derived-state prediction. Active terms select panels:
#' removed terms and random-only groups are excluded, while admitted quoted
#' names and fixed logical expressions retain their raw inputs. A logical pipe
#' within \code{I()} is a fixed operation, not a random-effect declaration.
#' The native unmarked parser may reject pipe-containing names or fixed
#' expressions before fitting. Precompute such fixed inputs and use safe
#' column names when that occurs; this plotting adapter cannot change a
#' formula rejected by the native fitting function.
#'
#' For models with interaction terms, a message is printed explaining that each
#' panel represents a conditional slice (at the mean of the interacting
#' variable) rather than the full interaction surface.
#'
#' \code{HLfit} (spaMM) and \code{zeroinfl} (pscl) models are handled with
#' class-specific prediction calls.
#' spaMM curves set random effects to zero, evaluate the linear predictor, and
#' transform native link-scale \code{fixefVar} limits through the inverse link.
#' Its response variance is not used as a confidence interval for the mean.

#' When \code{x$opt_mod} is a wrapper object around another fitted model,
#' \code{plot_marginal_effects()} uses the nested analysis model for predictor
#' discovery, data recovery, and prediction when possible.
#'
#' Mixed-model curves set response random effects to zero. glmmTMB curves use
#' native response predictions, including zero inflation or truncated-family
#' means, with native response-scale Wald limits when available. This is a
#' prediction at zero random effects, not an average over their distribution.
#' Frequentist limits condition on the selected spatial scales. Cox panels
#' default to log relative risk; exponentiation gives a relative risk score,
#' not an occurrence or survival probability.
#'
#' Bayesian panels use the stored marginal response-coefficient draws and the
#' original pooled reference normalization. They hold each landscape covariate
#' at a hypothetical grid of values while keeping the remaining covariates at
#' their reference values. Coefficient draws retain their fitted posterior
#' variation, but the landscape is not reevaluated at each draw's kernel scales.
#' Response random effects are zero and observation noise is excluded. Use
#' \code{\link{predict.multiScaleR_bayes}} with matched spatial inputs for
#' location-specific posterior expectations. A plotted credible interval does
#' not establish convergence or calibration; existing fit diagnostics remain
#' relevant and their stored status is shown in the caption.
#'
#' @examples
#' \donttest{
#' data('pts')
#' data('count_data')
#' hab <- terra::rast(system.file('extdata', 'hab.tif', package = 'multiScaleR'))
#'
#' kernel_inputs <- kernel_prep(pts = pts,
#'                              raster_stack = hab,
#'                              max_D = 250,
#'                              kernel = 'gaussian')
#'
#' mod <- glm(y ~ hab, family = poisson, data = count_data)
#' opt <- multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs)
#'
#' ## Default marginal effect plot on the response scale
#' plot_marginal_effects(opt)
#'
#' ## Custom y-axis label
#' plot_marginal_effects(opt, ylab = "Predicted abundance")
#'
#' ## Finer evaluation grid
#' plot_marginal_effects(opt, length.out = 200)
#' }
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
                                  link = FALSE, draws = NULL) {
  validate_scalar_logical(link, "link")
  if (inherits(x, "multiScaleR_bayes"))
    return(.msr_bayes_marginal_effects(x, ylab, length.out, draws))
  if (!is.null(draws)) stop("draws applies only to Bayesian effect curves.", call. = FALSE)
  if(!inherits(x, "multiScaleR")){
    stop("`x` must be a fitted `multiScaleR` object.")
  }
  validate_character_scalar(ylab, "ylab")
  validate_scalar_numeric(length.out, "length.out", integerish = TRUE, positive = TRUE)

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
  native_mixed_response <- inherits(analysis_mod, c("glmmTMB", "merMod"))
  if (missing(ylab) && inherits(analysis_mod, "coxph"))
    ylab <- if (link) "Relative risk" else "Log relative risk"

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
    if(is.null(type) || !type %in% c("lambda", "state")){
      stop("`type` must be specified as either 'lambda' or 'state'.")
    }

    component <- .msr_unmarked_marginal_formula(mod, type)
    dat <- unmarked::siteCovs(unmarked::getData(mod))
    active_terms <- stats::terms(component, data = dat)
    vars <- .msr_marginal_variables(active_terms, offsets = TRUE)
    c_vars <- .msr_marginal_variables(active_terms, offsets = FALSE)
    # Native model.frame can still evaluate variables removed by '-'. Supply
    # their original context without adding them or random groups as panels.
    frame_variables <- as.list(attr(active_terms, "variables"))[-1L]
    vars <- unique(c(vars, unlist(lapply(frame_variables, all.vars), use.names = FALSE)))
    missing_vars <- setdiff(vars, names(dat))
    if (length(missing_vars)) {
      stop("Could not recover site covariates for marginal effects: ",
           paste(missing_vars, collapse = ", "), call. = FALSE)
    }

    # Use column names, not coefficient names or column positions. Repeated
    # site-year keys may coexist with factors, offsets, and unscaled predictors.
    pred_defaults <- lapply(dat[vars], function(values) {
      if (is.numeric(values)) return(mean(values, na.rm = TRUE))
      values[which(!is.na(values))[1]]
    })
    factor_vars <- vars[vapply(dat[vars], is.factor, logical(1))]
    custom_contrasts <- any(vapply(dat[factor_vars], function(values) {
      fitted <- attr(values, "contrasts")
      if (is.null(fitted)) return(FALSE)
      attr(values, "contrasts") <- NULL
      !isTRUE(all.equal(fitted, stats::contrasts(values)))
    }, logical(1)))
    # Freeze learned fixed prediction calls once on the native training rows.
    basis <- .msr_unmarked_curve_basis(active_terms, dat)
    plot_list <- lapply(c_vars, function(v) {
      values <- dat[[v]]
      if (all(is.na(values))) {
        stop("Covariate '", v, "' has no observed values.", call. = FALSE)
      }
      if (is.numeric(values)) {
        x_seq <- seq(min(values, na.rm = TRUE), max(values, na.rm = TRUE),
                     length.out = length.out)
      } else if (is.factor(values)) {
        # Retain factor levels, ordering, and any fitted contrast attributes.
        x_seq <- values[match(levels(droplevels(values)), as.character(values))]
      } else {
        x_seq <- unique(values[!is.na(values)])
      }
      x_unscaled <- x_seq
      if (v %in% names(scl$mean)) {
        x_unscaled <- (x_seq * scl$sd[v]) + scl$mean[v]
      }

      newdata <- as.data.frame(lapply(pred_defaults, function(value) {
        rep(value, length(x_seq))
      }), stringsAsFactors = FALSE, check.names = FALSE)
      newdata[[v]] <- x_seq

      # unmarked's data-frame prediction path drops explicit contrasts in
      # model.frame(xlev = ...). Use its public linear-combination interface
      # for supported custom designs; otherwise native prediction is unchanged.
      for (factor_var in factor_vars) {
        attr(newdata[[factor_var]], "contrasts") <-
          if (custom_contrasts) attr(dat[[factor_var]], "contrasts") else NULL
      }
      preds <- tryCatch(
        if (custom_contrasts || basis$learned) {
          .msr_unmarked_contrast_prediction(mod, basis$terms, newdata, type)
        } else {
          unmarked::predict(mod, newdata = newdata, type = type, re.form = NA)
        },
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
    if (identical(x$backend, "conditional_logit")) {
      vars <- .msr_active_formula_variables(x$profile_fit$contract$terms)
      c_vars <- vars
    } else if (inherits(analysis_mod, c("glm", "lm"))) {
      # Native terms retain active variables and quoted names without a
      # provider's optional name-sanitizing interpretation.
      vars <- .msr_active_formula_variables(stats::terms(analysis_mod))
      c_vars <- vars
    } else {
      c_vars <- find_predictors(analysis_mod)[[1]]
      vars <- unlist(find_predictors(analysis_mod))
    }
    # Offset-only variables describe effort, not covariate effects. Keep them
    # in vars so prediction holds them at their mean, but omit their panels.
    # A variable that also enters as a model term remains a plotted covariate.
    offset_vars <- tryCatch(insight::find_offset(analysis_mod), error = function(e) NULL)
    if (length(offset_vars) > 0) {
      term_vars <- tryCatch({
        labels <- attr(stats::terms(formula(analysis_mod)), "term.labels")
        if (length(labels)) all.vars(stats::reformulate(labels)) else character(0)
      }, error = function(e) character(0))
      c_vars <- setdiff(c_vars, setdiff(offset_vars, term_vars))
    }
    dat_all <- if (identical(x$backend, "joint_mixed") &&
      is.data.frame(x$mixed_fit$prediction_data)) x$mixed_fit$prediction_data else
        .model_data(analysis_mod)
    if(is.null(dat_all)){
      dat_all <- extract_model_data(analysis_mod)
    }
    if(is.null(dat_all)){
      stop("Could not recover the original model data needed to plot marginal effects.")
    }

    design <- x$profile_fit$contract$objective$response_design
    if (!is.null(design$expression_programs)) {
      # Evaluated powers cannot identify the underlying signed kernel values.
      # Rebuild those values from the selected, verified spatial preparation.
      e <- x$profile_fit$contract$objective$engine
      if (inherits(e, "msr_surface_engine") && isTRUE(e$rebuild_from_inputs))
        e <- .msr_prepare_surface_engine(x$kernel_inputs, e$covariates, e$rows,
                                         e$standardize, e$implementation)
      value <- .msr_eval_kernel_engine(e,
        x$profile_fit$par[seq_along(e$parameter_names)])$value
      if (nrow(dat_all) != nrow(value)) stop("Effect-curve training rows do not match the fitted kernels.")
      for (nm in e$covariates) dat_all[[nm]] <- value[, nm]
      missing <- setdiff(vars, names(dat_all))
      if (length(missing)) stop("Raw fixed-transform inputs were not retained for effect curves: ",
        paste(missing, collapse = ", "), ". Refit with recoverable raw training data.", call. = FALSE)
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
    term_order <- tryCatch(attr(stats::terms(formula(analysis_mod)), "order"),
      error = function(e) integer())
    inter_terms <- fterms[term_order > 1L]
    poly_variables <- lapply(poly_terms, function(term) all.vars(str2lang(term)))
    inter_variables <- lapply(inter_terms, function(term) all.vars(str2lang(term)))

    # Which base variables appear inside an I() term?
    poly_vars <- c_vars[sapply(c_vars, function(v) {
      any(vapply(poly_variables, function(names) v %in% names, logical(1)))
    })]

    # Which base variables appear in at least one interaction term?
    inter_vars <- if(length(inter_terms) > 0) {
      c_vars[sapply(c_vars, function(v) {
        any(vapply(inter_variables, function(names) v %in% names, logical(1)))
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
      values <- dat_all[[v]]
      if (all(is.na(values)))
        stop("Covariate '", v, "' has no observed values.", call. = FALSE)
      if (is.numeric(values)) {
        x_seq <- seq(min(values, na.rm = TRUE), max(values, na.rm = TRUE),
                     length.out = length.out)
      } else if (is.factor(values)) {
        # Preserve fitted levels, ordering, and contrasts in native prediction.
        x_seq <- values[match(levels(droplevels(values)), as.character(values))]
      } else {
        x_seq <- unique(values[!is.na(values)])
      }
      x_unscaled <- x_seq
      if (v %in% names(scl$mean)) {
        x_unscaled <- (x_seq * scl$sd[v]) + scl$mean[v]
      }

      # Construct newdata: focal variable sweeps its range; all others at sample mean.
      # R's predict() evaluates I(x^2), log(x), etc. from the base column, so
      # newdata only needs to contain the base variable names.
      newdata <- as.data.frame(
        lapply(pred_defaults, function(value) rep(value, length(x_seq))),
        stringsAsFactors = FALSE, check.names = FALSE
      )
      newdata[[v]] <- x_seq

      if(inherits(analysis_mod, "HLfit")){
        preds_ <- tryCatch(
          predict(analysis_mod, newdata = newdata, type = "link",
            variances = list(fixefVar = TRUE), re.form = NA),
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
      } else if (native_mixed_response) {
        # Native response prediction also handles zero-inflation and truncated
        # families. Transforming only the conditional count link would not.
        preds <- tryCatch(
          predict(analysis_mod, newdata = newdata, type = "response",
                  re.form = NA, se.fit = TRUE),
          error = function(e) stop(sprintf(
            "Failed to compute marginal effects for covariate '%s': %s", v,
            e$message), call. = FALSE)
        )
      } else if (inherits(analysis_mod, "coxph")) {
        preds <- predict(analysis_mod, newdata = newdata, type = "lp",
                         reference = "zero", se.fit = TRUE)
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

      } else if(!native_mixed_response && !is.null(link_inverse(analysis_mod)) && (link || inherits(analysis_mod, "HLfit"))){
        if(!inherits(preds, 'list') || (is.list(preds) && length(preds) == 1)){
          fit_ <- link_inverse(analysis_mod)(as.data.frame(preds)[, 1])
          lwr  <- upr <- NA
        } else {
          fit_ <- link_inverse(analysis_mod)(preds[[1]])
          bounds <- .msr_native_link_ci(analysis_mod, preds[[1]], preds[[2]],
            link_inverse(analysis_mod))
          lwr <- bounds$lower
          upr <- bounds$upper
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
      rel_poly <- poly_terms[vapply(poly_variables, function(names) v %in% names, logical(1))]
      subtitle_parts <- c(subtitle_parts,
                          paste0("Includes: ", paste(rel_poly, collapse = ", ")))
    }

    # Note the mean-conditioning for interaction terms
    if(v %in% inter_vars) {
      partners <- inter_variables[vapply(inter_variables,
        function(names) v %in% names, logical(1))]
      others <- setdiff(unique(unlist(partners, use.names = FALSE)), v)
      subtitle_parts <- c(subtitle_parts,
                          paste0("At mean of: ", paste(others, collapse = ", ")))
    }

    subtitle <- if(length(subtitle_parts) > 0) paste(subtitle_parts, collapse = "; ") else NULL

    p <- ggplot(df, aes(x = x, y = fit)) +
      xlab(df$variable[1]) +
      ylab(ylab) +
      theme_cowplot()
    if (all(is.na(df$lwr)) && all(is.na(df$upr)))
      p <- p + labs(caption = paste(c(
        if (native_mixed_response || inherits(analysis_mod, "HLfit"))
          "Random effects are set to zero, not integrated.",
        "Confidence limits unavailable; point predictions condition on spatial scales."), collapse = " "))
    else if (native_mixed_response)
      p <- p + labs(caption = paste(
        "Random effects are set to zero, not integrated.",
        "Intervals use native response-scale standard errors and condition on spatial scales."))
    else if (inherits(analysis_mod, "HLfit"))
      p <- p + labs(caption = paste(
        "Random effects are set to zero, not integrated.",
        "Limits transform native link-scale fixed-effect uncertainty and condition on spatial scales."))
    else if (!all(is.na(df$lwr)))
      p <- p + labs(caption = "Confidence limits condition on the selected spatial scales.")

    discrete <- !is.numeric(df$x)
    if (discrete) {
      p <- p + ggplot2::geom_point()
    } else {
      p <- p + geom_line(linewidth = 1)
    }

    if(!is.null(subtitle))
      p <- p + labs(subtitle = subtitle)

    if(!all(is.na(df$lwr)) && !all(is.na(df$upr))) {
      if (discrete) {
        p <- p + ggplot2::geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.15)
      } else {
        p <- p + geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25)
      }
    }

    p
  })

  names(gg_list) <- c_vars

  lapply(gg_list, print)
  invisible(gg_list)
}


# Recover the state/abundance formula without depending on coefficient names
# or unmarked's private formula helpers. Older fits retain a double formula;
# newer classes may retain named component formulas in a formlist slot.
.msr_unmarked_marginal_formula <- function(mod, type) {
  slots <- methods::slotNames(mod)
  if ("formlist" %in% slots && inherits(mod@formlist[[type]], "formula")) {
    return(mod@formlist[[type]])
  }
  if ("formula" %in% slots && inherits(mod@formula, "formula")) {
    f <- mod@formula
    if (length(f) == 3L) {
      return(stats::as.formula(call("~", f[[3]]), env = environment(f)))
    }
    if (length(f) == 2L) return(f)
  }
  stop("Could not recover the '", type, "' formula for this unmarked model.",
       call. = FALSE)
}


# Formula variables are data columns; factor contrasts are coefficient names.
# Grouping variables describe random effects, and offset-only variables need
# prediction defaults but do not describe a separate effect to plot.
# A quoted name may contain a pipe. Logical pipes inside I() or another
# fixed function are not formula random-effect terms.
.msr_marginal_random_call <- function(expr) {
  while (is.call(expr) && identical(expr[[1L]], as.name("(")) && length(expr) == 2L)
    expr <- expr[[2L]]
  is.call(expr) && is.symbol(expr[[1L]]) &&
    as.character(expr[[1L]]) %in% c("|", "||")
}

.msr_marginal_variables <- function(expr, offsets) {
  if (inherits(expr, c("terms", "formula"))) {
    tt <- stats::terms(expr)
    expressions <- lapply(attr(tt, "term.labels"), str2lang)
    if (offsets) {
      variables <- as.list(attr(tt, "variables"))[-1L]
      expressions <- c(expressions, variables[attr(tt, "offset")])
    }
    return(as.character(unique(unlist(lapply(expressions,
      .msr_marginal_variables, offsets = offsets), use.names = FALSE))))
  }
  if (!is.call(expr)) return(all.vars(expr))
  fun <- if (is.symbol(expr[[1L]])) as.character(expr[[1L]]) else ""
  if (.msr_marginal_random_call(expr) || (!offsets && fun == "offset"))
    return(character())
  # Only formula operators can introduce a random-effect term. Preserve every
  # raw input in a fixed function, including logical |/|| inside I().
  if (!fun %in% c("~", "+", "-", "*", "/", ":", "^", "%in%", "("))
    return(all.vars(expr))
  as.character(unique(unlist(lapply(as.list(expr)[-1L], .msr_marginal_variables,
    offsets = offsets), use.names = FALSE)))
}


# Reconstruct the native fixed training frame before site removal. Store only
# its learned prediction terms; never relearn a basis on a hypothetical grid.
.msr_unmarked_curve_basis <- function(component, data) {
  tt <- stats::terms(component, data = data)
  labels <- attr(tt, "term.labels")
  keep <- !vapply(lapply(labels, str2lang), .msr_marginal_random_call, logical(1))
  variables <- as.list(attr(tt, "variables"))[-1L]
  offsets <- vapply(variables[attr(tt, "offset")], deparse1, character(1), backtick = TRUE)
  labels <- c(labels[keep], offsets)
  intercept <- attr(tt, "intercept") == 1L
  form <- if (length(labels)) stats::reformulate(labels, intercept = intercept,
    env = attr(tt, ".Environment")) else stats::as.formula(if (intercept) "~1" else "~0",
      env = attr(tt, ".Environment"))
  frame <- stats::model.frame(form, data, na.action = stats::na.pass,
    drop.unused.levels = FALSE)
  fixed <- attr(frame, "terms")
  # R's prediction-call methods historically recognize bare function heads.
  # Recover metadata for these identical, explicitly qualified base functions
  # from the evaluated training columns, retaining their qualified calls.
  prediction <- as.list(attr(fixed, "predvars"))
  for (j in seq_along(frame)) {
    ex <- prediction[[j + 1L]]
    if (!is.call(ex) || !is.call(head <- ex[[1L]]) ||
        !identical(head[[1L]], as.name("::")) || length(head) != 3L) next
    known <- (identical(head[[2L]], as.name("stats")) && identical(head[[3L]], as.name("poly"))) ||
      (identical(head[[2L]], as.name("base")) && identical(head[[3L]], as.name("scale")))
    if (!known) next
    bare <- ex; bare[[1L]] <- head[[3L]]
    learned_call <- stats::makepredictcall(frame[[j]], bare)
    learned_call[[1L]] <- head
    prediction[[j + 1L]] <- learned_call
  }
  attr(fixed, "predvars") <- as.call(prediction)
  list(terms = fixed,
    learned = !identical(attr(fixed, "predvars"), attr(fixed, "variables")))
}


# Preserve explicit contrasts when native data-frame prediction drops them.
# Restrict the fallback to models whose state is a single linked estimate;
# derived-state predictions (for example ZIP abundance) need their own method.
.msr_unmarked_contrast_prediction <- function(mod, component, newdata, type) {
  supported <- identical(class(mod)[1], "unmarkedFitOccu") ||
    (identical(class(mod)[1], "unmarkedFitPCount") && mod@mixture %in% c("P", "NB"))
  if (!supported) {
    stop("Explicit nondefault factor contrasts cannot safely be predicted for ",
         "this unmarked class by plot_marginal_effects(). Use its fitted design ",
         "matrix and class-specific prediction method.", call. = FALSE)
  }
  terms <- stats::terms(component)
  random <- which(vapply(lapply(attr(terms, "term.labels"), str2lang),
    .msr_marginal_random_call, logical(1)))
  if (length(random)) terms <- stats::drop.terms(terms, random, keep.response = FALSE)
  frame <- stats::model.frame(terms, newdata, na.action = stats::na.pass)
  X <- stats::model.matrix(terms, frame)
  offset <- stats::model.offset(frame)
  coefficient_names <- names(unmarked::coef(mod, type = type, altNames = FALSE))
  if (!identical(sub("^\\(Intercept\\)$", "Int", colnames(X)), coefficient_names)) {
    stop("The prediction design does not match the fitted unmarked coefficients.",
         call. = FALSE)
  }
  # One row at a time avoids constructing a dense covariance matrix over the
  # whole plotting grid. Confidence bounds retain the fitted link function.
  rows <- lapply(seq_len(nrow(X)), function(i) {
    lc <- unmarked::linearComb(mod, X[i, , drop = FALSE], type = type,
      offset = if (is.null(offset)) NULL else offset[i], re.form = NA)
    prediction <- unmarked::backTransform(lc)
    ci <- unmarked::confint(prediction, level = 0.95)
    data.frame(Predicted = as.numeric(unmarked::coef(prediction)),
      SE = as.numeric(unmarked::SE(prediction)), lower = ci[1], upper = ci[2])
  })
  do.call(rbind, rows)
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
    .msr_attach_package(pkg)
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
    .msr_attach_package(pkg)
  }

  invisible(pkg)
}

.msr_attach_package <- function(pkg) {
  suppressWarnings(
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  )
}

utils::globalVariables(c("fit", "lwr", "upr"))


# Transform native mean limits only on a valid branch. Sorting handles a
# decreasing inverse link; it cannot repair a pole, invalid domain or bad SE.
.msr_native_link_ci <- function(model, fit, se, inverse) {
  eta_lower <- fit + qnorm(.025) * se
  eta_upper <- fit + qnorm(.975) * se
  valid <- is.finite(fit) & is.finite(se) & se >= 0 &
    is.finite(eta_lower) & is.finite(eta_upper)
  family <- tryCatch(stats::family(model), error = function(e) NULL)
  if (is.function(family$valideta)) for (i in which(valid))
    valid[i] <- isTRUE(family$valideta(c(eta_lower[i], eta_upper[i])))
  if (isTRUE(family$link %in% c("inverse", "1/mu^2")))
    valid <- valid & !(eta_lower <= 0 & eta_upper >= 0)
  if (identical(family$link, "sqrt")) valid <- valid & eta_lower >= 0
  lower <- upper <- rep(NA_real_, length(fit))
  index <- which(valid)
  if (length(index)) {
    a <- as.numeric(inverse(eta_lower[index]))
    b <- as.numeric(inverse(eta_upper[index]))
    acceptable <- is.finite(a) & is.finite(b)
    if (is.function(family$validmu)) for (j in which(acceptable))
      acceptable[j] <- isTRUE(family$validmu(c(a[j], b[j])))
    valid[index[!acceptable]] <- FALSE
    lower[index[acceptable]] <- pmin(a[acceptable], b[acceptable])
    upper[index[acceptable]] <- pmax(a[acceptable], b[acceptable])
  }
  if (any(!valid)) warning(sprintf(
    "Native mean confidence limits were not shown for %d row(s): nonfinite SEs or limits outside the native link/family domain.",
    sum(!valid)), call. = FALSE)
  list(lower = lower, upper = upper)
}
