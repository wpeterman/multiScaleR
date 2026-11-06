# Conditional-choice CV has a set-level target, distinct from binary GLM loss.
.msr_choice_formula <- function(formula, data, spatial = character()) {
  if (!inherits(formula, "formula") || length(formula) != 3L || !is.symbol(formula[[2L]]))
    stop("Use a direct zero/one case column and one bare strata(choice_id) term.", call. = FALSE)
  # Expanding dot uses declared column names, never assessment values or a
  # fitted preprocessing transform. Missing spatial columns are placeholders.
  columns <- data
  for (nm in setdiff(spatial, names(columns))) columns[[nm]] <- 0
  tt <- stats::terms(formula, specials = "strata", data = columns)
  ex <- lapply(attr(tt, "term.labels"), str2lang)
  contains <- function(x, forbidden) {
    if (!is.call(x)) return(FALSE)
    head <- x[[1L]]
    if (is.call(head) && length(head) == 3L && as.character(head[[1L]]) %in% c("::", ":::")) head <- head[[3L]]
    (is.symbol(head) && as.character(head) %in% forbidden) ||
      any(vapply(as.list(x)[-1L], contains, logical(1), forbidden = forbidden))
  }
  sg <- which(vapply(ex, function(x) is.call(x) && identical(x[[1L]], as.name("strata")), logical(1)))
  if (length(sg) != 1L || length(ex[[sg]]) != 2L || !is.symbol(ex[[sg]][[2L]]) ||
      any(vapply(ex[-sg], contains, logical(1), forbidden = "strata")) ||
      any(vapply(ex, contains, logical(1), forbidden = c("cluster", "frailty", "tt", "pspline"))))
    stop("Choice CV requires one bare strata(column), without special/penalized Cox terms.", call. = FALSE)
  response <- as.character(formula[[2L]]); choice <- as.character(ex[[sg]][[2L]])
  if (!all(c(response, choice) %in% names(data)) || response %in% spatial || choice %in% spatial ||
      length(setdiff(all.vars(formula), c(names(data), spatial, "."))))
    stop("Supply explicit case, choice-set and formula variables; case/set IDs cannot be spatial predictors.", call. = FALSE)
  y <- data[[response]]; group <- as.character(data[[choice]])
  if ((!is.numeric(y) && !is.logical(y)) || !is.null(dim(y)) || anyNA(y) || any(!y %in% c(0, 1)) ||
      anyNA(group) || any(!nzchar(group)))
    stop("Choice CV requires finite zero/one cases and nonmissing, nonempty set IDs.", call. = FALSE)
  .msr_clogit_groups(group, as.numeric(y))
  list(response = response, choice = choice, formula = stats::formula(tt),
    targets = data.frame(row_id = rownames(data), observed = as.numeric(y), trials = 1,
      weight = 1, choice = group, stringsAsFactors = FALSE))
}

.msr_choice_splits <- function(folds, targets) {
  folds <- .msr_cv_folds(folds, targets$row_id)
  groups <- split(targets$row_id, targets$choice)
  for (name in names(folds)) {
    status <- stats::setNames(rep("excluded", nrow(targets)), targets$row_id)
    status[folds[[name]]$analysis] <- "analysis"; status[folds[[name]]$assessment] <- "assessment"
    if (any(vapply(groups, function(ids) length(unique(status[ids])) != 1L, logical(1))))
      stop("Every choice set must stay complete in analysis, assessment or exclusions for fold `", name, "`.", call. = FALSE)
  }
  folds
}

# The saved terms carry training polynomials/splines. Removing strata permits
# new complete sets without inventing an effect for their unseen labels.
.msr_choice_lp <- function(model, data) {
  native <- .msr_fixed_prediction_terms(stats::terms(model), model$model)
  labels <- attr(native, "term.labels")
  sg <- which(vapply(lapply(labels, str2lang), function(x)
    is.call(x) && identical(x[[1L]], as.name("strata")), logical(1)))
  if (length(sg) != 1L) stop("Cannot recover one native choice stratum.", call. = FALSE)
  offsets <- vapply(attr(native, "offset"), function(i)
    paste(deparse(attr(native, "variables")[[i + 1L]]), collapse = ""), "")
  formula <- stats::reformulate(c(labels[-sg], offsets), intercept = attr(native, "intercept") == 1L,
    env = environment(native))
  tt <- .msr_clogit_prediction_terms(native, formula)
  frame <- stats::model.frame(tt, data, na.action = stats::na.fail,
    xlev = model$xlevels[setdiff(names(model$xlevels), labels[sg])])
  X <- stats::model.matrix(tt, frame, contrasts.arg = model$contrasts)
  beta <- stats::coef(model)
  if (any(!is.finite(beta)) || any(!names(beta) %in% colnames(X)))
    stop("Native choice coefficients/design are not finite and aligned.", call. = FALSE)
  eta <- if (length(beta)) drop(X[, names(beta), drop = FALSE] %*% beta) else rep(0, nrow(data))
  offset <- stats::model.offset(frame)
  if (!is.null(offset)) eta <- eta + offset
  if (length(eta) != nrow(data) || any(!is.finite(eta))) stop("Invalid held-out choice linear predictors.", call. = FALSE)
  as.numeric(eta)
}

.msr_choice_native_frame <- function(fit) {
  m <- .analysis_model(if (inherits(fit, "multiScaleR")) fit$opt_mod else fit)
  # clogit can replace coxph.penal classes; retain native penalty markers too.
  if (!inherits(m, "coxph") || inherits(m, c("coxph.penal", "coxphms")) ||
      !is.data.frame(m$model) || !is.null(m$penalty) ||
      any(vapply(m$model, inherits, logical(1), "coxph.penalty")) ||
      !m$method %in% c("exact", "efron", "breslow") ||
      !is.null(m$na.action) || !is.null(m$call$subset) || !is.null(m$weights) ||
      !is.null(m$call$weights) || !is.null(m$naive.var) || !is.null(m$call$cluster) || !is.null(m$call$offset))
    stop("Require a complete, unweighted, unpenalized native conditional-logistic training fit.", call. = FALSE)
  y <- stats::model.response(m$model); tt <- stats::terms(m)
  labels <- attr(tt, "term.labels")
  sg <- which(vapply(lapply(labels, str2lang), function(x)
    is.call(x) && identical(x[[1L]], as.name("strata")), logical(1)))
  if (!inherits(y, "Surv") || ncol(y) != 2L || !identical(attr(y, "type"), "right") ||
      any(!is.finite(y)) || any(y[, 1L] != 1) || length(sg) != 1L ||
      any(lengths(attr(tt, "specials")[setdiff(names(attr(tt, "specials")), "strata")]) > 0L))
    stop("Require one constant-time, one-case conditional choice stratum.", call. = FALSE)
  group <- as.character(m$model[[labels[sg]]])
  .msr_clogit_groups(group, y[, 2L])
  list(model = m, targets = data.frame(row_id = rownames(m$model), observed = as.numeric(y[, 2L]),
    trials = 1, weight = 1, choice = group, stringsAsFactors = FALSE))
}

# survival stores a scalar null loglik but its ordinary logLik method can return
# NA. This CV-only subclass checks that scalar against the exact offset target.
#' @exportS3Method stats::logLik
#' @noRd
logLik.msr_cv_choice_null <- function(object, ...) {
  native <- .msr_choice_native_frame(object); t <- native$targets
  if (length(stats::coef(object)) || length(object$loglik) != 1L || !is.finite(object$loglik))
    stop("Missing native scalar zero-coefficient choice likelihood.", call. = FALSE)
  offset <- stats::model.offset(object$model)
  if (is.null(offset)) offset <- rep(0, nrow(t))
  z <- .msr_clogit_softmax(as.numeric(offset), .msr_clogit_groups(t$choice, t$observed))
  ll <- sum(z$log_probability[t$observed == 1])
  if (!is.finite(ll) || abs(ll - object$loglik) > 1e-8 * max(1, abs(ll), abs(object$loglik)))
    stop("Native null choice likelihood disagrees with its complete-set offset target.", call. = FALSE)
  structure(as.numeric(object$loglik), df = 0L, nobs = sum(t$observed), class = "logLik")
}

.msr_choice_fit <- function(formula, data, method, control) {
  coxph <- survival::coxph
  env <- new.env(parent = environment(formula)); env$strata <- survival::strata; env$Surv <- survival::Surv
  environment(formula) <- env
  m <- do.call(survival::clogit, list(formula = formula, data = data, method = method, control = control,
    na.action = stats::na.fail, model = TRUE, x = TRUE, y = TRUE))
  if (!length(stats::coef(m))) {
    class(m) <- c("msr_cv_choice_null", class(m))
    stats::logLik(m) # validate before permitting a scored fold
  }
  if (any(!is.finite(stats::coef(m))) || !is.finite(as.numeric(stats::logLik(m))))
    stop("Native training choice fit has nonfinite coefficients or likelihood.", call. = FALSE)
  m
}

.msr_choice_identity <- function(fit, split, targets, method = NULL) {
  native <- .msr_choice_native_frame(fit); actual <- native$targets
  expected <- targets[match(split$analysis, targets$row_id), , drop = FALSE]
  partition <- function(x) match(x, unique(x))
  if (!identical(actual$row_id, split$analysis) || !identical(actual$observed, expected$observed) ||
      !identical(partition(actual$choice), partition(expected$choice)))
    stop("Native training choice responses, strata or analysis row IDs disagree with targets.", call. = FALSE)
  if (!is.null(method) && !identical(native$model$method, method))
    stop("Native training partial-likelihood method disagrees with CV settings.", call. = FALSE)
  invisible(TRUE)
}

.msr_choice_score <- function(p, metrics) {
  if (!identical(metrics, "neg_log_score") || any(!is.finite(p$linear_predictor)))
    stop("Choice CV requires finite linear predictors and neg_log_score.", call. = FALSE)
  groups <- unlist(lapply(split(seq_len(nrow(p)), p$fold), function(i)
    lapply(.msr_clogit_groups(p$choice[i], p$observed[i]), function(g) i[g])), recursive = FALSE)
  z <- .msr_clogit_softmax(p$linear_predictor, groups)
  p$prediction <- z$probability; p$log_probability <- z$log_probability
  p$eligible <- p$observed == 1; p$neg_log_score <- ifelse(p$eligible, -z$log_probability, NA_real_)
  p
}

.msr_choice_metric <- function(p, complete) {
  keep <- p$eligible & !is.na(p$neg_log_score)
  data.frame(metric = "neg_log_score", estimate = if (complete && any(keep)) mean(p$neg_log_score[keep]) else NA_real_,
    unit = "choice_set", n_assessment = nrow(unique(p[c("fold", "choice")])),
    n_scored = sum(keep), weight_sum = sum(keep), complete = complete, n_alternatives = nrow(p), stringsAsFactors = FALSE)
}

#' Cross-Validate Conditional Choices in Complete Matched Sets
#'
#' Refit native conditional-logistic coefficients and, optionally, ecological
#' scales within explicit training folds. Score the probability of the used
#' alternative among every alternative in each held-out choice set.
#'
#' @param formula Two-sided formula with a direct numeric/logical zero/one case
#'   column and one bare \code{strata(choice_id)} term. Fixed factors,
#'   interactions, offsets and supported learned transformations use native terms.
#' @param data Data frame containing explicit formula variables and unique row
#'   IDs matching \code{folds} and, when supplied, prepared spatial inputs.
#' @param folds Explicit named analysis/assessment row-ID splits. Every choice
#'   set must remain wholly in training, assessment or excluded rows in each split.
#'   Use whole-animal groups when assessment concerns transfer to new animals.
#' @param kernel_inputs Optional raw-cell \code{multiScaleR_data}. Default
#'   \code{NULL} fits an ordinary native reference. Otherwise scales and learned
#'   normalization/bin references are refit only on analysis rows.
#' @param backend Spatial backend, \code{"conditional_logit"} (default) or
#'   \code{"refit"}. Existing kernel, metric and numerical-admission guards apply.
#' @param method Native one-case partial-likelihood method, \code{"exact"}
#'   (default), \code{"efron"} or \code{"breslow"}.
#' @param cox_control Named native \code{survival::coxph.control()} settings.
#' @inheritParams cross_validate
#' @return A \code{multiScaleR_choice_cv} object inheriting \code{multiScaleR_cv}.
#'   Standard CV fields retain their meanings with the following target-specific
#'   definitions. \code{predictions} retains every alternative, its \code{choice}
#'   ID, zero/one \code{observed} case, uncentered \code{linear_predictor},
#'   within-set \code{prediction} probability and natural-log \code{log_probability}.
#'   Only used alternatives are \code{eligible}; their \code{neg_log_score} is
#'   the negative log choice probability. \code{summary}/\code{metrics} average
#'   over choice sets, with \code{unit = "choice_set"}, set counts and a separate
#'   \code{n_alternatives} count. \code{fold_results} counts alternative rows.
#'   \code{targets} binds all used row IDs, cases and set membership.
#'   \code{training}, \code{diagnostics}, \code{warnings}, \code{errors},
#'   \code{raster_groups}, \code{fits} and opt-in \code{selection_records} retain
#'   their contracts described in \code{\link{cross_validate}}. Choice records
#'   use version 2, native partial-method/stratum signatures and event-based
#'   \code{nobs}. Failed folds retain every planned alternative and pooled
#'   estimates remain \code{NA}. Inspect original diagnostics before comparison.
#' @details Each set needs one used alternative and at least two alternatives.
#'   Weights, missing-row omission, robust clustering, penalized/special Cox terms
#'   and generic censored survival responses are outside this adapter. No rows
#'   are silently dropped. Predictions use saved training terms without assigning
#'   an effect to new set labels. Log scores are calculated directly in log space,
#'   including when a displayed probability underflows to zero.
#'
#'   A zero-coefficient native fit receives a CV-only likelihood adapter when
#'   survival retains a scalar null likelihood. Its value must agree with the
#'   exact complete-set offset likelihood before native K = 0 and event counts
#'   can be used. Other native Cox fits are not altered. Inconsistent or missing
#'   native metadata remains a failure.
#'   Every scored fold checks its unweighted, unpenalized native fit, analysis
#'   rows, complete strata and partial method before prediction. These checks
#'   apply even when fits and selection records are not retained.
#'
#'   Use \code{\link{cv_tab}} to compare this choice target and
#'   \code{\link{select_cv}} for training-only AIC/AICc/BIC selection. Penalties
#'   retain native event counts and full coefficient/scale/shape K. Ordinary
#'   Bernoulli probabilities and relative selection risks are different targets.
#'   Dot formulas expand declared column names before fitting and learn no
#'   assessment transformations. Exclude the set ID as an ordinary predictor,
#'   for example \code{case_ ~ . - choice_id + strata(choice_id)}. Predictors
#'   constant within every set have no identified conditional coefficient.
#'   The supported training-only transformations described in
#'   \code{\link{cross_validate}} also apply inside choice expressions and offsets.
#'   Their saved constants do not introduce stratum effects or change event K.
#'   Point scores do not integrate parameter/selection uncertainty, provide
#'   observation intervals or establish calibration. This adapter does not
#'   promote guarded spatial-scale or Bayesian methods.
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   d <- data.frame(case_ = rep(c(1, 0, 0), 6),
#'                   step = rep(1:6, each = 3), x = sin(1:18))
#'   f <- grouped_cv_folds(d$step, k = 2)
#'   cv <- cross_validate_clogit(case_ ~ x + strata(step), d, f,
#'                               keep_selection = TRUE)
#'   cv$summary # natural-log loss per held-out choice set
#' }
#' @seealso \code{\link{cross_validate}}, \code{\link{grouped_cv_folds}},
#'   \code{\link{spatial_cv_folds}}, \code{\link{cv_tab}}, \code{\link{select_cv}}
#' @export
cross_validate_clogit <- function(formula, data, folds, kernel_inputs = NULL,
    backend = c("conditional_logit", "refit"), method = c("exact", "efron", "breslow"),
    cox_control = survival::coxph.control(), optim_args = list(), keep_fits = FALSE, keep_selection = FALSE) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("Choice CV requires survival.", call. = FALSE)
  backend <- match.arg(backend); method <- match.arg(method)
  validate_scalar_logical(keep_fits, "keep_fits"); validate_scalar_logical(keep_selection, "keep_selection")
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data))) stop("Supply nonempty data with unique columns.", call. = FALSE)
  ids <- .msr_fold_ids(rownames(data), nrow(data)); inputs <- kernel_inputs
  spatial <- if (is.null(inputs)) character() else inputs$scale_vars$covariate
  specification <- .msr_choice_formula(formula, data, spatial)
  formula <- specification$formula
  folds <- .msr_choice_splits(folds, specification$targets)
  covariates <- intersect(spatial, all.vars(formula))
  if (!is.null(inputs)) {
    if (!inherits(inputs, "multiScaleR_data") || !is.list(inputs$d_list) || !is.list(inputs$raw_cov) ||
        length(inputs$raw_cov) != nrow(inputs$kernel_dat) || length(inputs$d_list) != nrow(inputs$kernel_dat) ||
        !.msr_valid_observation_ids(rownames(inputs$kernel_dat), nrow(inputs$kernel_dat)) ||
        !setequal(ids, rownames(inputs$kernel_dat)) || !length(covariates) ||
        !any(inputs$scale_vars$optimize[match(covariates, spatial)]))
      stop("Spatial choice CV requires matching raw-cell inputs and an optimized formula predictor.", call. = FALSE)
    if (!is.null(inputs$raster_group) && (!is.character(inputs$raster_group) ||
        !identical(names(inputs$raster_group), rownames(inputs$kernel_dat)) ||
        anyNA(inputs$raster_group) || any(!nzchar(inputs$raster_group))))
      stop("raster_group must align exactly with prepared row IDs.", call. = FALSE)
  }
  allowed <- setdiff(names(formals(multiScale_optim)), c("fitted_mod", "kernel_inputs", "join_by", "n_cores", "PSOCK", "refit_fn", "backend", "verbose"))
  if (!is.list(optim_args) || (length(optim_args) && (is.null(names(optim_args)) || anyNA(names(optim_args)) ||
      anyDuplicated(names(optim_args)) || any(!names(optim_args) %in% allowed))) ||
      (!is.null(optim_args$profile_control) && !is.list(optim_args$profile_control)) ||
      !is.null(optim_args$profile_control$worker_pool))
    stop("optim_args must name supported serial scale controls without data/worker overrides.", call. = FALSE)
  if (!is.list(cox_control) || is.null(names(cox_control)) || anyDuplicated(names(cox_control)) ||
      any(!names(cox_control) %in% names(formals(survival::coxph.control))))
    stop("cox_control must be a named coxph.control() list.", call. = FALSE)
  cox_control <- do.call(survival::coxph.control, cox_control)
  actual_backend <- if (is.null(inputs)) "clogit" else backend
  settings <- list(family = "conditional_choice", link = "softmax", method = method, backend = actual_backend,
    cox_control = cox_control, optim_args = optim_args, metrics = "neg_log_score", normalization = "training_only",
    prediction_target = "within_set_choice_probability", assessment_unit = "choice_set", recorded_row_unit = "alternative",
    parameter_uncertainty = "not_integrated", observation_intervals = "not_provided", calibration = "not_established", workers = 1L)
  output <- list(predictions = list(), metrics = list(), fold_results = list(), training = list(), diagnostics = list(),
    warnings = list(), errors = list(), prediction_contracts = list(), fits = if (keep_fits) list() else NULL,
    selection_records = if (keep_selection) list() else NULL)
  for (name in names(folds)) {
    split <- folds[[name]]; ai <- match(split$analysis, ids); ti <- match(split$assessment, ids)
    pred <- specification$targets[ti, , drop = FALSE]; pred$fold <- name
    pred$linear_predictor <- pred$prediction <- pred$log_probability <- pred$neg_log_score <- NA_real_; pred$eligible <- FALSE
    initial <- fit <- NULL; messages <- character(); stage <- "preparation"
    error <- tryCatch(withCallingHandlers({
      training <- data[ai, , drop = FALSE]; assessment <- data[ti, , drop = FALSE]
      if (!is.null(inputs)) {
        train <- .msr_cv_inputs(inputs, split$analysis, covariates)
        raw <- .msr_cv_initial_values(train); center <- colMeans(raw); spread <- apply(raw, 2L, stats::sd)
        if (any(!is.finite(spread) | spread <= 0)) stop("Training spatial predictors require positive finite spread.")
        train$kernel_dat <- as.data.frame(sweep(sweep(raw, 2L, center), 2L, spread, "/"), optional = TRUE)
        train$scl_params <- list(mean = center, sd = spread)
        for (nm in covariates) training[[nm]] <- train$kernel_dat[[nm]]
        output$training[[name]] <- list(initial_normalization = train$scl_params,
          bin_reference = if (is.null(train$binned)) NULL else train$binned$cv_reference)
      }
      stage <- "training_clogit"
      initial <- .msr_choice_fit(formula, training, method, cox_control); fit <- initial
      output$training[[name]]$initial_coefficients <- stats::coef(initial)
      if (length(messages)) stop("Native training choice fit emitted warnings; inspect the retained warnings.")
      if (!is.null(inputs)) {
        stage <- "scale_fit"
        fit <- NULL
        fit <- do.call(multiScale_optim, c(list(fitted_mod = initial, kernel_inputs = train, n_cores = 1L,
          backend = backend, verbose = FALSE), optim_args))
        if (!identical(as.integer(fit$optim_results$convergence), 0L)) stop("Training scale optimizer did not converge.")
        output$training[[name]] <- c(output$training[[name]], list(fitted_normalization = fit$scl_params,
          coefficients = stats::coef(.analysis_model(fit$opt_mod)), scale_est = fit$scale_est, shape_est = fit$shape_est,
          optimizer = fit$optim_results[c("par", "value", "convergence", "message")], start_results = fit$profile_fit$start_results))
      }
      stage <- "training_identity"
      .msr_choice_identity(fit, split, specification$targets, method)
      stage <- "prediction"
      fit <- .msr_cv_prediction_metadata(fit, training, if (is.null(inputs)) NULL else train)
      if (is.null(inputs)) pred$linear_predictor <- .msr_choice_lp(fit, assessment) else {
        assess <- .msr_cv_inputs(inputs, split$assessment, covariates, if (is.null(train$binned)) NULL else train$binned$cv_reference)
        raw <- .msr_cv_initial_values(assess)
        assess$kernel_dat <- as.data.frame(sweep(sweep(raw, 2L, center), 2L, spread, "/"), optional = TRUE); assess$scl_params <- train$scl_params
        for (nm in covariates) assessment[[nm]] <- assess$kernel_dat[[nm]]
        if (backend == "conditional_logit") pred$linear_predictor <- as.numeric(stats::predict(fit, assessment,
          kernel_inputs = assess, type = "link", se.fit = FALSE, uncertainty = "conditional", interval = "none")) else {
          nd <- .msr_refit_prediction_data(fit, assessment, assess, join_by = NULL, include_fixed = TRUE)
          pred$linear_predictor <- .msr_choice_lp(.analysis_model(fit$opt_mod), nd)
        }
      }
      stage <- "scoring"; pred <- .msr_choice_score(pred, "neg_log_score"); NULL
    }, warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) list(stage = stage, message = conditionMessage(e), start_results = e$start_results))
    complete <- is.null(error)
    if (!complete) {
      pred$linear_predictor <- pred$prediction <- pred$log_probability <- pred$neg_log_score <- NA_real_; pred$eligible <- FALSE
    }
    output$predictions[[name]] <- pred
    output$metrics[[name]] <- cbind(fold = name, .msr_choice_metric(pred, complete))
    output$fold_results[[name]] <- data.frame(fold = name, n_analysis = length(ai), n_assessment = length(ti),
      status = if (complete) "scored" else "failed", stringsAsFactors = FALSE)
    output$warnings[name] <- list(messages); output$errors[name] <- list(error)
    native <- if (is.null(fit)) NULL else .analysis_model(if (inherits(fit, "multiScaleR")) fit$opt_mod else fit)
    output$diagnostics[name] <- list(list(initial = if (is.null(initial)) NULL else list(iter = initial$iter, method = initial$method),
      native = if (is.null(native)) NULL else list(iter = native$iter, method = native$method),
      scale = if (inherits(fit, "multiScaleR")) fit$diagnostics else NULL))
    output$prediction_contracts[[name]] <- .msr_prediction_contract("conditional_likelihood", actual_backend,
      "choice_probability", "response", population = "fixed_at_estimate",
      kernel = if (is.null(inputs)) "not_applicable" else "fixed_at_estimate", conditioning = "training_fit_and_complete_assessment_choice_sets",
      effects = list(choice_sets = "supplied_complete_choice_sets"),
      omitted = c("parameter_uncertainty", "selection_uncertainty", "observation_noise"), fixed_inputs = TRUE)
    if (!is.null(inputs$raster_group)) output$raster_groups[[name]] <- lapply(split, function(rows) inputs$raster_group[rows])
    if (keep_fits) output$fits[name] <- list(fit)
    if (keep_selection) output$selection_records[name] <- list(.msr_choice_record(fit, output$diagnostics[[name]],
      .msr_cv_selection_controls(settings), split, specification$targets))
  }
  output <- .msr_cv_finish(output, folds, ids, "conditional_choice", "neg_log_score")
  used <- sort(unique(unlist(lapply(folds, function(f) c(f$analysis, f$assessment)), use.names = FALSE)))
  output$targets <- specification$targets[match(used, specification$targets$row_id), , drop = FALSE]; rownames(output$targets) <- NULL
  output$settings <- settings; output$call <- match.call(); output
}

#' @rdname cross_validate_clogit
#' @param x A conditional-choice CV result.
#' @param ... Reserved for future print controls.
#' @export
print.multiScaleR_choice_cv <- function(x, ...) {
  cat("Conditional-choice CV:", sum(x$fold_results$status == "scored"), "of", nrow(x$fold_results), "folds scored.\n")
  if (!x$complete) cat("Incomplete: inspect retained errors, warnings and diagnostics.\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}
