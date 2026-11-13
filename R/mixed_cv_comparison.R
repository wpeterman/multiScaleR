.msr_mixed_cv_output <- function(folds, targets, settings, keep_fits) {
  list(predictions = list(), joint_scores = list(), fold_results = list(),
    training = list(), selection_records = list(), prediction_contracts = list(), diagnostics = list(),
    warnings = list(), errors = list(), fits = if (keep_fits) list() else NULL,
    folds = folds, targets = targets, settings = settings)
}

.msr_mixed_cv_payload <- function(x) unclass(x)[c("predictions", "joint_scores", "fold_results",
  "training", "selection_records", "prediction_contracts", "diagnostics", "warnings", "errors",
  "folds", "targets", "settings", "selection")]

.msr_mixed_cv_metric <- function(pred, joint, metric, family, complete) {
  value <- NA_real_
  count <- nrow(pred)
  weight <- if (metric == "brier") sum(pred$trials) else count
  unit <- switch(metric, neg_log_score = "joint_group_density_per_row",
    brier = "trial", if (family == "binomial") "proportion" else "response")
  if (complete) {
    observed <- pred$observed
    estimate <- pred$prediction
    response <- if (family == "binomial") observed / pred$trials else observed
    value <- switch(metric,
      neg_log_score = sum(joint$neg_log_score) / count,
      rmse = sqrt(mean((response - estimate)^2)),
      mae = mean(abs(response - estimate)),
      brier = sum(observed * (1 - estimate)^2 + (pred$trials - observed) * estimate^2) / weight)
  }
  data.frame(metric = metric, estimate = value, unit = unit, n_assessment = count,
    n_scored = sum(pred$eligible), weight_sum = weight, complete = complete)
}

.msr_mixed_cv_finish <- function(x) {
  x$fold_results <- do.call(rbind, x$fold_results)
  rownames(x$fold_results) <- NULL
  x$complete <- all(x$fold_results$status == "scored")
  metrics <- function(p, j, complete) do.call(rbind, lapply(x$settings$metrics,
    function(metric) .msr_mixed_cv_metric(p, j, metric, x$settings$family, complete)))
  x$metrics <- do.call(rbind, lapply(names(x$folds), function(name) {
    row <- match(name, x$fold_results$fold)
    cbind(fold = name, metrics(x$predictions[[name]], x$joint_scores[[name]],
      x$fold_results$status[row] == "scored"))
  }))
  rownames(x$metrics) <- NULL
  x$predictions <- do.call(rbind, x$predictions)
  rownames(x$predictions) <- NULL
  x$summary <- metrics(x$predictions, do.call(rbind, x$joint_scores), x$complete)
  ids <- x$targets$row_id
  x$assessment_counts <- stats::setNames(tabulate(match(x$predictions$row_id, ids), length(ids)), ids)
  x$integrity_sha256 <- .msr_gdr_hash(.msr_mixed_cv_payload(x))
  class(x) <- c("multiScaleR_mixed_cv", "multiScaleR_cv")
  x
}

.msr_mixed_cv_signature <- function(mod_list) {
  if (!is.list(mod_list) || !length(mod_list) || is.null(names(mod_list)) ||
      anyNA(names(mod_list)) || any(!nzchar(names(mod_list))) || anyDuplicated(names(mod_list)))
    stop("Supply uniquely named mixed-model CV candidates.", call. = FALSE)
  keys <- lapply(mod_list, function(x) {
    if (!inherits(x, "multiScaleR_mixed_cv") || !isTRUE(x$complete) ||
        any(x$fold_results$status != "scored"))
      stop("Compare complete mixed-model CV candidates. Every fold must be scored.", call. = FALSE)
    if (!identical(x$integrity_sha256, .msr_gdr_hash(.msr_mixed_cv_payload(x))))
      stop("Mixed-model CV scores, controls or recorded targets changed.", call. = FALSE)
    splits <- lapply(x$folds[sort(names(x$folds))], function(split)
      lapply(split[c("analysis", "assessment")], sort))
    used <- sort(unique(unlist(splits, use.names = FALSE)))
    targets <- x$targets[match(used, x$targets$row_id), , drop = FALSE]
    rownames(targets) <- NULL
    support <- switch(x$settings$family, gaussian = "continuous", Gamma = "continuous",
      poisson = "count", nbinom2 = "count", binomial = "binomial")
    list(splits = splits, targets = targets, support = support,
      prediction = x$settings$prediction)
  })
  if (!all(vapply(keys, identical, logical(1), keys[[1L]])))
    stop("Mixed-model CV comparison requires the same observation measure, group IDs, prediction target and splits.", call. = FALSE)
  invisible(keys[[1L]])
}

.msr_mixed_cv_tab <- function(mod_list, metric, by_fold) {
  .msr_mixed_cv_signature(mod_list)
  if (any(!vapply(mod_list, function(x) metric %in% x$settings$metrics, logical(1))))
    stop("Every mixed-model candidate must retain the requested metric.", call. = FALSE)
  folds <- if (by_fold) names(mod_list[[1L]]$folds) else "pooled"
  result <- do.call(rbind, lapply(folds, function(fold) {
    tab <- do.call(rbind, lapply(names(mod_list), function(name) {
      x <- mod_list[[name]]
      p <- if (by_fold) x$predictions[x$predictions$fold == fold, , drop = FALSE] else x$predictions
      j <- if (by_fold) x$joint_scores[[fold]] else do.call(rbind, x$joint_scores)
      score <- .msr_mixed_cv_metric(p, j, metric, x$settings$family, TRUE)
      data.frame(Modnames = name, Score = score$estimate, Metric = metric, Unit = score$unit,
        N_assessment = score$n_assessment, N_scored = score$n_scored, Weight_sum = score$weight_sum,
        Warnings = length(unlist(if (by_fold) x$warnings[fold] else x$warnings, use.names = FALSE)))
    }))
    tab$Rank <- rank(tab$Score, ties.method = "min")
    tab$Delta <- tab$Score - min(tab$Score)
    if (by_fold) tab$Fold <- fold
    tab[order(tab$Rank), c("Modnames", "Rank", "Score", "Delta", "Metric", "Unit",
      "N_assessment", "N_scored", "Weight_sum", "Warnings", if (by_fold) "Fold"), drop = FALSE]
  }))
  rownames(result) <- NULL
  result
}

.msr_mixed_cv_selection_record <- function(fit, ids) {
  model <- if (inherits(fit, "multiScaleR")) fit$opt_mod else fit
  if (!inherits(model, "glmmTMB") || !identical(rownames(stats::model.frame(model)), ids))
    stop("The native mixed fit does not match its analysis row IDs.", call. = FALSE)
  n <- stats::nobs(model)
  ll <- .msr_model_loglik(model)
  k <- .msr_parameter_count(model)
  total <- .msr_selection_parameter_count(fit, k)
  if (!identical(as.integer(n), length(ids)) || any(!is.finite(c(ll, k, total))) ||
      k < 0 || total < k || any(c(k, total) != floor(c(k, total))))
    stop("The mixed-model training likelihood, parameter counts or nobs are invalid.", call. = FALSE)
  list(logLik = ll, K = total, native_K = k, ecological_K = total - k,
    nobs = length(ids), analysis = ids)
}

.msr_mixed_select_cv <- function(mod_list, criterion, ties, keep_fits) {
  .msr_mixed_cv_signature(mod_list)
  first <- mod_list[[1L]]
  if (!all(vapply(mod_list, function(x) setequal(x$settings$metrics, first$settings$metrics), logical(1))))
    stop("Mixed-model selection candidates must retain the same metrics.", call. = FALSE)
  if (keep_fits && any(vapply(mod_list, function(x) is.null(x$fits), logical(1))))
    stop("Retaining selected fits requires keep_fits = TRUE in every candidate.", call. = FALSE)
  out <- .msr_mixed_cv_output(first$folds, first$targets, first$settings, keep_fits)
  candidates <- names(mod_list)
  out$selection <- list(criterion = criterion, ties = ties, candidates = candidates,
    chosen = stats::setNames(rep(NA_character_, length(first$folds)), names(first$folds)),
    tables = list(), candidate_warnings = lapply(mod_list, function(x) x$warnings),
    candidate_diagnostics = lapply(mod_list, function(x) x$diagnostics))
  out$settings$backend <- "training_selection"
  out$settings$selection_uncertainty <- "not_integrated"
  for (name in names(first$folds)) {
    winner <- NULL
    error <- tryCatch({
      records <- lapply(mod_list, function(x) {
        record <- x$selection_records[[name]]
        if (!is.null(x$fits) && !identical(record,
            .msr_mixed_cv_selection_record(x$fits[[name]], x$folds[[name]]$analysis)))
          stop("Mixed-model selection records disagree with retained fits.", call. = FALSE)
        record
      })
      k <- vapply(records, function(x) x$K, numeric(1))
      ll <- vapply(records, function(x) x$logLik, numeric(1))
      n <- vapply(records, function(x) x$nobs, numeric(1))
      if (criterion == "AICc" && any(n <= k + 1))
        stop("Training mixed-model AICc requires n > K + 1.", call. = FALSE)
      value <- switch(criterion, AIC = -2 * ll + 2 * k, BIC = -2 * ll + log(n) * k,
        AICc = -2 * ll + 2 * k + 2 * k * (k + 1) / (n - k - 1))
      if (any(!is.finite(value))) stop("Training information criteria must all be finite.", call. = FALSE)
      tab <- data.frame(Modnames = candidates, K = k, logLik = ll, criterion_value = value)
      names(tab)[4L] <- criterion
      out$selection$tables[[name]] <- tab
      best <- candidates[value == min(value)]
      if (length(best) > 1L && ties == "error")
        stop("Exact training-criterion tie. Predeclare ties = 'first' to resolve it.", call. = FALSE)
      winner <- best[1L]
      NULL
    }, error = function(e) list(stage = "training_selection", message = conditionMessage(e)))
    x <- if (is.null(error)) mod_list[[winner]] else first
    p <- x$predictions[x$predictions$fold == name, , drop = FALSE]
    row <- x$fold_results[x$fold_results$fold == name, , drop = FALSE]
    if (is.null(error)) {
      out$selection$chosen[name] <- winner
      for (item in c("joint_scores", "training", "selection_records", "prediction_contracts", "diagnostics"))
        out[[item]][name] <- x[[item]][name]
      if (keep_fits) out$fits[name] <- x$fits[name]
      out$settings$fold_settings[[name]] <- x$settings
    } else {
      p$prediction <- NA_real_
      p$eligible <- FALSE
      row$status <- "failed"
    }
    out$predictions[[name]] <- p
    out$fold_results[[name]] <- row
    out$warnings[name] <- x$warnings[name]
    out$errors[name] <- list(error)
  }
  .msr_mixed_cv_finish(out)
}
