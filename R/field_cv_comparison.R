# Compact score metadata is bound without retaining native models or dense
# conditional covariance arrays. This is an integrity check, not provenance
# authentication or a substitute for independent numerical qualification.
.msr_field_cv_payload <- function(x) unclass(x)[c("predictions", "joint_scores", "fold_results",
  "training", "selection_records", "diagnostics", "warnings", "errors",
  "prediction_contracts", "targets", "coordinates", "folds", "settings", "selection")]

.msr_field_cv_metric <- function(x, metric, fold = NULL) {
  p <- x$predictions
  if (!is.null(fold)) p <- p[p$fold == fold, , drop = FALSE]
  groups <- if (is.null(fold)) names(x$folds) else fold
  value <- if (metric == "rmse") sqrt(mean(p$squared_error)) else
    if (metric == "mae") mean(p$absolute_error) else
      sum(vapply(x$joint_scores[groups], function(v) v$neg_log_score, 0)) / nrow(p)
  list(estimate = value, n = nrow(p), unit = if (metric == "neg_log_score")
    "joint_block_density_per_observation" else "response")
}

.msr_field_cv_finish <- function(x, ids) {
  x$predictions <- do.call(rbind, x$predictions); rownames(x$predictions) <- NULL
  x$fold_results <- do.call(rbind, x$fold_results); rownames(x$fold_results) <- NULL
  x$complete <- all(x$fold_results$status == "scored")
  x$assessment_counts <- stats::setNames(tabulate(match(x$predictions$row_id, ids), length(ids)), ids)
  x$summary <- do.call(rbind, lapply(x$settings$metrics, function(metric) {
    value <- if (x$complete) .msr_field_cv_metric(x, metric) else
      list(estimate = NA_real_, n = nrow(x$predictions), unit = if (metric == "neg_log_score")
        "joint_block_density_per_observation" else "response")
    data.frame(metric = metric, estimate = value$estimate, n_assessment = value$n,
      unit = value$unit, complete = x$complete)
  }))
  x$integrity_sha256 <- .msr_gdr_hash(.msr_field_cv_payload(x))
  class(x) <- c("multiScaleR_field_cv", "multiScaleR_cv")
  x
}

.msr_field_cv_comparison <- function(mod_list) {
  if (!is.list(mod_list) || !length(mod_list) || is.null(names(mod_list)) ||
      anyNA(names(mod_list)) || any(!nzchar(names(mod_list))) || anyDuplicated(names(mod_list)))
    stop("Supply a uniquely named nonempty list of spatial-field CV candidates.", call. = FALSE)
  signatures <- lapply(mod_list, function(x) {
    if (!inherits(x, "multiScaleR_field_cv") || !isTRUE(x$complete) ||
        !is.data.frame(x$fold_results) || any(x$fold_results$status != "scored"))
      stop("Compare complete matching spatial-field CV targets; failed folds and independent GLM scores cannot be mixed or dropped.", call. = FALSE)
    if (!identical(x$integrity_sha256, .msr_gdr_hash(.msr_field_cv_payload(x))))
      stop("Recorded field CV targets, controls, scores or training metadata changed.", call. = FALSE)
    if (!identical(x$settings$family, "gaussian") ||
        !identical(x$settings$prediction_target, "held_out_Gaussian_response_joint_block") ||
        !identical(x$settings$parameter_uncertainty, "not_integrated") ||
        !identical(x$settings$ecological_scale_uncertainty, "not_integrated") ||
        !isTRUE(x$settings$observation_noise_included))
      stop("Unqualified spatial-field CV scoring or uncertainty target.", call. = FALSE)
    t <- x$targets; xy <- x$coordinates
    if (!is.data.frame(t) || !identical(names(t), c("row_id", "observed")) ||
        anyDuplicated(t$row_id) || anyNA(t$row_id) || !all(is.finite(t$observed)) ||
        !is.data.frame(xy) || !identical(names(xy), c("row_id", "x", "y")) ||
        !identical(xy$row_id, t$row_id) || any(!is.finite(as.matrix(xy[c("x", "y")]))))
      stop("Invalid retained field CV observation/coordinate identities.", call. = FALSE)
    .msr_cv_folds(x$folds, t$row_id)
    for (name in names(x$folds)) {
      split <- x$folds[[name]]; p <- x$predictions[x$predictions$fold == name, , drop = FALSE]
      if (!identical(p$row_id, split$assessment) || !all(p$eligible) || any(!is.finite(p$prediction)) ||
          !identical(as.numeric(p$observed), as.numeric(t$observed[match(p$row_id, t$row_id)])) ||
          !isTRUE(all.equal(p$squared_error, (p$observed - p$prediction)^2, tolerance = 0)) ||
          !isTRUE(all.equal(p$absolute_error, abs(p$observed - p$prediction), tolerance = 0)))
        stop("Recorded field assessment rows or mean losses disagree with targets.", call. = FALSE)
      j <- x$joint_scores[[name]]
      if (!is.list(j) || !identical(j$n, nrow(p)) || !is.finite(j$quadratic) || j$quadratic < 0 ||
          !is.finite(j$log_determinant) || !is.finite(j$neg_log_score) ||
          !isTRUE(all.equal(j$neg_log_score,
            .5 * (j$n * log(2 * pi) + j$log_determinant + j$quadratic), tolerance = 0)))
        stop("Recorded joint Gaussian block density is inconsistent.", call. = FALSE)
      d <- x$selection_records[[name]]
      if (!is.list(d) || !identical(d$training_ids, split$analysis) ||
          !identical(d$training_response_sha256, .msr_gdr_hash(t$observed[match(split$analysis, t$row_id)])) ||
          !identical(d$training_coordinate_sha256,
            .msr_gdr_hash(unname(as.matrix(xy[match(split$analysis, xy$row_id), c("x", "y")])))) ||
          !is.finite(d$logLik) || d$nobs != length(split$analysis) ||
          any(!is.finite(c(d$K, d$native_K, d$ecological_K))) ||
          any(c(d$K, d$native_K, d$ecological_K) < 0) ||
          any(c(d$K, d$native_K, d$ecological_K) != floor(c(d$K, d$native_K, d$ecological_K))) ||
          d$K != d$native_K + d$ecological_K)
        stop("Field CV training selection records disagree with analysis targets or parameter counts.", call. = FALSE)
    }
    used <- sort(unique(unlist(lapply(x$folds, function(s) c(s$analysis, s$assessment)), use.names = FALSE)))
    folds <- lapply(x$folds[sort(names(x$folds))], function(s) lapply(s[c("analysis", "assessment")], sort))
    attributes(folds) <- list(names = names(folds))
    targets <- t[match(used, t$row_id), , drop = FALSE]; rownames(targets) <- NULL
    coordinates <- xy[match(used, xy$row_id), , drop = FALSE]; rownames(coordinates) <- NULL
    list(folds = folds, targets = targets, coordinates = coordinates, crs = x$settings$crs)
  })
  if (!all(vapply(signatures, identical, logical(1), signatures[[1L]])))
    stop("Field CV comparison requires identical splits, observations and projected coordinates.", call. = FALSE)
  invisible(signatures[[1L]])
}

.msr_field_cv_tab <- function(mod_list, metric, by_fold) {
  .msr_field_cv_comparison(mod_list)
  if (!metric %in% c("rmse", "mae", "neg_log_score") ||
      any(!vapply(mod_list, function(x) metric %in% x$settings$metrics, logical(1))))
    stop("Every field candidate must record the same requested metric.", call. = FALSE)
  groups <- if (by_fold) names(mod_list[[1L]]$folds) else "pooled"
  out <- lapply(groups, function(group) {
    fold <- if (by_fold) group else NULL
    rows <- lapply(names(mod_list), function(name) {
      x <- mod_list[[name]]; score <- .msr_field_cv_metric(x, metric, fold)
      data.frame(Modnames = name, Score = score$estimate, Metric = metric, Unit = score$unit,
        N_assessment = score$n, N_scored = score$n, Weight_sum = score$n,
        Warnings = length(unlist(if (by_fold) x$warnings[fold] else x$warnings, use.names = FALSE)))
    })
    tab <- do.call(rbind, rows); tab$Rank <- rank(tab$Score, ties.method = "min")
    tab$Delta <- tab$Score - min(tab$Score)
    if (by_fold) tab$Fold <- group
    tab[order(tab$Rank), c("Modnames", "Rank", "Score", "Delta", "Metric", "Unit",
      "N_assessment", "N_scored", "Weight_sum", "Warnings", if (by_fold) "Fold"), drop = FALSE]
  })
  result <- do.call(rbind, out); rownames(result) <- NULL; result
}

.msr_field_select_cv <- function(mod_list, criterion, ties, keep_fits) {
  .msr_field_cv_comparison(mod_list)
  first <- mod_list[[1L]]
  if (!all(vapply(mod_list, function(x) setequal(x$settings$metrics, first$settings$metrics), logical(1))))
    stop("Field CV selection candidates must retain matching metrics.", call. = FALSE)
  if (keep_fits && any(vapply(mod_list, function(x) is.null(x$fits), logical(1))))
    stop("Selected native fits require keep_fits = TRUE in every candidate.", call. = FALSE)
  out <- list(predictions = list(), joint_scores = list(), fold_results = list(), training = list(),
    selection_records = list(), diagnostics = list(), warnings = list(), errors = list(),
    prediction_contracts = list(), fits = if (keep_fits) list() else NULL)
  selection <- list(criterion = criterion, ties = ties, candidates = names(mod_list),
    chosen = stats::setNames(rep(NA_character_, length(first$folds)), names(first$folds)), tables = list())
  selected_settings <- list()
  for (name in names(first$folds)) {
    winner <- NULL
    error <- tryCatch({
      records <- lapply(mod_list, function(x) x$selection_records[[name]])
      k <- vapply(records, `[[`, 0, "K"); ll <- vapply(records, `[[`, 0, "logLik")
      n <- vapply(records, `[[`, 0, "nobs")
      if (criterion == "AICc" && any(n <= k + 1)) stop("Training field AICc requires n > K + 1.", call. = FALSE)
      value <- if (criterion == "BIC") -2 * ll + log(n) * k else -2 * ll + 2 * k +
        if (criterion == "AICc") 2 * k * (k + 1)/(n - k - 1) else 0
      if (any(!is.finite(value))) stop("Training field information criteria must be finite.", call. = FALSE)
      tab <- data.frame(Modnames = names(mod_list), K = k, logLik = ll, criterion_value = value)
      names(tab)[4L] <- criterion; selection$tables[[name]] <- tab
      best <- names(mod_list)[value == min(value)]
      if (length(best) > 1L && ties == "error") stop("Exact training-criterion tie; predeclare ties = 'first'.", call. = FALSE)
      winner <- best[1L]; NULL
    }, error = function(e) list(stage = "training_selection", message = conditionMessage(e)))
    x <- if (is.null(error)) mod_list[[winner]] else first
    pred <- x$predictions[x$predictions$fold == name, , drop = FALSE]
    if (!is.null(error)) {
      pred$prediction <- pred$squared_error <- pred$absolute_error <- NA_real_; pred$eligible <- FALSE
    } else {
      selection$chosen[name] <- winner
      for (item in c("joint_scores", "training", "selection_records", "diagnostics", "prediction_contracts"))
        out[[item]][name] <- x[[item]][name]
      if (keep_fits) out$fits[name] <- x$fits[name]
      selected_settings[name] <- list(x$settings)
    }
    out$predictions[[name]] <- pred
    out$fold_results[[name]] <- data.frame(fold = name, n_analysis = length(first$folds[[name]]$analysis),
      n_assessment = length(first$folds[[name]]$assessment), status = if (is.null(error)) "scored" else "failed")
    out$warnings[name] <- x$warnings[name]; out$errors[name] <- list(error)
  }
  out$selection <- selection; out$targets <- first$targets; out$coordinates <- first$coordinates
  out$folds <- first$folds; out$settings <- first$settings
  out$settings$backend <- "training_selection"; out$settings$covariance <- "per_fold_training_selection"
  out$settings$fold_settings <- selected_settings; out$settings$selection_uncertainty <- "not_integrated"
  .msr_field_cv_finish(out, first$targets$row_id)
}
