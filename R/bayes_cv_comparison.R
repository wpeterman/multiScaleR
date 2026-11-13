# Freeze score, split and observation records without retaining posterior fits.
.msr_bayes_cv_payload <- function(x) unclass(x)[c("predictions", "joint_scores", "fold_results",
  "training", "diagnostics", "warnings", "errors", "prediction_contracts",
  "targets", "folds", "settings")]

.msr_bayes_cv_metric <- function(x, fold = NULL) {
  p <- x$predictions; j <- x$joint_scores
  if (!is.null(fold)) {
    p <- p[p$fold == fold, , drop = FALSE]
    j <- j[j$fold == fold, , drop = FALSE]
  }
  n <- nrow(p)
  complete <- all(p$eligible)
  data.frame(metric = "neg_log_score",
    estimate = if (complete) sum(j$neg_log_score) / n else NA_real_,
    unit = "joint_block_density_per_observation",
    n_assessment = n, n_scored = sum(p$eligible), weight_sum = sum(p$eligible),
    complete = complete)
}

.msr_bayes_cv_finish <- function(x, ids) {
  for (name in c("predictions", "joint_scores", "fold_results")) {
    x[[name]] <- do.call(rbind, x[[name]])
    rownames(x[[name]]) <- NULL
  }
  x$complete <- all(x$fold_results$status == "scored")
  x$assessment_counts <- stats::setNames(tabulate(match(x$predictions$row_id, ids), length(ids)), ids)
  x$summary <- .msr_bayes_cv_metric(x)
  x$metrics <- do.call(rbind, lapply(names(x$folds), function(f)
    cbind(fold = f, .msr_bayes_cv_metric(x, f))))
  x$integrity_sha256 <- digest::digest(.msr_bayes_cv_payload(x), algo = "sha256")
  class(x) <- c("multiScaleR_bayes_cv", "multiScaleR_cv")
  x
}

.msr_bayes_cv_tab <- function(mod_list, metric, by_fold) {
  if (is.null(names(mod_list)) || anyNA(names(mod_list)) ||
      any(!nzchar(names(mod_list))) || anyDuplicated(names(mod_list)))
    stop("Supply a uniquely named list of Bayesian CV candidates.", call. = FALSE)
  if (!identical(metric, "neg_log_score"))
    stop("Bayesian CV comparison uses neg_log_score for joint predictive blocks.", call. = FALSE)
  signatures <- lapply(mod_list, function(x) {
    if (!inherits(x, "multiScaleR_bayes_cv") || !isTRUE(x$complete) ||
        any(x$fold_results$status != "scored"))
      stop("Compare complete Bayesian CV results; failed folds and other scoring targets cannot be mixed or dropped.", call. = FALSE)
    if (!identical(x$integrity_sha256, digest::digest(.msr_bayes_cv_payload(x), algo = "sha256")))
      stop("Recorded Bayesian CV targets, scores or diagnostics changed.", call. = FALSE)
    folds <- lapply(x$folds[sort(names(x$folds))], function(f)
      lapply(f[c("analysis", "assessment")], sort))
    settings <- x$settings[c("family", "link", "prediction_target", "new_scale_levels",
      "new_response_levels")]
    # Poisson and NB2 both score normalized mass on nonnegative counts.
    if (settings$family %in% c("poisson", "nbinom2")) settings$family <- "count"
    list(folds = folds, targets = x$targets, settings = settings)
  })
  if (!all(vapply(signatures, identical, logical(1), signatures[[1L]])))
    stop("Bayesian CV requires identical splits, observation targets, response measures, joint blocks and new-level policies.", call. = FALSE)
  groups <- if (by_fold) names(mod_list[[1L]]$folds) else "pooled"
  out <- lapply(groups, function(group) {
    fold <- if (by_fold) group else NULL
    rows <- lapply(names(mod_list), function(name) {
      x <- mod_list[[name]]; score <- .msr_bayes_cv_metric(x, fold)
      warnings <- if (by_fold) x$warnings[fold] else x$warnings
      data.frame(Modnames = name, Score = score$estimate, Metric = metric, Unit = score$unit,
        N_assessment = score$n_assessment, N_scored = score$n_scored,
        Weight_sum = score$weight_sum, Warnings = length(unlist(warnings, use.names = FALSE)))
    })
    tab <- do.call(rbind, rows)
    tab$Rank <- rank(tab$Score, ties.method = "min")
    best <- min(tab$Score)
    tab$Delta <- if (is.finite(best)) tab$Score - best else NA_real_
    if (by_fold) tab$Fold <- fold
    tab[order(tab$Rank), c("Modnames", "Rank", "Score", "Delta", "Metric", "Unit",
      "N_assessment", "N_scored", "Weight_sum", "Warnings", if (by_fold) "Fold")]
  })
  result <- do.call(rbind, out)
  rownames(result) <- NULL
  result
}
