.msr_choice_record <- function(fit, diagnostics, controls, split, targets) {
  record <- list(version = 2L, controls = controls, diagnostics = diagnostics); messages <- character()
  record$data <- tryCatch(withCallingHandlers({
    .msr_choice_identity(fit, split, targets, controls$method)
    m <- .msr_choice_native_frame(fit)$model
    native_k <- .msr_parameter_count(m)
    t <- targets[match(split$analysis, targets$row_id), , drop = FALSE]; rownames(t) <- NULL
    list(logLik = .msr_model_loglik(m), K = .msr_selection_parameter_count(fit, native_k),
      native_K = native_k, ecological_K = .msr_selection_parameter_count(fit, 0L),
      nobs = .msr_selection_nobs(m), observation_ids = rownames(m$model), targets = t,
      signature = .msr_selection_signature(m), formula = paste(deparse(stats::formula(m)), collapse = " "),
      method = m$method, family = "conditional_choice", link = "softmax")
  }, warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) { record$error <<- conditionMessage(e); NULL })
  record$warnings <- messages; record
}

.msr_choice_record_identity <- function(record, x, name, targets) {
  if (!is.list(record) || !identical(record$version, 2L) ||
      !identical(record$controls, .msr_cv_selection_controls(x$settings)) ||
      !identical(record$diagnostics, x$diagnostics[[name]]) || !is.character(record$warnings))
    stop("Invalid choice selection record controls, diagnostics or schema.", call. = FALSE)
  if (!is.null(record$error)) {
    if (!is.character(record$error) || length(record$error) != 1L || is.na(record$error) || !is.null(record$data))
      stop("Malformed choice selection metadata failure.", call. = FALSE)
    return(invisible(TRUE))
  }
  d <- record$data
  expected <- targets[match(x$folds[[name]]$analysis, targets$row_id), , drop = FALSE]; rownames(expected) <- NULL
  scalar <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)
  if (!is.list(d) || !identical(d$family, "conditional_choice") || !identical(d$link, "softmax") ||
      !identical(d$method, x$settings$method) || !identical(d$targets, expected) ||
      !identical(d$observation_ids, expected$row_id))
    stop("Choice record training targets, row IDs or partial method disagree.", call. = FALSE)
  if (!all(c("K", "native_K", "ecological_K", "nobs") %in% names(d)) || !scalar(d$logLik) ||
      !all(vapply(d[c("K", "native_K", "ecological_K", "nobs")],
        function(v) scalar(v) && v >= 0 && v == floor(v), logical(1))) ||
      d$K != d$native_K + d$ecological_K || d$nobs != length(unique(expected$choice)) ||
      !is.character(d$formula) || length(d$formula) != 1L || is.na(d$formula) || !nzchar(d$formula))
    stop("Invalid choice native likelihood, event count, parameter count or formula.", call. = FALSE)
  t <- expected[order(expected$row_id), , drop = FALSE]
  signature <- list(response = unname(survival::Surv(rep(1, nrow(t)), t$observed)), weights = rep(1, nrow(t)),
    measure = "partial_likelihood", trials = NULL, reml = FALSE, design = NULL,
    partial_method = d$method, strata = match(t$choice, unique(t$choice)))
  agreement <- all.equal(d$signature, signature, tolerance = 0, check.attributes = FALSE)
  if (!isTRUE(agreement))
    stop("Choice record partial-likelihood signature disagrees with training targets: ", paste(agreement, collapse = "; "), call. = FALSE)
  invisible(TRUE)
}

# No fallback from alternative-level binary records: set membership is required.
.msr_choice_comparison <- function(mod_list, metric = NULL) {
  signatures <- lapply(mod_list, function(x) {
    s <- x$settings
    if (!inherits(x, "multiScaleR_choice_cv") || !isTRUE(x$complete) ||
        !is.data.frame(x$fold_results) || any(x$fold_results$status != "scored") ||
        !identical(s$family, "conditional_choice") || !identical(s$link, "softmax") ||
        !identical(s$prediction_target, "within_set_choice_probability") ||
        !identical(s$parameter_uncertainty, "not_integrated") ||
        !identical(s$observation_intervals, "not_provided") || !identical(s$metrics, "neg_log_score") ||
        !s$method %in% c("exact", "efron", "breslow") || (!is.null(metric) && !identical(metric, "neg_log_score")))
      stop("Compare complete qualified conditional-choice CV targets using neg_log_score only; failed folds cannot be dropped.", call. = FALSE)
    used <- sort(unique(unlist(lapply(x$folds, function(f) c(f$analysis, f$assessment)), use.names = FALSE)))
    t <- x$targets
    columns <- c("row_id", "observed", "trials", "weight", "choice")
    if (!is.data.frame(t) || !all(columns %in% names(t)) || anyNA(t$row_id) || anyDuplicated(t$row_id) ||
        !setequal(t$row_id, used) || !is.numeric(t$observed) || anyNA(t$observed) || any(!t$observed %in% c(0, 1)) ||
        !is.character(t$choice) || anyNA(t$choice) || any(!nzchar(t$choice)) ||
        !is.numeric(t$trials) || anyNA(t$trials) || any(t$trials != 1) ||
        !is.numeric(t$weight) || anyNA(t$weight) || any(t$weight != 1))
      stop("Invalid choice targets or unsupported trial/frequency definitions.", call. = FALSE)
    t <- t[match(used, t$row_id), columns, drop = FALSE]; rownames(t) <- NULL
    .msr_clogit_groups(t$choice, t$observed); .msr_choice_splits(x$folds, t)
    folds <- lapply(x$folds[sort(names(x$folds))], function(f) lapply(f[c("analysis", "assessment")], sort))
    attributes(folds) <- list(names = names(folds))
    p <- x$predictions
    if (!is.data.frame(p) || !all(c(columns, "fold", "linear_predictor", "prediction", "log_probability", "eligible", "neg_log_score") %in% names(p)))
      stop("Missing recorded choice assessment targets or scores.", call. = FALSE)
    if (anyNA(p$fold) || any(!p$fold %in% names(folds))) stop("Unknown recorded choice folds.", call. = FALSE)
    for (name in names(folds)) {
      ids <- p$row_id[p$fold == name]
      if (anyNA(ids) || anyDuplicated(ids) || !setequal(ids, folds[[name]]$assessment))
        stop("Recorded choice alternatives disagree with assessment splits.", call. = FALSE)
    }
    expected <- t[match(p$row_id, t$row_id), columns, drop = FALSE]
    actual <- p[columns]; rownames(actual) <- rownames(expected) <- NULL
    if (!identical(actual, expected)) stop("Recorded choice responses or set membership disagree with targets.", call. = FALSE)
    check <- .msr_choice_score(p, "neg_log_score")
    for (field in c("prediction", "log_probability", "eligible", "neg_log_score"))
      if (!isTRUE(all.equal(p[[field]], check[[field]], tolerance = 0, check.attributes = FALSE)))
        stop("Recorded choice probabilities or losses disagree with linear predictors.", call. = FALSE)
    list(folds = folds, targets = t, family = s$family, link = s$link, method = s$method)
  })
  first <- signatures[[1L]]
  if (!all(vapply(signatures, identical, logical(1), first)))
    stop("Choice CV requires identical complete sets, alternatives, cases, splits and native partial methods.", call. = FALSE)
  invisible(first)
}
