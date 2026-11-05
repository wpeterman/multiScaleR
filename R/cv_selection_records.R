# Capture from actual native fits before releasing their large working objects.
# Formula text has no enclosing environment; records contain no model or raster.
.msr_cv_selection_controls <- function(settings) {
  settings[intersect(c("backend", "optim_args", "glm_control", "cox_control", "method", "distribution_control"), names(settings))]
}

.msr_cv_selection_record <- function(fit, diagnostics, controls) {
  if (!is.null(controls$distribution_control)) return(.msr_cv_distribution_record(fit, diagnostics, controls))
  messages <- character()
  record <- list(version = 1L, controls = controls, diagnostics = diagnostics)
  record$data <- tryCatch(withCallingHandlers({
    model <- if (inherits(fit, "multiScaleR")) fit$opt_mod else fit
    if (!inherits(model, "glm") || inherits(model, "gam") || !is.data.frame(model$model))
      stop("Selection records require a retained native canonical GLM.", call. = FALSE)
    dat <- data.frame(row.names = rownames(model$model))
    dat$outcome <- stats::model.response(model$model)
    y <- .msr_cv_outcome(outcome ~ 1, dat, model$family$family)
    w <- stats::model.weights(model$model)
    if (is.null(w)) w <- rep(1, nrow(dat))
    native_k <- .msr_parameter_count(model)
    list(logLik = .msr_model_loglik(model), K = .msr_selection_parameter_count(fit, native_k),
      native_K = native_k, ecological_K = .msr_selection_parameter_count(fit, 0L),
      nobs = .msr_selection_nobs(model), observation_ids = .msr_observation_ids(model),
      targets = data.frame(row_id = rownames(dat), observed = y$observed,
        trials = y$trials, weight = as.numeric(w), stringsAsFactors = FALSE),
      response_layout = if (is.matrix(dat$outcome)) "counts" else "single",
      signature = .msr_selection_signature(model),
      formula = paste(deparse(stats::formula(model)), collapse = " "),
      family = model$family$family, link = model$family$link)
  }, warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) { record$error <<- conditionMessage(e); NULL })
  record$warnings <- messages
  record
}

# Bind every analysis row, including zero-weight rows, before checking the
# effective likelihood signature. Count input keeps case frequencies separate.
.msr_cv_record_identity <- function(record, x, name, targets) {
  if (inherits(x, "multiScaleR_choice_cv")) return(.msr_choice_record_identity(record, x, name, targets))
  settings <- .msr_cv_fold_settings(x, name)
  if (.msr_cv_distribution_family(settings$family)) return(.msr_cv_distribution_record_identity(record, x, name, targets))
  x$settings <- settings
  if (!is.list(record) || !identical(record$version, 1L) ||
      !identical(record$controls, .msr_cv_selection_controls(x$settings)) ||
      !identical(record$diagnostics, x$diagnostics[[name]]) ||
      !is.character(record$warnings))
    stop("Invalid selection record controls, diagnostics or schema; rerun CV.", call. = FALSE)
  if (!is.null(record$error)) {
    if (!is.character(record$error) || length(record$error) != 1L || is.na(record$error) || !is.null(record$data))
      stop("Malformed selection metadata failure.", call. = FALSE)
    return(invisible(TRUE))
  }
  d <- record$data; family <- x$settings$family
  if (!is.list(d) || !identical(d$family, family) || !identical(d$link, x$settings$link))
    stop("Selection record family and link must match the CV scoring target.", call. = FALSE)
  expected <- targets[match(x$folds[[name]]$analysis, targets$row_id), , drop = FALSE]
  rownames(expected) <- NULL
  if (!is.data.frame(d$targets) || !identical(d$targets, expected))
    stop("Selection record training responses, trials, weights or analysis row IDs disagree with targets.", call. = FALSE)
  scalar <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)
  if (!all(c("K", "native_K", "ecological_K", "nobs") %in% names(d)) ||
      !scalar(d$logLik) || !all(vapply(d[c("K", "native_K", "ecological_K", "nobs")],
      function(v) scalar(v) && v >= 0 && v == floor(v), logical(1))) ||
      d$K != d$native_K + d$ecological_K || !is.character(d$formula) || length(d$formula) != 1L ||
      is.na(d$formula) || !nzchar(d$formula) || !is.character(d$response_layout) ||
      length(d$response_layout) != 1L || is.na(d$response_layout) || !d$response_layout %in% c("counts", "single"))
    stop("Invalid native likelihood, parameter counts or formula in selection record.", call. = FALSE)
  used <- expected$weight > 0 & (family != "binomial" | expected$trials > 0)
  if (!identical(d$observation_ids, expected$row_id[used]) || d$nobs != sum(used))
    stop("Selection record effective likelihood row IDs or native observation count disagree.", call. = FALSE)
  e <- expected[used, , drop = FALSE]; e <- e[order(e$row_id), , drop = FALSE]
  counts <- identical(d$response_layout, "counts")
  if ((family != "binomial" && counts) || (family == "binomial" && !counts && any(expected$trials != 1)))
    stop("Selection record response layout disagrees with trial definitions.", call. = FALSE)
  signature <- list(response = if (family == "binomial" && !counts) e$observed * e$weight else e$observed,
    weights = if (family == "binomial" && !counts) rep(1, nrow(e)) else e$weight,
    measure = if (family == "gaussian") "continuous" else "discrete",
    trials = if (family != "binomial") NULL else if (counts) e$trials else e$weight,
    reml = FALSE, design = NULL, partial_method = NULL, strata = NULL)
  if (!isTRUE(all.equal(d$signature, signature, tolerance = 0, check.attributes = FALSE)))
    stop("Selection record native likelihood signature disagrees with training targets.", call. = FALSE)
  invisible(TRUE)
}

.msr_cv_record_table <- function(records, criterion, names) {
  for (record in records) {
    for (message in record$warnings) warning(message, call. = FALSE)
    if (!is.null(record$error)) stop(record$error, call. = FALSE)
  }
  data <- lapply(records, `[[`, "data")
  n <- vapply(data, `[[`, integer(1), "nobs")
  ids <- lapply(data, function(d) sort(d$observation_ids))
  if (length(unique(n)) != 1L || !all(vapply(ids, identical, logical(1), ids[[1L]])))
    stop("Training records must describe the same effective observation sets.", call. = FALSE)
  .msr_assert_selection_signatures(lapply(data, `[[`, "signature"))
  k <- vapply(data, `[[`, numeric(1), "K")
  ll <- vapply(data, `[[`, numeric(1), "logLik")
  if (criterion == "AICc" && any(n <= k + 1L))
    stop("AICc requires n > K + 1 for every model; use `AICc = FALSE` for AIC.", call. = FALSE)
  if (criterion == "BIC")
    AICcmodavg::bictabCustom(logL = ll, K = k, modnames = names, nobs = n, sort = TRUE) else
    AICcmodavg::aictabCustom(logL = ll, K = k, modnames = names, second.ord = criterion == "AICc", nobs = n, sort = TRUE)
}
