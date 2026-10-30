# Prediction reuses the fitted model basis and standardization. It must not
# call the fitting preparer, which would estimate a new reference transform.
.msr_bayes_design_info <- function(frame, X, response = FALSE) {
  if (response) frame <- frame[-1L]
  levels <- lapply(frame, function(x) if (is.factor(x)) levels(x) else
    if (is.character(x)) levels(as.factor(x)) else NULL)
  levels <- levels[!vapply(levels, is.null, TRUE)]
  contrasts <- attr(X, "contrasts")
  # Realize the standard deterministic contrasts now, so later changes to
  # options(contrasts) cannot change the coefficient basis of a saved fit.
  standard <- c("contr.treatment", "contr.sum", "contr.helmert", "contr.poly", "contr.SAS")
  for (nm in names(contrasts)) {
    value <- contrasts[[nm]]
    if (is.character(value) && length(value) == 1L && value %in% standard)
      contrasts[[nm]] <- getExportedValue("stats", value)(levels[[nm]])
  }
  list(xlevels = levels, contrasts = contrasts)
}

.msr_bayes_prediction_matrix <- function(model, data, training, response = FALSE) {
  info <- model$prediction
  if (is.null(info)) {
    # Older compact fits have no prediction metadata. Accept reconstruction
    # only when it reproduces their actual retained training design.
    original <- stats::model.matrix(model$terms, model$frame)
    keep <- colnames(training)
    same <- all(keep %in% colnames(original)) &&
      isTRUE(all.equal(as.numeric(original[, keep, drop = FALSE]),
        as.numeric(training), tolerance = 1e-12))
    if (!same)
      stop("Cannot reconstruct this saved fit's contrast coding. Use its original contrast settings or refit with prediction metadata.", call. = FALSE)
    info <- .msr_bayes_design_info(model$frame, original, response)
  }
  tt <- if (response) stats::delete.response(model$terms) else model$terms
  # Fitting strips caller environments to avoid serializing whole workspaces.
  # Base/stats transforms and namespace-qualified functions remain available.
  attr(tt, ".Environment") <- asNamespace("stats")
  frame <- stats::model.frame(tt, data, na.action = stats::na.fail,
    xlev = info$xlevels, drop.unused.levels = FALSE)
  X <- stats::model.matrix(tt, frame, contrasts.arg = info$contrasts)
  keep <- colnames(training)
  if (!all(keep %in% colnames(X)) || any(!is.finite(X)))
    stop("New fixed-effect design does not match the fitted model or contains nonfinite values.", call. = FALSE)
  list(X = X[, keep, drop = FALSE], frame = frame)
}

.msr_bayes_empty_new_groups <- function() data.frame(component = character(),
  covariate = character(), group = character(), level = character(), action = character())

.msr_bayes_prediction_group <- function(term, data, policy, start, component, covariate = "") {
  g <- term$group
  value <- data[[g]]
  if (is.null(value) || !is.atomic(value) || !is.null(dim(value)) ||
      anyNA(value) || any(!nzchar(as.character(value))))
    stop("New data need nonmissing group IDs in '", g, "' for included ", component,
      " effects.", call. = FALSE)
  value <- as.character(value)
  added <- sort(setdiff(unique(value), term$levels))
  if (length(added) && policy == "error")
    stop("New ", component, " group levels in '", g, "': ",
      paste(utils::head(added, 6L), collapse = ", "),
      ". Choose the corresponding new-level policy explicitly.", call. = FALSE)
  term$levels <- c(term$levels, added)
  term$index <- match(value, term$levels)
  old_n <- length(term$z_index)
  term$z_index <- c(term$z_index, start + seq_along(added))
  metadata <- if (!length(added)) .msr_bayes_empty_new_groups() else
    data.frame(component = component, covariate = covariate, group = g,
      level = added, action = policy)
  list(term = term, n = length(added), new_rows = term$index > old_n, metadata = metadata)
}

.msr_bayes_prediction_scales <- function(spec, data, policy, include = TRUE) {
  if (.msr_has_gaussian_blocks(spec))
    stop("Gaussian projection blocks require the dedicated private spatial prediction adapter.", call. = FALSE)
  count <- 0L; metadata <- .msr_bayes_empty_new_groups(); flags <- list()
  for (nm in spec$covariates) {
    sm <- spec$scales[[nm]]
    sm$X <- .msr_bayes_prediction_matrix(sm, data, sm$X)$X
    new_rows <- rep(FALSE, nrow(data))
    if (include) for (k in seq_along(sm$random)) {
      mapped <- .msr_bayes_prediction_group(sm$random[[k]], data, policy,
        length(spec$parameters$scale_z) + count, "scale", nm)
      sm$random[[k]] <- mapped$term
      count <- count + mapped$n; new_rows <- new_rows | mapped$new_rows
      metadata <- rbind(metadata, mapped$metadata)
    } else sm$random <- list()
    variables <- union(all.vars(sm$fixed), if (include) names(sm$random) else character())
    labels <- if (!length(variables)) rep("shared", nrow(data)) else
      apply(as.data.frame(lapply(data[variables], as.character)), 1L,
        function(z) paste(paste(variables, z, sep = "="), collapse = ", "))
    representatives <- if (!length(variables)) 1L else which(!duplicated(data[variables]))
    sm$representatives <- representatives; sm$labels <- labels[representatives]
    spec$scales[[nm]] <- sm; flags[[nm]] <- new_rows
  }
  list(spec = spec, n_scale = count, scale_flags = flags, new_groups = metadata,
    scale_policy = policy)
}

.msr_bayes_prediction_re_form <- function(x, name) {
  if (!is.null(x) && !(length(x) == 1L && is.atomic(x) && is.na(x)))
    stop(name, " must be NULL or NA.", call. = FALSE)
  is.null(x)
}

.msr_bayes_prediction_prepare <- function(object, newdata, kernel_inputs, row_id,
    scale_policy, response_policy, include_scale, include_response, include_fixed = FALSE) {
  validate_scalar_logical(include_fixed, "include_fixed")
  spec <- object$specification
  if (.msr_has_gaussian_blocks(spec))
    stop("Gaussian projection blocks require the dedicated private spatial prediction adapter.", call. = FALSE)
  if (is.null(newdata)) {
    if (!is.null(kernel_inputs) || !is.null(row_id))
      stop("kernel_inputs and row_id require newdata.", call. = FALSE)
    if (!include_scale) for (nm in spec$covariates) spec$scales[[nm]]$random <- list()
    return(list(spec = spec, n_scale = 0L, n_response = 0L,
      scale_policy = scale_policy, response_policy = response_policy,
      scale_flags = lapply(spec$covariates, function(x) rep(FALSE, length(spec$y))),
      response_flags = rep(FALSE, length(spec$y)), new_groups = .msr_bayes_empty_new_groups()))
  }
  if (!is.data.frame(newdata) || nrow(newdata) < 1L || anyDuplicated(names(newdata)))
    stop("newdata must be a nonempty data frame with distinct column names.", call. = FALSE)
  has_surface <- length(spec$surface$buffers) > 0L
  if (is.null(kernel_inputs) || (is.null(kernel_inputs$binned) && !has_surface))
    stop("Supply matching binned kernel_inputs for newdata.", call. = FALSE)
  if (!identical(kernel_inputs$kernel, spec$kernel))
    stop("Prediction kernel_inputs must use the fitted kernel.", call. = FALSE)
  validate_scalar_numeric(kernel_inputs$unit_conv, "kernel_inputs$unit_conv", positive = TRUE)
  for (nm in c("min_D", "max_D"))
    validate_scalar_numeric(kernel_inputs[[nm]], paste0("kernel_inputs$", nm), positive = TRUE)
  if (kernel_inputs$min_D >= kernel_inputs$max_D)
    stop("Prediction min_D must be smaller than max_D.", call. = FALSE)
  if (is.null(row_id)) row_id <- spec$row_id
  keys <- if (is.null(row_id)) rownames(newdata) else {
    if (!is.character(row_id) || length(row_id) != 1L || !row_id %in% names(newdata))
      stop("row_id must name a prepared point-ID column in newdata.", call. = FALSE)
    newdata[[row_id]]
  }
  points <- .msr_scale_point_ids(kernel_inputs)
  surface_names <- names(spec$surface$buffers)
  mean_names <- setdiff(spec$covariates, surface_names)
  if (!all(mean_names %in% names(kernel_inputs$binned$vsum)))
    stop("Prediction bins must contain every fitted landscape predictor.", call. = FALSE)
  variables <- kernel_inputs$scale_vars
  if (!is.null(variables)) {
    used <- variables[variables$covariate %in% spec$covariates, , drop = FALSE]
    if (any(used$type != "kernel" & !used$covariate %in% surface_names))
      stop("Prediction inputs must supply kernel means for the fitted landscape predictors.", call. = FALSE)
  }
  if (has_surface) {
    expected <- spec$surface$specs
    if (is.null(variables)) stop("Prediction surface specifications are required.", call. = FALSE)
    actual <- variables[match(expected$covariate, variables$covariate), , drop = FALSE]
    fields <- c("covariate", "source", "type", "metric", "weighted", "optimize")
    if (!isTRUE(all.equal(actual[, fields], expected[, fields], check.attributes = FALSE)))
      stop("Prediction surface source and metric specifications must match the fitted model.", call. = FALSE)
  }
  if (!is.atomic(keys) || !is.null(dim(keys)) || length(keys) != nrow(newdata) ||
      anyNA(keys) || any(!nzchar(as.character(keys))))
    stop("Each new observation needs a nonmissing prepared point ID.", call. = FALSE)
  rows <- match(as.character(keys), as.character(points))
  if (anyNA(rows)) stop("New observation IDs do not match prediction binned point_ids; supply row_id explicitly.", call. = FALSE)
  if (include_fixed && length(spec$fixed_reference$covariates)) {
    reference <- spec$fixed_reference
    # Earlier fixed-feature specifications already retained this fitted bound.
    if (is.null(reference$max_D)) reference$max_D <- spec$support[["max_D"]]
    newdata[spec$fixed_reference$covariates] <- as.data.frame(
      .msr_fixed_reference_values(reference, kernel_inputs, rows))
  }
  if (has_surface) {
    prepared <- kernel_inputs
    ratio <- prepared$unit_conv / spec$unit_conv
    if (!is.null(prepared$binned)) prepared$binned$bin_rep <- prepared$binned$bin_rep * ratio
    prepared$d_list <- lapply(prepared$d_list, function(d) d * ratio)
    engine <- .msr_prepare_surface_engine(prepared, spec$covariates, rows, standardize = FALSE)
    spec$surface$buffers <- engine$buffers
    spec$surface$memory <- .msr_scale_surface_size(spec$surface)
    engine <- engine$means
    distance <- engine$distance
  } else {
    engine <- .msr_prepare_kernel_engine(kernel_inputs$binned, spec$covariates,
      rows, spec$kernel, standardize = FALSE)
    distance <- engine$distance * kernel_inputs$unit_conv / spec$unit_conv
  }
  spec$V <- engine$V; spec$C <- engine$C
  spec$delta <- lapply(mean_names, function(nm) {
    d <- matrix(rep(distance, each = nrow(spec$C[[nm]])), nrow(spec$C[[nm]]))
    nearest <- apply(replace(d, spec$C[[nm]] <= 0, Inf), 1L, min)
    pmax(if (spec$kernel == "gaussian") d^2 - nearest^2 else d - nearest, 0)
  }); names(spec$delta) <- mean_names
  model_data <- newdata
  for (nm in spec$covariates) model_data[[nm]] <- 0
  fixed <- .msr_bayes_prediction_matrix(spec$response, model_data, spec$X_fixed, TRUE)
  spec$X_fixed <- fixed$X
  spec$offset <- stats::model.offset(fixed$frame)
  if (is.null(spec$offset)) spec$offset <- rep(0, nrow(newdata))
  if (any(!is.finite(spec$offset))) stop("Prediction offsets must be finite.", call. = FALSE)
  context <- .msr_bayes_prediction_scales(spec, newdata, scale_policy, include_scale)
  spec <- context$spec; response_count <- 0L; response_flags <- rep(FALSE, nrow(newdata))
  if (include_response) for (k in seq_along(spec$response$random)) {
    mapped <- .msr_bayes_prediction_group(spec$response$random[[k]], newdata,
      response_policy, length(spec$parameters$response_z) + response_count, "response")
    spec$response$random[[k]] <- mapped$term
    response_count <- response_count + mapped$n
    response_flags <- response_flags | mapped$new_rows
    context$new_groups <- rbind(context$new_groups, mapped$metadata)
  } else spec$response$random <- list()
  spec$data <- newdata; spec$observation_ids <- as.character(keys); spec$rows <- rows
  spec$y <- rep(0, nrow(newdata)); spec$trials <- rep(1, nrow(newdata))
  context$spec <- spec; context$n_response <- response_count
  context$response_policy <- response_policy; context$response_flags <- response_flags
  context
}

.msr_bayes_prediction_parameters <- function(p, context) {
  for (component in c("scale", "response")) {
    n <- context[[paste0("n_", component)]]
    if (n > 0L) {
      name <- paste0(component, "_z")
      p[[name]] <- c(p[[name]], if (context[[paste0(component, "_policy")]] == "sample")
        stats::rnorm(n) else rep(0, n))
    }
  }
  p
}

.msr_bayes_prediction_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  validate_scalar_numeric(seed, "seed", integerish = TRUE, lower = 1, upper = .Machine$integer.max)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed); force(code)
}
