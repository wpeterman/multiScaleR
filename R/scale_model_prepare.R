# Freeze row identity, formula designs and reference standardization without
# choosing population priors or an inference target. The compatibility wrappers
# retain their own prior settings and numerical initialization policies.
.msr_scale_prepare_inputs <- function(formula, inputs, data, row_id) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("Supply a two-sided response formula.", call. = FALSE)
  if (!is.data.frame(data) || nrow(data) < 2L || anyDuplicated(names(data)))
    stop("data must be a data frame with at least two observations and distinct columns.",
         call. = FALSE)
  if (!inputs$kernel %in% c("gaussian", "exp"))
    stop("Bayesian scales require Gaussian or exponential kernel means or weighted SQ.", call. = FALSE)
  for (nm in c("unit_conv", "min_D", "max_D"))
    validate_scalar_numeric(inputs[[nm]], nm, positive = TRUE)
  if (inputs$min_D >= inputs$max_D)
    stop("min_D must be smaller than max_D.", call. = FALSE)
  ids <- .msr_scale_point_ids(inputs)
  keys <- if (is.null(row_id)) rownames(data) else {
    if (!is.character(row_id) || length(row_id) != 1L || !row_id %in% names(data))
      stop("row_id must name a column of prepared point IDs in data.", call. = FALSE)
    data[[row_id]]
  }
  if (length(keys) != nrow(data) || anyNA(keys))
    stop("Each observation needs a nonmissing prepared point ID.", call. = FALSE)
  rows <- match(as.character(keys), as.character(ids))
  if (anyNA(rows))
    stop("Observation IDs do not match binned point_ids; supply row_id explicitly.",
         call. = FALSE)
  surface <- .msr_scale_surface_names(inputs)
  covs <- intersect(all.vars(formula[[3L]]), c(names(inputs$binned$vsum), surface))
  if (!length(covs))
    stop("The response formula must include at least one binned kernel predictor.", call. = FALSE)
  if (length(setdiff(all.vars(formula), c(names(data), covs))))
    stop("All response formula variables must be columns in data or prepared kernel predictors.",
         call. = FALSE)
  specs <- inputs$scale_vars
  if (!is.null(specs)) {
    used <- specs[specs$covariate %in% all.vars(formula), , drop = FALSE]
    optimized <- used[used$optimize, , drop = FALSE]
    if (any(!(optimized$type == "kernel" | optimized$covariate %in% surface)) || any(!optimized$covariate %in% covs))
      stop("Only kernel means and optimized weighted SQ can be spatial predictors in the Bayesian response formula.",
           call. = FALSE)
  }
  list(rows = rows, keys = keys, covariates = covs)
}

.msr_scale_reference <- function(reference_scale, covariates) {
  if (!is.numeric(reference_scale) || any(!is.finite(reference_scale) | reference_scale <= 0))
    stop("reference_scale must contain positive distances in map units.", call. = FALSE)
  if (length(reference_scale) == 1L && is.null(names(reference_scale)))
    reference_scale <- rep(reference_scale, length(covariates))
  else if (!is.null(names(reference_scale)) && !anyDuplicated(names(reference_scale)) &&
           setequal(names(reference_scale), covariates)) reference_scale <- reference_scale[covariates]
  else stop("reference_scale must be a scalar or name every kernel predictor.", call. = FALSE)
  stats::setNames(reference_scale, covariates)
}

.msr_scale_prepare <- function(formula, inputs, data, scales, family, row_id,
                                reference_scale, prepared = NULL) {
  if (is.null(prepared)) prepared <- .msr_scale_prepare_inputs(formula, inputs, data, row_id)
  rows <- prepared$rows; keys <- prepared$keys; covs <- prepared$covariates
  reference_scale <- .msr_scale_reference(reference_scale, covs)
  engine <- if (any(covs %in% .msr_scale_surface_names(inputs)))
    .msr_prepare_surface_engine(inputs, covs, rows, standardize = FALSE) else
      .msr_prepare_kernel_engine(inputs$binned, covs, rows, inputs$kernel, standardize = FALSE)
  means <- if (inherits(engine, "msr_surface_engine")) engine$means else engine
  # Keep the same pooled reference transformation for every inference target.
  # Do not restandardize separately within years or at each fitted scale.
  if (inherits(engine, "msr_surface_engine")) {
    raw <- matrix(0, length(rows), length(covs), dimnames = list(NULL, covs))
    if (!is.null(means)) raw[, means$covariates] <- .msr_eval_kernel_engine(means,
      unname(log(reference_scale[means$covariates] / inputs$unit_conv)))$raw
    for (nm in names(engine$buffers)) raw[, nm] <- .msr_scale_surface_raw(
      engine$buffers[[nm]], rep(log(reference_scale[[nm]] / inputs$unit_conv), length(rows)),
      inputs$kernel)
  } else raw <- .msr_eval_kernel_engine(engine,
    unname(log(reference_scale / inputs$unit_conv)))$raw
  center <- colMeans(raw)
  spread <- apply(raw, 2L, stats::sd)
  if (any(!is.finite(spread) | spread <= 1e-12))
    stop("Kernel predictors must vary across fitting observations at the prior median scale.",
         call. = FALSE)
  model_data <- data
  for (j in seq_along(covs))
    model_data[[covs[j]]] <- (raw[, j] - center[j]) / spread[j]
  rf <- .msr_bayes_formula(formula, model_data, response = TRUE)
  tt <- rf$terms
  for (label in attr(tt, "term.labels")) {
    expr <- str2lang(label)
    hit <- intersect(all.vars(expr), covs)
    if (length(hit) && !(length(hit) == 1L && identical(expr, as.name(hit))))
      stop("Kernel predictors must be bare additive terms in the response formula.",
           call. = FALSE)
  }
  for (i in attr(tt, "offset"))
    if (length(intersect(all.vars(attr(tt, "variables")[[i + 1L]]), covs)))
      stop("Kernel predictors cannot occur in offsets.", call. = FALSE)
  X <- stats::model.matrix(tt, rf$frame)
  if (any(!is.finite(X)) || anyDuplicated(colnames(X)) ||
      qr(X)$rank < ncol(X) || !all(covs %in% colnames(X)))
    stop("The response formula needs a finite, full-rank fixed-effect design.", call. = FALSE)
  fixed_columns <- setdiff(colnames(X), covs)
  Xfixed <- X[, fixed_columns, drop = FALSE]
  rf$prediction <- .msr_bayes_design_info(rf$frame, X, response = TRUE)
  y <- stats::model.response(rf$frame)
  trials <- rep(1, nrow(data))
  if (family == "binomial") {
    if (is.matrix(y) && ncol(y) == 2L) {
      if (any(!is.finite(y) | y < 0 | y != floor(y)))
        stop("Binomial successes and failures must be nonnegative integers.", call. = FALSE)
      trials <- rowSums(y)
      y <- y[, 1L]
    } else if (is.factor(y) && nlevels(y) == 2L) y <- as.numeric(y) - 1
    if (any(trials < 1) || any(!y %in% 0:1) && all(trials == 1))
      stop("Use binary responses or cbind(successes, failures).", call. = FALSE)
  }
  if (!is.numeric(y) || !is.null(dim(y)) || length(y) != nrow(data) ||
      any(!is.finite(y)) || (family != "gaussian" && any(y < 0 | y != floor(y))))
    stop("The response must be a finite numeric vector appropriate to the family.",
         call. = FALSE)
  offset <- stats::model.offset(rf$frame)
  if (is.null(offset)) offset <- rep(0, length(y))
  if (any(!is.finite(offset))) stop("Offsets must be finite.", call. = FALSE)
  if (is.null(scales)) scales <- stats::setNames(rep(list(~1), length(covs)), covs)
  if (!is.list(scales) || is.null(names(scales)) || anyDuplicated(names(scales)) ||
      !setequal(names(scales), covs))
    stop("scales must be a named list with one formula for every kernel predictor.",
         call. = FALSE)
  scale_models <- vector("list", length(covs))
  names(scale_models) <- covs
  scale_beta <- scale_z <- log_scale_sd <- numeric()
  for (nm in covs) {
    if (length(intersect(all.vars(scales[[nm]]), covs)))
      stop("Scale formulas cannot use kernel predictors.", call. = FALSE)
    if (length(setdiff(all.vars(scales[[nm]]), names(data))))
      stop("All scale formula variables must be columns in data.", call. = FALSE)
    sf <- .msr_bayes_formula(scales[[nm]], data)
    if (length(attr(sf$terms, "offset")))
      stop("Scale formulas cannot contain offsets.", call. = FALSE)
    S <- stats::model.matrix(sf$terms, sf$frame)
    if (any(!is.finite(S)) || qr(S)$rank < ncol(S))
      stop("Scale formulas need finite, full-rank fixed-effect designs.", call. = FALSE)
    intercept <- attr(sf$terms, "intercept") == 1L
    cell_means <- !intercept && all(S %in% c(0, 1)) && all(rowSums(S) == 1)
    if (!intercept && !cell_means)
      stop("An intercept-free scale formula must define disjoint factor levels, such as ~0 + year.",
           call. = FALSE)
    sf$X <- S
    sf$prediction <- .msr_bayes_design_info(sf$frame, S)
    sf$beta_index <- length(scale_beta) + seq_len(ncol(S))
    # Identify distance locations separately from log-distance contrasts.
    # Prior attachment and likelihood starts both use this fixed parameter map.
    sf$distance_location <- if (cell_means) rep(TRUE, ncol(S)) else
      c(TRUE, rep(FALSE, ncol(S) - 1L))
    scale_beta <- c(scale_beta, ifelse(sf$distance_location,
      log(reference_scale[[nm]] / inputs$unit_conv), 0))
    for (k in seq_along(sf$random)) {
      term <- sf$random[[k]]
      term$z_index <- length(scale_z) + seq_along(term$levels)
      term$sd_index <- length(log_scale_sd) + 1L
      scale_z <- c(scale_z, rep(0, length(term$levels)))
      log_scale_sd <- c(log_scale_sd, log(.2))
      sf$random[[k]] <- term
    }
    variables <- all.vars(scales[[nm]])
    labels <- if (!length(variables)) rep("shared", nrow(data)) else
      apply(as.data.frame(lapply(data[variables], as.character)), 1L,
        function(z) paste(paste(variables, z, sep = "="), collapse = ", "))
    sf$representatives <- which(!duplicated(labels))
    sf$labels <- labels[sf$representatives]
    scale_models[[nm]] <- sf
  }
  response_z <- log_response_sd <- numeric()
  for (k in seq_along(rf$random)) {
    term <- rf$random[[k]]
    term$z_index <- length(response_z) + seq_along(term$levels)
    term$sd_index <- length(log_response_sd) + 1L
    response_z <- c(response_z, rep(0, length(term$levels)))
    log_response_sd <- c(log_response_sd, log(.2))
    rf$random[[k]] <- term
  }
  # Shift exponents by each row's nearest supported bin. The constant cancels
  # from the ratio and prevents all weights underflowing at small scales.
  mean_names <- if (is.null(means)) character() else means$covariates
  delta <- lapply(mean_names, function(nm) {
    C <- means$C[[nm]]
    d <- matrix(rep(means$distance, each = nrow(C)), nrow(C))
    nearest <- apply(replace(d, C <= 0, Inf), 1L, min)
    ans <- if (inputs$kernel == "gaussian") d^2 - nearest^2 else d - nearest
    pmax(ans, 0)
  })
  names(delta) <- mean_names
  # All evaluated designs are already numeric. Keeping a caller's formula
  # environment could otherwise serialize unrelated rasters or an entire fit.
  freeze_formula <- function(x) {
    environment(x$formula) <- baseenv()
    environment(x$fixed) <- baseenv()
    attr(x$terms, ".Environment") <- baseenv()
    attr(x$frame, "terms") <- x$terms
    x
  }
  rf <- freeze_formula(rf)
  scale_models <- lapply(scale_models, freeze_formula)
  out <- list(data = data, rows = rows, observation_ids = as.character(keys), row_id = row_id,
       covariates = covs, kernel = inputs$kernel, unit_conv = inputs$unit_conv,
       support = c(min_D = inputs$min_D, max_D = inputs$max_D),
       V = means$V, C = means$C, delta = delta, center = center, spread = spread,
       response = rf, X_fixed = Xfixed, beta_names = c(fixed_columns, covs),
       scales = scale_models, family = family, y = as.numeric(y), trials = as.numeric(trials),
       offset = offset,
       parameters = list(beta = rep(0, ncol(Xfixed) + length(covs)),
         scale_beta = scale_beta, scale_z = scale_z, log_scale_sd = log_scale_sd,
         response_z = response_z, log_response_sd = log_response_sd,
         log_dispersion = if (family %in% c("gaussian", "nbinom2"))
           0 else numeric()))
  if (inherits(engine, "msr_surface_engine")) {
    out$surface <- list(buffers = engine$buffers, specs = engine$surface_specs)
    out$surface$memory <- .msr_scale_surface_size(out$surface)
  }
  fixed_variables <- unique(c(all.vars(formula[[3L]]),
    unlist(lapply(scales, all.vars), use.names = FALSE)))
  fixed_reference <- .msr_fixed_reference_contract(inputs, data, rows, fixed_variables)
  if (length(fixed_reference$covariates)) out$fixed_reference <- fixed_reference
  out
}
