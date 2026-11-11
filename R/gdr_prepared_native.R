# Bind the fitted frame's physical survey keys to its native design rows.
.msr_gdr_prepared_rows <- function(frame, event_keys, event_id, interval_units) {
  if (!all(event_id %in% union(names(frame@siteCovs), names(frame@yearlySiteCovs))))
    stop("Survey keys must be retained in the original site or primary covariates.", call. = FALSE)
  audit <- validate_unmarked_data(frame, event_keys, event_id, interval_units = interval_units)
  if (!isTRUE(audit$valid) || anyNA(frame@yDistance) || anyNA(frame@yRemoval))
    stop("Prepared GDR requires verified survey order and complete paired histories.", call. = FALSE)
  n <- nrow(frame@yDistance)
  primary <- frame@numPrimary
  intervals <- ncol(frame@yRemoval) / primary
  event <- seq_len(nrow(event_keys))
  list(lambda = data.frame(site_row = seq_len(n)),
    phi = data.frame(event = event, site_row = rep(seq_len(n), each = primary),
      primary = rep(seq_len(primary), n)),
    dist = data.frame(event = event),
    rem = data.frame(event = rep(event, each = intervals), interval = rep(seq_len(intervals), length(event))),
    protocol = audit$protocol, survey_findings = audit$findings)
}

# Use observed membership to verify the native column order. Display labels
# alone cannot identify the latent coordinate of a group.
.msr_gdr_prepared_groups <- function(formula, raw, Z) {
  labels <- attr(stats::terms(formula, data = raw), "term.labels")
  terms <- Filter(.msr_marginal_random_call, lapply(labels, str2lang))
  if (!length(terms) && !ncol(Z)) return(NULL)
  term <- if (length(terms) == 1L) terms[[1L]] else NULL
  while (is.call(term) && identical(term[[1L]], as.name("("))) term <- term[[2L]]
  if (is.null(term) || !identical(term[[1L]], as.name("|")) ||
      !identical(term[[2L]], 1) || !is.symbol(term[[3L]]))
    stop("Prepared GDR group identity requires one bare random intercept per component.", call. = FALSE)
  group <- as.character(term[[3L]])
  observed <- raw[[group]]
  if ((!is.factor(observed) && !is.character(observed) && !is.numeric(observed)) || !is.null(dim(observed)))
    stop("Prepared GDR groups require original IDs in a single column.", call. = FALSE)
  ids <- colnames(Z)
  coordinate <- .msr_gdr_field_map(Z, nrow(raw))
  levels <- if (is.factor(observed)) levels(observed) else unique(as.character(observed))
  if (nrow(Z) != nrow(raw) || !.msr_valid_observation_ids(ids, ncol(Z)) ||
      !all(ids %in% levels) || anyNA(observed) || !identical(unname(ids[coordinate]), unname(as.character(observed))))
    stop("Original group IDs do not match the native random design columns.", call. = FALSE)
  list(variable = group, ids = ids, row_coordinate = coordinate)
}

# This private snapshot contains values, verified recipes and position maps.
# The native fitter and its uncertainty and admission checks remain separate.
.msr_gdr_prepare_native <- function(model, event_keys, event_id, interval_units, tmb_settings) {
  if (!identical(class(model)[1L], "unmarkedFitGDR") || length(model@sitesRemoved))
    stop("Prepare an ordinary native GDR model with complete retained rows.", call. = FALSE)
  settings <- .msr_gdr_TMB_settings(model, tmb_settings)
  original_policy <- attr(model, "msr_gdr_tmb_diagnostics")$workflow$policy
  inner <- c("random_start", "inner_control")
  if (!is.null(original_policy) && !identical(settings$policy[inner], original_policy[inner]))
    stop("Prepared GDR reconstruction requires the original inner-mode policy.", call. = FALSE)
  frame <- model@data
  rows <- .msr_gdr_prepared_rows(frame, event_keys, event_id, interval_units)
  native <- .msr_gdr_TMB_preparation(model)
  if (!identical(native$K, model@K) || !identical(native$y_dist, as.vector(t(frame@yDistance))) ||
      !identical(native$y_rem, as.vector(t(frame@yRemoval))))
    stop("Prepared GDR requires the exact native data, types and response order.", call. = FALSE)

  components <- names(.msr_gdr_active_formulas(model))
  ns <- asNamespace("unmarked")
  recipes <- groups <- stats::setNames(vector("list", length(components)), components)
  for (type in components) {
    raw <- get("get_covariates", ns)(model, type)
    formula <- model@formlist[[type]]
    basis <- .msr_unmarked_prediction_basis(model, type)
    basis$terms <- .msr_cv_prediction_terms(basis$terms, raw, qualify = TRUE, allow_learning = FALSE)
    environment(basis$terms) <- baseenv()
    X <- native[[paste0("X_", type)]]
    Z <- native[[paste0("Z_", type)]]
    rebuilt <- .msr_unmarked_fixed_design(basis, raw)$X
    if (nrow(raw) != nrow(rows[[type]]) || any(!is.finite(X)) ||
        !identical(typeof(X), typeof(rebuilt)) ||
        !isTRUE(all.equal(rebuilt, X, tolerance = 1e-12, check.attributes = FALSE)))
      stop("Trained GDR recipes must reproduce the original component rows and design.", call. = FALSE)
    recipes[[type]] <- basis
    groups[type] <- list(.msr_gdr_prepared_groups(formula, raw, Z))
  }
  # TMB stores its original constructor seeds in parameters. Use the native
  # random coordinate indices to recover the ordered parameter block names.
  env <- model@TMB$env
  parameters <- env$parameters
  random <- unique(names(env$last.par)[env$random])
  active_random <- paste0("b_", names(groups)[!vapply(groups, is.null, logical(1))], recycle0 = TRUE)
  if (!setequal(random, active_random) ||
      any(vapply(names(groups), function(type)
        length(parameters[[paste0("b_", type)]]) != length(groups[[type]]$ids), logical(1))) ||
      length(model@TMB$par) != settings$parameter_count)
    stop("Prepared GDR group coordinates disagree with the native parameter layout.", call. = FALSE)
  payload <- list(native = native, parameters = parameters, random = random,
    random_index = env$random, policy = settings$policy, fixed_start = model@TMB$par,
    parameter_count = settings$parameter_count,
    recipes = recipes, event_keys = event_keys[event_id], rows = rows, groups = groups,
    # TMB converts constructor integers and sparse storage internally. Bind
    # both the original typed input and the representation on the fitted tape.
    tape_data_sha256 = .msr_gdr_hash(env$data), provider = .msr_gdr_field_provider(),
    r_version = R.version.string)
  structure(list(format_version = 1L, payload = payload, fingerprint = .msr_gdr_hash(payload)),
    class = "msr_gdr_prepared_native")
}

.msr_gdr_rebuild_native <- function(prepared) {
  if (!inherits(prepared, "msr_gdr_prepared_native") || !identical(prepared$format_version, 1L) ||
      !identical(.msr_gdr_hash(prepared$payload), prepared$fingerprint))
    stop("Prepared GDR input identity changed.", call. = FALSE)
  input <- prepared$payload
  invisible(loadNamespace("unmarked"))
  invisible(loadNamespace("TMB"))
  if (!identical(input$r_version, R.version.string) || !identical(input$provider, .msr_gdr_field_provider()))
    stop("Prepared GDR input requires its original native provider.", call. = FALSE)
  old_threads <- TMB::openmp(DLL = "unmarked_TMBExports")
  on.exit(TMB::openmp(old_threads, DLL = "unmarked_TMBExports"))
  TMB::openmp(1L, DLL = "unmarked_TMBExports")
  objective <- .msr_gdr_TMB_adfun(input$native$model, input$native[names(input$native) != "model"],
    input$parameters, input$random, input$policy)
  if (!identical(.msr_gdr_hash(objective$env$data), input$tape_data_sha256) ||
      !identical(objective$par, input$fixed_start) ||
      !identical(objective$env$parameters, input$parameters) ||
      !identical(objective$env$random, input$random_index) ||
      length(objective$par) != input$parameter_count)
    stop("Reconstructed GDR objective changed native inputs or coordinates.", call. = FALSE)
  list(objective = objective, input_fingerprint = prepared$fingerprint,
    parameter_count = input$parameter_count, inference_eligible = FALSE,
    uncertainty = "none_objective_reconstruction_only")
}
