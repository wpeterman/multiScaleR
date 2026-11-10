# Numerical policy is separate from the likelihood and population priors.
# Keep the native cached-start defaults by omitting unrequested MakeADFun args.
.msr_gdr_TMB_policy <- function(settings = list()) {
  requested <- settings[intersect(c("random.start", "inner.control", "record_memory"), names(settings))]
  policy <- if (is.null(settings$random.start)) "cached_best" else settings$random.start
  if (!is.character(policy) || !is.null(attributes(policy)) || length(policy) != 1L || is.na(policy) ||
      !policy %in% c("cached_best", "zero"))
    stop("random.start must be 'cached_best' or 'zero'; arbitrary expressions and latent vectors require separate qualification.", call. = FALSE)
  control <- settings$inner.control
  if (!is.null(control)) {
    required <- c("maxit", "tol", "grad.tol", "step.tol")
    if (!is.list(control) || !identical(attributes(control), list(names = names(control))) || is.null(names(control)) ||
        anyNA(names(control)) || anyDuplicated(names(control)) ||
        !setequal(names(control), required) ||
        any(!vapply(control, function(x) is.numeric(x) && is.null(attributes(x)) &&
          length(x) == 1L && is.finite(x) && x > 0, logical(1))) ||
        control$maxit != floor(control$maxit) || control$maxit > .Machine$integer.max)
      stop("inner.control must explicitly contain positive finite maxit, tol, grad.tol and step.tol; maxit is an integer.", call. = FALSE)
  }
  memory <- if (is.null(settings$record_memory)) FALSE else settings$record_memory
  validate_scalar_logical(memory, "record_memory")
  if (!is.null(attributes(memory))) stop("record_memory must be one plain logical value.", call. = FALSE)
  if (memory && !requireNamespace("ps", quietly = TRUE))
    stop("record_memory = TRUE requires the optional package 'ps'.", call. = FALSE)
  list(random_start = policy, inner_control = control, record_memory = memory,
    requested = requested, defaults = "native_TMB_provider_bound_when_not_overridden")
}

.msr_gdr_TMB_adfun <- function(model, data, params, random, policy) {
  args <- list(data = c(model = model, data), parameters = params,
    random = random, silent = TRUE, DLL = "unmarked_TMBExports")
  if (identical(policy$random_start, "zero"))
    args$random.start <- expression(rep(0, length(random)))
  if (!is.null(policy$inner_control)) args$inner.control <- policy$inner_control
  # Evaluate the random-start expression only inside TMB's complete joint
  # environment, where random is the native coordinate index vector.
  do.call(TMB::MakeADFun, args, quote = TRUE)
}

# These are process observations at boundaries, not isolated native-tape costs.
# maxrss is a lifetime maximum and cannot be assigned to the current phase.
.msr_gdr_TMB_phase <- function(name, started, policy) {
  memory <- list(requested = policy$record_memory, available = FALSE,
    scope = "whole_current_R_process_at_phase_boundary_not_phase_peak_or_isolated_tape",
    rss_bytes = NA_real_, lifetime_maxrss_bytes = NA_real_)
  if (policy$record_memory) memory <- tryCatch({
    observed <- ps::ps_memory_full_info(ps::ps_handle())
    memory$available <- is.finite(observed[["rss"]])
    memory$rss_bytes <- unname(observed[["rss"]])
    memory$lifetime_maxrss_bytes <- if ("maxrss" %in% names(observed)) unname(observed[["maxrss"]]) else NA_real_
    memory$provider <- list(package = "ps", version = as.character(utils::packageVersion("ps")))
    memory
  }, error = function(e) {memory$reason <- conditionMessage(e); memory})
  observed_elapsed <- unname(proc.time()[["elapsed"]])
  list(name = name, timestamp_UTC = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
    PID = Sys.getpid(), elapsed_seconds = observed_elapsed - started,
    process_elapsed_seconds = observed_elapsed, memory = memory)
}

# Fixed positions preserve repeated native display names. They are not a
# physical coefficient/state identity map or constructor-coordinate order.
.msr_gdr_TMB_point <- function(opt, model, data, params, random, policy, phases) {
  structure(list(schema_version = 1L, type = "native_TMB_GDR_optimizer_point",
    status = "optimization_returned_before_covariance", inference_eligible = FALSE,
    covariance_status = "not_started", qualification_status = "not_qualified",
    optimizer = .msr_gdr_optimizer_checkpoint(opt),
    coordinates = data.frame(position = seq_along(opt$par),
      native_display_label = if (is.null(names(opt$par))) rep("", length(opt$par)) else names(opt$par)),
    coordinate_scale = "native_fixed_optimizer_coordinates_not_constructor_or_physical_state_order",
    native_data_sha256 = .msr_gdr_hash(c(model = model, data)),
    finite_cutoff = list(value = data$K, storage = typeof(data$K)),
    parameter_layout_sha256 = .msr_gdr_hash(list(names(params), lengths(params), random)),
    policy = policy, phases = phases, models_or_tapes_saved = FALSE,
    latent_states_saved = FALSE, covariance_saved = FALSE), class = "msr_gdr_optimizer_point")
}
