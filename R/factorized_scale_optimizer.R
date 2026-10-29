# A coarse optimizing integral may fail at otherwise valid parameter values.
# Retry only the same vector with a fixed, finer partition; the independent
# final audits and all acceptance thresholds remain the fitter's original ones.
.msr_factorized_optimizer <- function(spec, parameters = spec$parameters, step = 1) {
  if (!is.numeric(step) || length(step) != 1L || is.na(step) ||
      !step %in% c(1, 2, 4, 8)) stop("Invalid optimizing integration step.", call. = FALSE)
  if (step == 1 || !length(spec$scales[[1L]]$random))
    return(.msr_factorized_adfun(spec, parameters, step = step))
  .msr_factorized_validate(spec)
  parameters$scale_z <- parameters$response_z <- numeric()
  sizes <- lengths(parameters)
  initial <- unlist(parameters, use.names = FALSE)
  names(initial) <- rep(names(parameters), sizes)
  unpack <- function(x) {
    p <- parameters; at <- 0L
    for (i in seq_along(p)) {
      n <- sizes[i]
      if (n) p[[i]][] <- as.numeric(x[at + seq_len(n)])
      at <- at + n
    }
    p
  }
  steps <- c(8L, 4L, 2L, 1L); steps <- steps[steps <= step]
  state <- new.env(parent = emptyenv())
  state$current <- 1L; state$active <- NULL
  state$integration_events <- list(); state$integration_generation <- 0L
  state$integration_current_step <- steps[1L]
  state$last.par <- initial
  state$parList <- function(x, par = NULL) unpack(x)
  evaluate <- function(x, phase) {
    if (length(x) != length(initial) || any(!is.finite(x)))
      stop("Invalid factorized parameter vector.", call. = FALSE)
    repeat {
      construction <- is.null(state$active)
      value <- tryCatch({
        if (construction) {
          state$active <- .msr_factorized_adfun(spec, unpack(x), step = steps[state$current])
          state$integration_generation <- state$integration_generation + 1L
        }
        answer <- state$active[[phase]](x)
        if (any(!is.finite(answer))) stop(paste("Nonfinite optimizing", phase), call. = FALSE)
        list(answer = answer)
      }, error = function(e) list(error = conditionMessage(e)))
      if (!length(value$error)) {
        state$last.par <- x
        return(value$answer)
      }
      next_step <- if (state$current < length(steps)) steps[state$current + 1L] else NA_integer_
      state$integration_events[[length(state$integration_events) + 1L]] <- list(
        from = steps[state$current], to = next_step, phase = phase,
        during_construction = construction, parameters = x, error = value$error)
      # An exception can leave integration internals at the failed point.
      # Discard the object before constructing at the identical vector.
      state$active <- NULL
      if (is.na(next_step)) stop(value$error, call. = FALSE)
      gc(FALSE)
      state$current <- state$current + 1L
      state$integration_current_step <- steps[state$current]
    }
  }
  list(par = initial, env = state,
    fn = function(x) evaluate(x, "fn"), gr = function(x) evaluate(x, "gr"),
    he = function(x) evaluate(x, "he"))
}
