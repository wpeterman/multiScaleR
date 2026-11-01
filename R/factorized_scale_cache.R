# Standalone numeric tape for a single recoded group. Its population layout
# and immutable group payload belong to the caller; it is not a fitted model.
.msr_factorized_group_log_tape <- function(group, initial, unpack) {
  .msr_factorized_group_check(group, TRUE)
  if (!is.numeric(initial) || !length(initial) || any(!is.finite(initial)) || !is.function(unpack))
    stop("Supply a finite numeric group layout and an unpack function.", call. = FALSE)
  p <- unpack(initial)
  sizes <- lengths(group$parameters)
  sizes[c("scale_z", "response_z")] <- c(1L, 0L)
  if (!is.list(p) || is.null(names(p)) || anyDuplicated(names(p)) ||
      !identical(names(p), names(group$parameters)) ||
      !identical(lengths(p), sizes) ||
      !all(vapply(p, is.numeric, TRUE)) || sum(lengths(p)) != length(initial) ||
      !identical(as.numeric(unlist(p, use.names = FALSE)), as.numeric(initial)) ||
      length(p$scale_z) != 1L || length(p$response_z))
    stop("Supply a matching one-latent integrated-response group layout.", call. = FALSE)
  force(group); force(unpack)
  configuration <- RTMB::TapeConfig()
  on.exit(RTMB::TapeConfig(configuration))
  RTMB::MakeTape(function(x) {
    p <- unpack(x)
    -.msr_factorized_group_joint_nll(group, p, TRUE)
  }, initial)
}

# Numeric evaluation of a reusable, small joint-density tape. These tapes are
# computational caches, not portable fit results or independent R references.
# Reuse one cache only with its unchanged prepared specification. Supplying it
# to conditional reconstruction avoids repeated R evaluation inside quadrature;
# quadrature bounds, tolerances and the default ordinary audit stay unchanged.
.msr_factorized_conditional_cache <- function(spec) {
  .msr_factorized_validate(spec)
  .msr_scale_surface_ram(spec)
  if (!requireNamespace("RTMB", quietly = TRUE))
    stop("Conditional density caching requires optional package 'RTMB'.", call. = FALSE)
  term <- spec$scales[[1L]]$random
  if (length(term) != 1L)
    stop("Conditional density caching requires one annual random-effect term.", call. = FALSE)
  template <- spec$parameters
  template$scale_z <- 0; template$response_z <- numeric()
  sizes <- lengths(template)
  blocks <- rep(names(template), sizes)
  initial <- unlist(template, use.names = FALSE)
  unpack <- function(x) {
    p <- template; at <- 0L
    for (j in seq_along(p)) {
      n <- sizes[j]
      p[[j]] <- if (n) x[at + seq_len(n)] else numeric()
      at <- at + n
    }
    p
  }
  configuration <- RTMB::TapeConfig()
  on.exit(RTMB::TapeConfig(configuration))
  tapes <- lapply(seq_along(term[[1L]]$levels), function(i) {
    group <- .msr_factorized_group(spec, i); force(group)
    .msr_factorized_group_log_tape(group, initial, unpack)
  })
  log_density <- function(parameters, group) {
    if (!is.list(parameters) || is.null(names(parameters)) ||
        anyDuplicated(names(parameters)) || any(!names(parameters) %in% names(template)) ||
        length(parameters$scale_z) ||
        length(group) != 1L || !is.numeric(group) || !is.finite(group) ||
        group != floor(group) || group < 1L || group > length(tapes))
      stop("Supply a compatible population parameter layout and annual group.", call. = FALSE)
    p <- lapply(seq_along(template), function(j) {
      nm <- names(template)[j]
      if (nm == "scale_z") 0 else if (sizes[j] == 0L && is.null(parameters[[nm]]))
        numeric() else parameters[[nm]]
    })
    if (!identical(lengths(p), unname(sizes)) ||
        !all(vapply(p, is.numeric, TRUE)) || any(!is.finite(unlist(p))))
      stop("Supply a compatible finite population parameter layout.", call. = FALSE)
    input <- unlist(p, use.names = FALSE)
    latent <- which(blocks == "scale_z"); tape <- tapes[[group]]
    function(z) {
      if (!is.numeric(z) || length(z) != 1L || !is.finite(z))
        stop("Supply one finite latent annual effect.", call. = FALSE)
      input[latent] <- z
      as.numeric(tape(input))
    }
  }
  structure(list(specification = spec, log_density = log_density,
    parameter_sizes = sizes, groups = length(tapes),
    backend = "Numeric RTMB joint-density tape evaluation"),
    class = "msr_factorized_conditional_cache")
}
