# Run from the worktree root: Rscript run.R grid|point OUTPUT
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, args[1L] %in% c("grid", "point"))
panel <- args[1L]
source_root <- normalizePath(".", mustWork = TRUE)
output <- normalizePath(args[2L], mustWork = FALSE)
n_each <- if (panel == "grid") 16L else 32L
B <- 99L
level <- .90
fixed_grid <- c(90, 140, 200, 280, 400, 550, 680)
coarse_nodes <- c(50, 140, 280, 550, 800)
scenarios <- c("interior", "weak", "near_bound")
plan <- expand.grid(scenario = scenarios, replicate = seq_len(n_each),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
plan$seed <- (if (panel == "grid") 1700000L else 1800000L) + seq_len(nrow(plan)) * 100L
plan$bootstrap_seed <- plan$seed + 1L
plan$truth_seed <- plan$seed + 2L
plan$id <- sprintf("%s-%02d", plan$scenario, plan$replicate)

capture <- function(expr) {
  messages <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  list(value = value, warnings = unique(messages))
}

atomic <- function(x, file) {
  tmp <- tempfile(tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  if (!file.rename(tmp, file)) stop("Cannot publish ", file)
}

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
stopifnot(file.exists(script), requireNamespace("digest", quietly = TRUE),
  requireNamespace("MASS", quietly = TRUE))
package_files <- c("DESCRIPTION", "NAMESPACE",
  file.path("R", list.files(file.path(source_root, "R"), pattern = "[.]R$")),
  file.path("src", list.files(file.path(source_root, "src"),
    pattern = "[.](c|cpp|h|dll|so)$")))
hash_files <- function() {
  files <- file.path(source_root, package_files)
  stats::setNames(vapply(files, digest::digest, "", file = TRUE,
    algo = "sha256"), package_files)
}

# This is the pilot's independent, fixed generator with new dataset seeds.
simulate_case <- function(case) {
  set.seed(case$seed)
  n <- 110L
  sigma <- if (case$scenario == "near_bound") .73 else .24
  effect <- if (case$scenario == "weak") .12 else .95
  size <- if (case$scenario == "near_bound") Inf else 2
  d <- seq(.015, 1, length.out = 21L)
  count <- matrix(sample(1:5, n * length(d), replace = TRUE), n)
  radial <- outer(rnorm(n), sin(7 * d)) + outer(rnorm(n), cos(12 * d)) +
    outer(rnorm(n), d) + matrix(rnorm(n * length(d), sd = .25), n)
  V <- radial * count
  habitat <- function(scale) {
    w <- exp(-(d^2 - min(d)^2) / (2 * scale^2))
    as.numeric(base::scale(drop(V %*% w) / drop(count %*% w)))
  }
  data <- data.frame(hab = habitat(.23), z = rnorm(n), effort = runif(n, .8, 1.2),
    site = factor(rep(seq_len(20L), each = 5L, length.out = n)))
  mu <- exp(.4 + effect * habitat(sigma) + .3 * data$z + log(data$effort))
  data$y <- if (is.infinite(size)) rpois(n, mu) else rnbinom(n, size = size, mu = mu)
  inputs <- list(raw_cov = NULL, d_list = NULL, min_D = 50, max_D = 800,
    unit_conv = 1000, kernel = "gaussian",
    binned = list(bin_rep = d, vsum = list(hab = V), csum = list(hab = count),
      point_ids = rownames(data), n_points = n))
  list(data = data, inputs = inputs, sigma = sigma * 1000, effect = effect, size = size)
}

ci95 <- function(k, n) if (n == 0L) rep(NA_real_, 2L) else
  as.numeric(stats::binom.test(k, n)$conf.int)


# Use an independent reconstruction of components and tail probabilities.
components <- function(sigma, accept) {
  if (anyNA(accept)) return(NULL)
  accepted <- which(accept)
  if (!length(accepted)) return(data.frame(lower = numeric(), upper = numeric(),
    lower_outside = numeric(), upper_outside = numeric(),
    lower_bound_limited = logical(), upper_bound_limited = logical()))
  groups <- split(accepted, cumsum(c(TRUE, diff(accepted) > 1L)))
  do.call(rbind, lapply(groups, function(i) {
    l <- min(i); h <- max(i)
    data.frame(lower = sigma[l], upper = sigma[h],
      lower_outside = if (l > 1L) sigma[l - 1L] else NA_real_,
      upper_outside = if (h < length(sigma)) sigma[h + 1L] else NA_real_,
      lower_bound_limited = l == 1L, upper_bound_limited = h == length(sigma))
  }))
}

tails <- function(result, keep) {
  grid <- result$grid
  pairs <- result$bootstrap_records
  out <- lapply(grid$sigma, function(s) {
    p <- pairs[pairs$generating_sigma == s & pairs$replicate <= keep, , drop = FALSE]
    if (!nrow(p)) return(rep(NA_real_, 5L))
    stopifnot(nrow(p) == keep, identical(p$replicate, seq_len(keep)))
    if (any(!p$ok)) return(rep(NA_real_, 5L))
    k <- sum(p$lr >= grid$lr[grid$sigma == s] - 1e-8)
    p_value <- (1 + k) / (keep + 1)
    mc <- (1 + keep * ci95(k, keep)) / (keep + 1)
    c(p_value, mc, 10 * (1 + k) > keep + 1, mc[1L] <= .1 && mc[2L] > .1)
  })
  out <- as.data.frame(do.call(rbind, out))
  names(out) <- c("p", "mc_lower", "mc_upper", "accept", "ambiguous")
  out$accept <- as.logical(out$accept)
  out$ambiguous <- as.logical(out$ambiguous)
  if (keep == B) {
    stopifnot(isTRUE(all.equal(out$p, grid$bootstrap_p, tolerance = 1e-14)),
      isTRUE(all.equal(out$mc_lower, grid$mc_lower, tolerance = 1e-14)),
      isTRUE(all.equal(out$mc_upper, grid$mc_upper, tolerance = 1e-14)),
      identical(out$accept, grid$bootstrap_accept),
      identical(out$ambiguous, grid$mc_ambiguous))
  }
  cbind(sigma = grid$sigma, out)
}

contains <- function(set, sigma) {
  if (is.null(set)) return(NA)
  any(set$lower <= sigma & set$upper >= sigma)
}

set_metrics <- function(sigma, tail, truth, comparison) {
  set <- components(sigma, tail$accept)
  absent <- is.null(set)
  span <- function(s) sum(log(s$upper) - log(s$lower))
  lower <- function(s) log(s$lower)
  upper <- function(s) log(s$upper)
  common <- if (absent) NA_real_ else sum(pmax(0,
    outer(upper(set), upper(comparison), pmin) -
      outer(lower(set), lower(comparison), pmax)))
  union <- if (absent) NA_real_ else span(set) + span(comparison) - common
  bracket <- if (absent) NA_real_ else c(log(set$lower / set$lower_outside),
    log(set$upper_outside / set$upper))
  clear <- !tail$ambiguous
  certain <- if (absent) NULL else components(sigma, tail$accept & clear)
  possible <- if (absent) NULL else components(sigma, tail$accept | !clear)
  list(available = !absent, covers = contains(set, truth),
    certain_covers = contains(certain, truth), possible_covers = contains(possible, truth),
    components = if (absent) NA_integer_ else nrow(set),
    empty = if (absent) NA else !nrow(set),
    disconnected = if (absent) NA else nrow(set) > 1L,
    singleton_component = if (absent) NA else any(set$lower == set$upper),
    any_bound = if (absent) NA else any(set$lower_bound_limited | set$upper_bound_limited),
    full_range = if (absent) NA else nrow(set) == 1L &&
      set$lower_bound_limited && set$upper_bound_limited,
    log_span = if (absent) NA_real_ else span(set) / diff(log(range(sigma))),
    log_endpoint_bracket_max = if (all(is.na(bracket))) NA_real_ else max(bracket, na.rm = TRUE),
    profile_jaccard_log_length = if (!is.finite(union) || union == 0) NA_real_ else common / union,
    ambiguous_nodes = sum(tail$ambiguous, na.rm = TRUE),
    grid_nodes = length(sigma), tested_nodes = sum(!is.na(tail$accept)), set = set)
}

fit_case <- function(case) {
  clock <- proc.time()[["elapsed"]]
  f <- simulate_case(case)
  x <- list(case = case, true_sigma = f$sigma, effect = f$effect, size = f$size,
    status = "initial_fit_failed", initial = NULL, fit_starts = NULL,
    fitted_sigma = NA_real_, exact_poisson = NA, fit_warnings = NULL,
    grid_error = NULL, point_error = NULL, grid_result = NULL, point_result = NULL,
    sets = NULL, point = NULL, seconds = NA_real_)
  done <- function() {
    x$seconds <- proc.time()[["elapsed"]] - clock
    x
  }
  initial <- capture(stats::glm(y ~ hab + z + offset(log(effort)), data = f$data,
    family = MASS::negative.binomial(2)))
  x$initial <- list(warnings = initial$warnings,
    error = if (inherits(initial$value, "error")) conditionMessage(initial$value))
  if (inherits(initial$value, "error") || !isTRUE(initial$value$converged)) return(done())
  x$status <- "point_fit_failed"
  fitted <- capture(multiScale_optim(initial$value, f$inputs, backend = "joint_nb2",
    n_cores = 1L, verbose = FALSE,
    profile_control = list(starts = matrix(c(.075, .20, .43, .71), ncol = 1L))))
  x$fit_warnings <- fitted$warnings
  if (inherits(fitted$value, "error")) {
    x$fit_error <- conditionMessage(fitted$value)
    x$fit_starts <- fitted$value$start_results
    return(done())
  }
  fit <- fitted$value
  x$fitted_sigma <- fit$scale_est$Mean[1L]
  x$exact_poisson <- fit$dispersion_est$exact_poisson
  x$fit_starts <- fit$profile_fit$start_results

  if (panel == "grid") {
    grid <- capture(scale_inference(fit, sigma_values = fixed_grid, level = level,
      bootstrap = B, seed = case$bootstrap_seed, n_starts = 2L))
    x$grid_warnings <- grid$warnings
    if (inherits(grid$value, "error")) {
      x$grid_error <- conditionMessage(grid$value)
    } else {
      x$grid_result <- grid$value
      g <- grid$value
      full <- tails(g, B)
      reduced <- tails(g, 49L)
      stopifnot(all(fixed_grid %in% full$sigma), all(g$bounds %in% full$sigma))
      coarse <- vapply(full$sigma, function(s) any(abs(s - c(coarse_nodes,
        g$estimate)) < 1e-9), logical(1))
      stopifnot(sum(coarse) >= 5L, sum(coarse) < nrow(full))
      x$sets <- list(
        fine99 = set_metrics(full$sigma, full, f$sigma, g$profile_set),
        coarse99 = set_metrics(full$sigma[coarse], full[coarse, ], f$sigma, g$profile_set),
        fine49 = set_metrics(reduced$sigma, reduced, f$sigma, g$profile_set),
        coarse49 = set_metrics(reduced$sigma[coarse], reduced[coarse, ], f$sigma, g$profile_set))
      stopifnot(isTRUE(all.equal(x$sets$fine99$set, g$bootstrap_set,
        tolerance = 1e-14, check.attributes = FALSE)))
      x$tails49 <- reduced
      x$coarse <- coarse
    }
  }

  # The direct test has a separately seeded stream and does not contribute
  # a node to the full-grid set.
  x$status <- "point_inference_failed"
  point <- capture(scale_inference(fit, sigma_values = f$sigma, level = level,
    bootstrap = B, bootstrap_values = f$sigma, seed = case$truth_seed, n_starts = 2L))
  x$point_warnings <- point$warnings
  if (inherits(point$value, "error")) {
    x$point_error <- conditionMessage(point$value)
    return(done())
  }
  x$point_result <- point$value
  point99 <- tails(point$value, B)
  point49 <- tails(point$value, 49L)
  truth <- point99$sigma == f$sigma
  stopifnot(sum(truth) == 1L, is.null(point$value$bootstrap_set),
    nrow(point$value$bootstrap_records) == B)
  x$point <- rbind(cbind(draws = B, point99[truth, ]),
    cbind(draws = 49L, point49[truth, ]))
  x$status <- if (is.na(point99$accept[truth])) "point_pair_failed" else "point_complete"
  done()
}

main <- function() {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  cases <- file.path(output, "cases")
  dir.create(cases, showWarnings = FALSE)
  lock <- file.path(output, "run.lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("The output directory is already in use.")
  stopifnot(identical(dirname(normalizePath(lock, winslash = "/")),
    normalizePath(output, winslash = "/")))
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  sources <- hash_files()
  pkgload::load_all(path = source_root, compile = FALSE, quiet = TRUE)
  stopifnot(identical(sources, hash_files()))
  config <- list(panel = panel, B = B, n_each = n_each, level = level,
    grid = fixed_grid, coarse = coarse_nodes, plan = plan, R = R.version.string,
    source = sources, script = digest::digest(script, file = TRUE, algo = "sha256"))
  metadata <- file.path(output, "configuration.rds")
  if (file.exists(metadata)) stopifnot(identical(readRDS(metadata), config)) else {
    atomic(config, metadata)
    utils::write.csv(plan, file.path(output, "plan.csv"), row.names = FALSE)
    utils::write.csv(data.frame(path = names(sources), sha256 = unname(sources)),
      file.path(output, "source_sha256.csv"), row.names = FALSE)
    writeLines(capture.output(sessionInfo()), file.path(output, "session.txt"))
  }
  paths <- file.path(cases, paste0(plan$id, ".rds"))
  for (i in seq_len(nrow(plan))) {
    if (file.exists(paths[i])) next
    if (file.exists(file.path(output, "STOP"))) break
    id <- plan$id[i]
    cat(format(Sys.time()), panel, "starting", id, "\n")
    flush.console()
    x <- fit_case(plan[i, , drop = FALSE])
    atomic(x, paths[i])
    cat(format(Sys.time()), panel, "finished", id, x$status,
      "seconds", x$seconds, "grid", x$sets$fine99$available,
      "covered", x$sets$fine99$covers, "point", x$point$accept[1L],
      "error", x$error, "\n")
    flush.console()
  }
  stopifnot(identical(sources, hash_files()))
}

main()
