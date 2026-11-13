# Analyze archived checkpoints without loading or refitting the package.
# Rscript summarize.R GRID_OUTPUT POINT_OUTPUT REPORT_OUTPUT [GRID_SHARD ...]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3L)
out <- normalizePath(args[3L], mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ci95 <- function(k, n) if (n == 0L) rep(NA_real_, 2L) else
  as.numeric(stats::binom.test(k, n)$conf.int)
metric <- function(x, n = length(x)) {
  k <- sum(x, na.rm = TRUE)
  limits <- ci95(k, n)
  c(successes = k, denominator = n, rate = if (n) k/n else NA_real_,
    lower95 = limits[1L], upper95 = limits[2L])
}
num <- function(x, f) if (!length(x) || all(is.na(x))) NA_real_ else f(x, na.rm = TRUE)
or <- function(x, value) if (is.null(x) || !length(x)) value else x
same_components <- function(a, b) {
  if (is.null(a) || is.null(b)) return(NA)
  identical(unname(as.matrix(a[c("lower", "upper")])),
    unname(as.matrix(b[c("lower", "upper")])))
}
union_overlap <- function(a, b) {
  if (is.null(a) || is.null(b)) return(NA_real_)
  span <- function(s) sum(log(s$upper) - log(s$lower))
  common <- sum(pmax(0, outer(log(a$upper), log(b$upper), pmin) -
    outer(log(a$lower), log(b$lower), pmax)))
  total <- span(a) + span(b) - common
  if (total <= 0) NA_real_ else common / total
}

load_panel <- function(path, shards = character()) {
  config <- readRDS(file.path(path, "configuration.rds"))
  plan <- config$plan
  paths <- file.path(path, "cases", paste0(plan$id, ".rds"))
  assigned <- plan$id[file.exists(paths)]
  for (shard in shards) {
    meta <- readRDS(file.path(shard, "shard.rds"))
    stopifnot(identical(meta$source, config$source), identical(meta$script, config$script),
      all(meta$wanted %in% plan$id), !anyDuplicated(meta$wanted),
      !any(meta$wanted %in% assigned))
    idx <- match(meta$wanted, plan$id)
    paths[idx] <- file.path(shard, "cases", paste0(meta$wanted, ".rds"))
    assigned <- c(assigned, meta$wanted)
  }
  present <- file.exists(paths)
  outcomes <- rep(list(NULL), nrow(plan))
  outcomes[present] <- lapply(paths[present], readRDS)
  names(outcomes) <- plan$id
  stopifnot(all(vapply(which(present), function(i)
    identical(outcomes[[i]]$case$id, plan$id[i]), logical(1))))
  complete_wave <- vapply(split(present, plan$replicate), all, logical(1))
  prefix_end <- sum(cumprod(as.integer(complete_wave)))
  plan$prefix <- present & plan$replicate <= prefix_end
  plan$present <- present
  plan$panel <- config$panel
  list(plan = plan, cases = outcomes, configuration = config)
}

extract_datasets <- function(panel) {
  records <- lapply(seq_len(nrow(panel$plan)), function(i) {
    p <- panel$plan[i, ]
    x <- panel$cases[[i]]
    g <- x$grid_result
    z <- x$point_result
    data.frame(p, status = or(x$status, "pending"), true_sigma = or(x$true_sigma, NA_real_),
      fitted_sigma = or(x$fitted_sigma, NA_real_), exact_poisson = or(x$exact_poisson, NA),
      seconds = or(x$seconds, NA_real_), profile_status = or(g$diagnostics$profile_status, NA_character_),
      grid_status = or(g$diagnostics$bootstrap_status, NA_character_),
      point_status = or(z$diagnostics$bootstrap_status, NA_character_),
      initial_error = or(x$initial$error, ""), fit_error = or(x$fit_error, ""),
      grid_error = or(x$grid_error, ""), point_error = or(x$point_error, ""),
      initial_warnings = paste(x$initial$warnings, collapse = " | "),
      fit_warnings = paste(x$fit_warnings, collapse = " | "),
      grid_warnings = paste(x$grid_warnings, collapse = " | "),
      point_warnings = paste(x$point_warnings, collapse = " | "),
      grid_pairs = or(nrow(g$bootstrap_records), 0L),
      grid_failed_pairs = sum(!as.logical(g$bootstrap_records$ok)),
      grid_warned_pairs = sum(nzchar(g$bootstrap_records$warnings)),
      point_pairs = or(nrow(z$bootstrap_records), 0L),
      point_failed_pairs = sum(!as.logical(z$bootstrap_records$ok)),
      point_warned_pairs = sum(nzchar(z$bootstrap_records$warnings)))
  })
  do.call(rbind, records)
}

extract_sets <- function(panel) {
  methods <- c("fine99", "coarse99", "fine49", "coarse49")
  blank <- list(available = FALSE, covers = NA, certain_covers = NA,
    possible_covers = NA, components = NA_integer_, empty = NA,
    disconnected = NA, singleton_component = NA, any_bound = NA,
    full_range = NA, log_span = NA_real_, log_endpoint_bracket_max = NA_real_,
    profile_jaccard_log_length = NA_real_, ambiguous_nodes = NA_integer_,
    grid_nodes = NA_integer_, tested_nodes = NA_integer_)
  records <- lapply(seq_len(nrow(panel$plan)), function(i) {
    p <- panel$plan[i, ]
    x <- panel$cases[[i]]
    y <- lapply(methods, function(method) {
      value <- x$sets[[method]]
      if (is.null(value)) value <- blank else value$set <- NULL
      data.frame(p, method, as.data.frame(value))
    })
    do.call(rbind, y)
  })
  do.call(rbind, records)
}

extract_points <- function(panel) {
  records <- lapply(seq_len(nrow(panel$plan)), function(i) {
    p <- panel$plan[i, ]
    x <- panel$cases[[i]]
    if (is.null(x$point)) x$point <- data.frame(draws = c(99L, 49L),
      sigma = or(x$true_sigma, NA_real_), p = NA_real_, mc_lower = NA_real_,
      mc_upper = NA_real_, accept = NA, ambiguous = NA)
    z <- x$point_result$grid
    profile <- if (is.null(z)) NA else z$profile_accept[z$sigma == x$true_sigma]
    stopifnot(length(profile) == 1L)
    p <- p[rep(1L, nrow(x$point)), , drop = FALSE]
    rownames(p) <- rownames(x$point) <- NULL
    cbind(p, x$point, chi_square_accept = profile)
  })
  do.call(rbind, records)
}

extract_components <- function(panel) {
  rows <- lapply(panel$cases, function(x) {
    if (is.null(x)) return(NULL)
    result <- lapply(names(x$sets), function(method) {
      set <- x$sets[[method]]$set
      if (is.null(set) || !nrow(set)) return(NULL)
      data.frame(id = x$case$id, method, component = seq_len(nrow(set)), set)
    })
    do.call(rbind, result)
  })
  do.call(rbind, rows)
}

extract_nodes <- function(panel) {
  rows <- lapply(panel$cases, function(x) {
    g <- x$grid_result
    if (is.null(g)) return(NULL)
    data.frame(id = x$case$id, g$grid, nested_coarse = x$coarse,
      p49 = x$tails49$p, mc_lower49 = x$tails49$mc_lower,
      mc_upper49 = x$tails49$mc_upper, accept49 = x$tails49$accept,
      ambiguous49 = x$tails49$ambiguous)
  })
  do.call(rbind, rows)
}

point_comparison <- function(data, scope) {
  data$scope <- scope
  grouped(data, c("scope", "panel", "scenario", "draws"), function(x) {
    y <- x[x$present & (scope == "all_attempted" | x$prefix), ]
    z <- y[!is.na(y$chi_square_accept) & !is.na(y$accept), ]
    chi <- y$chi_square_accept[!is.na(y$chi_square_accept)]
    outcomes <- list(chi_square_acceptance = metric(chi),
      paired_disagreement = metric(z$chi_square_accept != z$accept),
      bootstrap_only_acceptance = metric(z$accept & !z$chi_square_accept),
      chi_square_only_acceptance = metric(z$chi_square_accept & !z$accept))
    do.call(rbind, lapply(names(outcomes), function(nm)
      cbind(event = nm, as.data.frame(as.list(outcomes[[nm]])))))
  })
}

extract_comparisons <- function(panel) {
  records <- lapply(seq_len(nrow(panel$plan)), function(i) {
    p <- panel$plan[i, ]
    x <- panel$cases[[i]]
    z <- x$sets
    pairs <- list(coarse99_vs_fine99 = c("coarse99", "fine99"),
      coarse49_vs_fine49 = c("coarse49", "fine49"),
      fine49_vs_fine99 = c("fine49", "fine99"),
      coarse49_vs_coarse99 = c("coarse49", "coarse99"))
    result <- lapply(names(pairs), function(key) {
      a <- z[[pairs[[key]][1L]]]
      b <- z[[pairs[[key]][2L]]]
      data.frame(p, comparison = key, same_components = same_components(a$set, b$set),
        first_covers = or(a$covers, NA),
        second_covers = or(b$covers, NA),
        log_jaccard = union_overlap(a$set, b$set))
    })
    a <- z$fine99
    truth_accept <- if (is.null(x$point)) NA else x$point$accept[x$point$draws == 99L]
    result[[length(result) + 1L]] <- data.frame(p, comparison = "fine99_vs_direct99",
      same_components = NA, first_covers = or(a$covers, NA),
      second_covers = truth_accept, log_jaccard = NA_real_)
    rows <- do.call(rbind, result)
    rows$coverage_changed <- rows$first_covers != rows$second_covers
    rows
  })
  do.call(rbind, records)
}

extract_pair_audit <- function(panel) {
  extract_one <- function(x, target) {
    if (is.null(x)) return(NULL)
    g <- x[[paste0(target, "_result")]]
    pairs <- g$bootstrap_records
    if (is.null(pairs) || !nrow(pairs)) return(NULL)
    pairs <- pairs[!pairs$ok | nzchar(pairs$warnings), , drop = FALSE]
    if (!nrow(pairs)) return(NULL)
    data.frame(panel = panel$plan$panel[1L], id = x$case$id, target, pairs)
  }
  do.call(rbind, lapply(panel$cases, function(x)
    rbind(extract_one(x, "grid"), extract_one(x, "point"))))
}

grouped <- function(data, columns, fun) {
  group <- interaction(data[columns], drop = TRUE)
  do.call(rbind, lapply(split(data, group), function(x) {
    result <- fun(x)
    keys <- x[rep(1L, nrow(result)), columns, drop = FALSE]
    rownames(keys) <- rownames(result) <- NULL
    cbind(keys, result)
  }))
}

set_summary <- function(data, scope) {
  data$scope <- scope
  grouped(data, c("scope", "scenario", "method"), function(x) {
    y <- x[x$present & (scope == "all_attempted" | x$prefix), ]
    z <- y[y$available, ]
    data.frame(planned = nrow(x), attempted = nrow(y), available = nrow(z),
      coverage_sensitivity_min = if (nrow(y)) sum(z$covers) / nrow(y) else NA_real_,
      coverage_sensitivity_max = if (nrow(y)) (sum(z$covers) + nrow(y) - nrow(z)) / nrow(y) else NA_real_,
      certain_coverage_min = if (nrow(y)) sum(z$certain_covers) / nrow(y) else NA_real_,
      possible_coverage_max = if (nrow(y)) (sum(z$possible_covers) + nrow(y) - nrow(z)) / nrow(y) else NA_real_,
      any_ambiguity = sum(y$ambiguous_nodes > 0L, na.rm = TRUE),
      median_log_span = num(z$log_span, stats::median),
      median_profile_jaccard = num(z$profile_jaccard_log_length, stats::median),
      max_log_bracket = num(z$log_endpoint_bracket_max, max))
  })
}

binary_set <- function(data, scope) {
  data$scope <- scope
  grouped(data, c("scope", "scenario", "method"), function(x) {
    y <- x[x$present & (scope == "all_attempted" | x$prefix), ]
    z <- y[y$available, ]
    outcomes <- list(available = metric(y$available),
      unconditional_coverage = metric(y$covers), conditional_coverage = metric(z$covers),
      empty = metric(z$empty), zero_log_length = metric(z$log_span == 0),
      disconnected = metric(z$disconnected),
      singleton_component = metric(z$singleton_component), any_bound = metric(z$any_bound),
      full_range = metric(z$full_range))
    do.call(rbind, lapply(names(outcomes), function(nm)
      cbind(event = nm, as.data.frame(as.list(outcomes[[nm]])))))
  })
}

binary_point <- function(data, scope) {
  data$scope <- scope
  grouped(data, c("scope", "panel", "scenario", "draws"), function(x) {
    y <- x[x$present & (scope == "all_attempted" | x$prefix), ]
    z <- y[!is.na(y$accept), ]
    outcomes <- list(available = metric(!is.na(y$accept)),
      unconditional_acceptance = metric(y$accept), conditional_acceptance = metric(z$accept),
      ambiguous = metric(z$ambiguous))
    do.call(rbind, lapply(names(outcomes), function(nm)
      cbind(event = nm, as.data.frame(as.list(outcomes[[nm]])))))
  })
}

comparison_summary <- function(data, scope) {
  data$scope <- scope
  grouped(data, c("scope", "scenario", "comparison"), function(x) {
    y <- x[x$present & (scope == "all_attempted" | x$prefix), ]
    observed <- !is.na(y$coverage_changed)
    a <- y$first_covers[observed]
    b <- y$second_covers[observed]
    events <- list(coverage_changed = metric(a != b),
      first_only = metric(a & !b), second_only = metric(b & !a),
      same_components = metric(y$same_components, sum(!is.na(y$same_components))))
    do.call(rbind, lapply(names(events), function(nm)
      cbind(event = nm, as.data.frame(as.list(events[[nm]])),
        median_jaccard = num(y$log_jaccard, stats::median))))
  })
}

grid <- load_panel(args[1L], args[-seq_len(3L)])
point <- load_panel(args[2L])
stopifnot(identical(grid$configuration$source, point$configuration$source),
  identical(grid$configuration$script, point$configuration$script),
  identical(grid$configuration$B, point$configuration$B),
  identical(grid$configuration$level, point$configuration$level),
  !any(grid$plan$seed %in% point$plan$seed))
datasets <- rbind(extract_datasets(grid), extract_datasets(point))
sets <- extract_sets(grid)
points <- rbind(extract_points(grid), extract_points(point))
comparisons <- extract_comparisons(grid)
scopes <- c("balanced_prefix", "all_attempted")
bind <- function(f, x) do.call(rbind, lapply(scopes, function(scope) f(x, scope)))
tables <- list(datasets = datasets, sets = sets, points = points, paired_comparisons = comparisons,
  set_summary = bind(set_summary, sets), set_exact_binomial = bind(binary_set, sets),
  point_exact_binomial = bind(binary_point, points),
  comparison_summary = bind(comparison_summary, comparisons),
  point_chi_square_comparison = bind(point_comparison, points),
  returned_components = extract_components(grid), grid_nodes = extract_nodes(grid),
  pair_warnings_and_failures = rbind(extract_pair_audit(grid), extract_pair_audit(point)))
for (nm in names(tables))
  utils::write.csv(tables[[nm]], file.path(out, paste0(nm, ".csv")), row.names = FALSE)
print(with(datasets, table(panel, scenario, present, useNA = "ifany")))
