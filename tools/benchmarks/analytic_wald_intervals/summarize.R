args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply the run directory and one output directory.")
input <- args[1L]
output <- args[2L]
dir.create(output, showWarnings = FALSE, recursive = TRUE)
dirs <- file.path(input, paste0("shard", 1:2))
configs <- lapply(file.path(dirs, "config.rds"), readRDS)
stopifnot(identical(configs[[1L]], configs[[2L]]))
config <- configs[[1L]]
plan <- config$plan
files <- file.path(dirs[ifelse(plan$replicate %% 2L == 1L, 1L, 2L)], paste0(plan$id, ".rds"))
available <- file.exists(files)
attempts <- unlist(lapply(files[available], readRDS), recursive = FALSE, use.names = FALSE)
extract <- function(name) do.call(rbind, lapply(attempts, function(x) x[[name]]))
metrics <- extract("metrics")
sets <- extract("sets")
package_profile <- extract("package_profile")
starts <- lapply(attempts, function(x) {
  if (is.null(x$starts)) return(NULL)
  cbind(id = x$metrics$id[1L], backend = x$metrics$backend[1L], x$starts)
})
start_fields <- unique(unlist(lapply(starts, names), use.names = FALSE))
starts <- do.call(rbind, lapply(starts, function(x) {
  if (is.null(x)) return(NULL)
  for (name in setdiff(start_fields, names(x))) x[[name]] <- NA
  x[start_fields]
}))
reference_warnings <- do.call(rbind, lapply(attempts, function(x) {
  if (is.null(x$failed_reference_points) || !nrow(x$failed_reference_points)) return(NULL)
  cbind(id = x$metrics$id[1L], backend = x$metrics$backend[1L], x$failed_reference_points)
}))
if (!is.null(package_profile)) package_profile$max_endpoint_gap_map <- with(package_profile,
  pmax(abs(lower - connected_lower), abs(upper - connected_upper)))
for (name in c("plan", "metrics", "sets", "package_profile"))
  utils::write.csv(get(name), file.path(output, paste0(name, ".csv")), row.names = FALSE)
for (name in c("starts", "reference_warnings")) {
  value <- get(name)
  if (is.null(value)) value <- data.frame(id = character(), backend = character())
  connection <- gzfile(file.path(output, paste0(name, ".csv.gz")), "wt")
  utils::write.csv(value, connection, row.names = FALSE)
  close(connection)
}
utils::write.csv(data.frame(path = names(config$sources), sha256 = unname(config$sources)),
  file.path(output, "source-sha256.csv"), row.names = FALSE)
writeLines(c(config$R, paste("Locale", config$locale), paste("MASS", config$MASS)),
  file.path(output, "runtime.txt"))
plan$completed <- available
utils::write.csv(plan, file.path(output, "completion.csv"), row.names = FALSE)
stopifnot(!anyDuplicated(metrics[c("id", "backend", "level")]),
  all(metrics$id %in% plan$id[available]), all(is.finite(metrics$n)),
  all(metrics$profile_cover[metrics$profile_available] ==
    metrics$profile_set_cover[metrics$profile_available]))
stopifnot(all(table(metrics$id) %in% c(2L, 4L)),
  identical(sort(unique(metrics$id)), sort(plan$id[available])))
expected <- 2L + 2L * (plan$family == "nbinom2")
stopifnot(nrow(metrics) == sum(expected[available]))

coverage <- function(covered, denominator, estimand) {
  count <- sum(covered %in% TRUE)
  z <- stats::qnorm(.975)
  p <- if (denominator) count / denominator else NA_real_
  center <- if (denominator) (p + z^2 / (2 * denominator)) / (1 + z^2 / denominator) else NA_real_
  half <- if (denominator) z * sqrt(p * (1 - p) / denominator + z^2 / (4 * denominator^2)) /
    (1 + z^2 / denominator) else NA_real_
  data.frame(estimand = estimand, covered = count, denominator = denominator,
    rate = p, lo95 = pmax(0, center - half), hi95 = pmin(1, center + half))
}
sum_known <- function(x) sum(x %in% TRUE)
describe <- function(d, group) {
  counts <- data.frame(n_done = nrow(d), fit = sum_known(d$fit_available),
    wald = sum_known(d$wald_available), wald_outside_search = sum_known(d$wald_outside_search),
    profile = sum_known(d$profile_available), scale_mle_boundary = sum_known(d$fit_scale_boundary),
    any_mle_boundary = sum_known(d$fit_boundary), poisson_face_mle = sum_known(d$fit_poisson_face),
    profile_bound = sum_known(d$profile_lower_bound | d$profile_upper_bound),
    profile_full = sum_known(d$profile_full), profile_disconnected = sum_known(d$profile_components > 1L),
    profile_no_effect_inside = sum_known(d$null_effect_lr <= stats::qchisq(d$level, 1)),
    wald_median_log_width = if (any(d$wald_available)) median(log(d$wald_upper[d$wald_available] /
      d$wald_lower[d$wald_available])) else NA_real_,
    profile_median_log_width = if (any(d$profile_available))
      median(d$profile_width_log[d$profile_available]) else NA_real_)
  both <- d$wald_available & d$profile_available
  matched <- data.frame(paired = sum(both),
    both_cover = sum_known(both & d$wald_cover & d$profile_cover),
    only_wald = sum_known(both & d$wald_cover & !d$profile_cover),
    only_profile = sum_known(both & !d$wald_cover & d$profile_cover),
    neither = sum_known(both & !d$wald_cover & !d$profile_cover),
    wald_missing_profile_covers = sum_known(!d$wald_available & d$profile_cover))
  rates <- rbind(
    coverage(d$wald_cover, sum(d$wald_available), "wald_conditional"),
    coverage(d$wald_cover, nrow(d), "wald_unconditional"),
    coverage(d$profile_cover, sum(d$profile_available), "profile_conditional"),
    coverage(d$profile_cover, nrow(d), "profile_unconditional"),
    coverage(d$profile_cover & d$wald_available, sum(both), "profile_matched_to_wald"))
  list(counts = cbind(group, counts, matched), rates = cbind(group[rep(1L, nrow(rates)), ], rates))
}
for (type in c("cell", "pooled")) {
  keys <- if (type == "cell") c("family", "backend", "scenario", "n", "level") else
    c("family", "backend", "scenario", "level")
  groups <- split(metrics, do.call(paste, c(metrics[keys], sep = "/")))
  result <- lapply(groups, function(d) describe(d, d[1L, keys, drop = FALSE]))
  for (item in c("counts", "rates")) utils::write.csv(do.call(rbind, lapply(result, function(x) x[[item]])),
    file.path(output, paste0(type, "-", item, ".csv")), row.names = FALSE)
}

# Descriptive, outcome-blind warning flags. Evaluate what users can learn from
# the reported interval geometry without applying a new eligibility rule.
geometry <- metrics[metrics$wald_available & metrics$profile_available, , drop = FALSE]
geometry$wald_geometry <- ifelse(geometry$wald_outside_search, "crosses_bound", "inside_search")
keys <- c("family", "backend", "scenario", "n", "level", "wald_geometry")
grouped_geometry <- split(geometry, do.call(paste, c(geometry[keys], sep = "/")))
geometry_summary <- do.call(rbind, lapply(grouped_geometry, function(d) {
  counts <- data.frame(n = nrow(d), wald_covers = sum_known(d$wald_cover),
    profile_covers = sum_known(d$profile_cover),
    profile_bound = sum_known(d$profile_lower_bound | d$profile_upper_bound),
    profile_full = sum_known(d$profile_full), profile_disconnected = sum_known(d$profile_components > 1L))
  cbind(d[1L, keys, drop = FALSE], counts)
}))
utils::write.csv(geometry_summary, file.path(output, "wald-geometry.csv"), row.names = FALSE)

# Separate statistical shape from convergence. Compare endpoint-crossing Wald
# flags to the union's topology, and retain fit/profile failure denominators.
scenarios <- split(metrics, interaction(metrics$scenario, metrics$n, metrics$level, drop = TRUE))
scenario_counts <- do.call(rbind, lapply(scenarios, function(d) {
  counts <- data.frame(n = nrow(d), wald = sum_known(d$wald_available),
    wald_covers = sum_known(d$wald_cover), profile = sum_known(d$profile_available),
    profile_covers = sum_known(d$profile_cover), profile_bound = sum_known(d$profile_lower_bound | d$profile_upper_bound),
    profile_full = sum_known(d$profile_full), profile_disconnected = sum_known(d$profile_components > 1L))
  cbind(d[1L, c("scenario", "n", "level"), drop = FALSE], stats::setNames(counts, sub("^n$", "backend_attempts", names(counts))))
}))
utils::write.csv(scenario_counts, file.path(output, "scenario-by-n.csv"), row.names = FALSE)
# Paired backend contrasts retain the same generating NB2 dataset.
nb <- metrics[metrics$family == "nbinom2", , drop = FALSE]
profile <- nb[nb$backend == "profile_glm", ]
joint <- nb[nb$backend == "joint_nb2", ]
stopifnot(identical(profile$id, joint$id), identical(profile$level, joint$level))
paired_nb <- data.frame(profile[c("id", "scenario", "n", "replicate", "level")],
  profile_glm_wald = profile$wald_available, joint_nb2_wald = joint$wald_available,
  profile_glm_wald_cover = profile$wald_cover, joint_nb2_wald_cover = joint$wald_cover,
  profile_glm_profile_cover = profile$profile_cover, joint_nb2_profile_cover = joint$profile_cover,
  joint_exact_poisson = joint$fit_poisson_face,
  profile_glm_theta = profile$fit_theta, joint_nb2_theta = joint$fit_theta)
utils::write.csv(paired_nb, file.path(output, "nb2-matched-backends.csv"), row.names = FALSE)


m <- metrics[metrics$level == .95, ]
checks <- do.call(rbind, lapply(split(m, paste(m$family, m$backend)), function(x) {
  errors <- abs(x$same_scale_nll_gap[is.finite(x$same_scale_nll_gap)])
  gaps <- abs(x$global_nll_gap[is.finite(x$global_nll_gap)])
  data.frame(family = x$family[1L], backend = x$backend[1L], attempted = nrow(x),
    fit = sum(x$fit_available), profile = sum(x$profile_available),
    wald = sum(x$wald_available), Wald_covers = sum_known(x$wald_cover),
    profile_covers = sum_known(x$profile_cover),
    same_scale_median = if (length(errors)) median(errors) else NA_real_,
    same_scale_q99 = if (length(errors)) quantile(errors, .99) else NA_real_,
    same_scale_max = if (length(errors)) max(errors) else NA_real_,
    same_scale_above_1e5 = sum(errors > 1e-5),
    global_median = if (length(gaps)) median(gaps) else NA_real_,
    global_q99 = if (length(gaps)) quantile(gaps, .99) else NA_real_,
    global_max = if (length(gaps)) max(gaps) else NA_real_,
    global_above_1e5 = sum(gaps > 1e-5))
}))
utils::write.csv(checks, file.path(output, "likelihood-checks.csv"), row.names = FALSE)
print(checks, row.names = FALSE)
cat("Completed datasets", sum(available), "of", length(available), "\n")
cat("Mismatched coarse/fine component counts", sum_known(m$profile_refined_count_change), "\n")
cat("Maximum coarse/fine log endpoint discrepancy", max(m$profile_refined_endpoint_max, na.rm = TRUE), "\n")
if (!is.null(package_profile)) {
  cat("Public package profile checks returned", sum(package_profile$returned), "of", nrow(package_profile), "\n")
  cat("Public profile maximum finite endpoint discrepancy", max(package_profile$max_endpoint_gap_map, na.rm = TRUE), "\n")
}
