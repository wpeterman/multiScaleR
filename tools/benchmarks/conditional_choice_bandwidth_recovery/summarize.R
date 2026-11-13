args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1L] else "tools/benchmarks/conditional_choice_bandwidth_recovery/results"
dirs <- file.path(root, c("part1", "part2"))
read <- function(nm) do.call(rbind, lapply(file.path(dirs, paste0(nm, ".csv")),
  read.csv, stringsAsFactors = FALSE))
cases <- read("scenarios")
starts <- read("starts")
predictions <- read("assessment")
profiles <- read("profiles")
refined <- read("refined")
planned <- read("design")
stopifnot(nrow(cases) == 384L, setequal(cases$id, planned$id), !anyDuplicated(cases$id),
  !anyDuplicated(starts[c("id", "start")]), !anyDuplicated(predictions[c("id", "choice")]))
cases <- cases[order(cases$id), ]
stopifnot(identical(as.integer(table(factor(predictions$id, levels = cases$id))), rep(160L, nrow(cases))))
cases$paired <- with(cases, fitted_status == "scored" & reference_status == "returned")
cases$objective_gap <- with(cases, fitted_nll - independent_nll)
cases$profile_support_ratio <- with(cases, reference_support_upper/reference_support_lower)
cases$bandwidth_ratio <- with(cases, fitted_bandwidth/true_bandwidth)
cases$reference_bandwidth_ratio <- with(cases, independent_bandwidth/true_bandwidth)
cases$boundary_side <- with(cases, ifelse(is.na(fitted_bandwidth), NA_character_,
  ifelse(fitted_bandwidth <= 60 * exp(1e-5), "lower",
    ifelse(fitted_bandwidth >= 800 * exp(-1e-5), "upper", "interior"))))
stopifnot(all(cases$training_K[cases$fitted_status == "scored"] == 4L),
  all(cases$training_nobs[cases$fitted_status == "scored"] == cases$n_analysis[cases$fitted_status == "scored"]),
  !any(cases$fitted_scale_se_available), all(cases$fitted_boundary[cases$fitted_status == "scored"] ==
    (cases$boundary_side[cases$fitted_status == "scored"] != "interior")))
sets <- merge(predictions, cases[c("id", "fitted_status")], by = "id", sort = FALSE)
scored <- subset(sets, fitted_status == "scored")
p <- as.matrix(scored[c("p1", "p2", "p3", "p4")])
stopifnot(all(p >= 0 & p <= 1), max(abs(rowSums(p) - 1)) < 1e-12)
observed_loss <- rowsum(-log(p[cbind(seq_len(nrow(p)), scored$used)]), scored$id)
out <- cases[cases$fitted_status == "scored", ]
stopifnot(isTRUE(all.equal(as.numeric(observed_loss[match(out$id, rownames(observed_loss)), ]),
  out$fitted_assessment_nll, tolerance = 1e-10)))
groups <- split(cases, list(cases$kernel, cases$n_analysis, cases$signal), drop = TRUE)
set.seed(25120)
ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L) return(c(estimate = if (n) x else NA_real_, lower = NA_real_, upper = NA_real_))
  boot <- vapply(seq_len(4000L), function(i) mean(x[sample.int(n, replace = TRUE)]), 0)
  c(estimate = mean(x), lower = unname(quantile(boot, .025)), upper = unname(quantile(boot, .975)))
}
qs <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) rep(NA_real_, 3L) else unname(quantile(x, c(.1, .5, .9)))
}
key <- function(x) x[1L, c("kernel", "n_analysis", "signal")]
recovery <- do.call(rbind, lapply(groups, function(x) {
  p <- x[x$paired, ]
  q <- qs(p$bandwidth_ratio)
  ref <- qs(x$reference_bandwidth_ratio)
  span <- qs(x$profile_support_ratio)
  cbind(key(x), data.frame(attempted = nrow(x), public_scored = sum(x$fitted_status == "scored"),
    reference_returned = sum(x$reference_status == "returned"), paired = nrow(p),
    public_lower = sum(p$boundary_side == "lower"), public_upper = sum(p$boundary_side == "upper"),
    reference_lower = sum(x$independent_bandwidth <= 60 * exp(1e-5), na.rm = TRUE),
    reference_upper = sum(x$independent_bandwidth >= 800 * exp(-1e-5), na.rm = TRUE),
    fitted_ratio_p10 = q[1L], fitted_ratio_median = q[2L], fitted_ratio_p90 = q[3L],
    reference_ratio_p10 = ref[1L], reference_ratio_median = ref[2L], reference_ratio_p90 = ref[3L],
    recovered_within_25pct = if (all(is.na(p$bandwidth_ratio))) NA_integer_ else sum(abs(p$bandwidth_ratio - 1) <= .25),
    median_support_upper_to_lower = span[2L],
    support_reaches_both_bounds = sum(x$reference_support_lower <= 60 * exp(1e-5) &
      x$reference_support_upper >= 800 * exp(-1e-5), na.rm = TRUE),
    disconnected_support = sum(x$reference_support_regions > 1L, na.rm = TRUE),
    better_reference_cases = sum(p$objective_gap > .0001),
    max_reference_objective_advantage = if (nrow(p)) max(p$objective_gap) else NA_real_))
}))
prediction <- do.call(rbind, lapply(groups, function(x) {
  p <- x[x$paired, ]
  do.call(rbind, lapply(c("fixed_truth", "oracle", "no_habitat", "independent"), function(method) {
    loss <- p[[paste0(method, "_assessment_nll")]]
    cbind(key(x), data.frame(comparator = method, attempted = nrow(x), paired = nrow(p),
      fitted_mean_nll = mean(p$fitted_assessment_nll/160),
      comparator_mean_nll = mean(loss/160)),
      as.data.frame(as.list(ci((p$fitted_assessment_nll - loss)/160))))
  }))
}))
agreement <- do.call(rbind, lapply(groups, function(x) {
  p <- x[x$paired, ]
  metric <- list(public_vs_own_profile = abs(p$public_vs_own_profile_nll),
    probability_gap = p$max_probability_gap, raw_coefficient_gap = p$max_raw_coefficient_gap,
    independent_projected_bandwidth_score = abs(p$reference_constrained_score),
    independent_nuisance_score = p$reference_nuisance_score)
  do.call(rbind, lapply(names(metric), function(nm) {
    m <- metric[[nm]]
    q <- if (length(m)) unname(quantile(m, c(.5, .95, 1))) else rep(NA_real_, 3L)
    cbind(key(x), data.frame(metric = nm, paired = nrow(p), median = q[1L], p95 = q[2L], max = q[3L]))
  }))
}))
flagged <- subset(cases, reference_status == "returned" & (fitted_status == "failed" | objective_gap > .0001))
for (nm in c("cases", "starts", "recovery", "prediction", "agreement", "flagged"))
  write.csv(get(nm), file.path(root, paste0(nm, "-combined.csv")), row.names = FALSE)
write.csv(data.frame(id = names(table(profiles$id)), grid_points = as.integer(table(profiles$id))),
  file.path(root, "profile-counts.csv"), row.names = FALSE)
write.csv(data.frame(id = names(table(refined$id)), candidates = as.integer(table(refined$id))),
  file.path(root, "refined-counts.csv"), row.names = FALSE)
print(recovery, row.names = FALSE)
print(subset(prediction, comparator %in% c("fixed_truth", "no_habitat")), row.names = FALSE)
cat("Total scored/independent reference\n")
print(with(cases, table(fitted_status, reference_status)))
cat("Boundary-start diagnostic triggers:", nrow(flagged), "\n")
