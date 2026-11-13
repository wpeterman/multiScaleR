# Render a report from the tables created by summarize.R.
# Rscript grid_report.R SUMMARY_DIRECTORY OUTPUT_MD
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- args[1L]
read <- function(file) utils::read.csv(file.path(root, paste0(file, ".csv")),
  stringsAsFactors = FALSE)
datasets <- read("datasets")
sets <- read("sets")
points <- read("points")
binary <- read("set_exact_binomial")
point_binary <- read("point_exact_binomial")
set_summary <- read("set_summary")
paired <- read("comparison_summary")
nodes <- read("grid_nodes")
scenario <- c("interior", "weak", "near_bound")
labels <- c("Interior", "Weak effect", "Near upper bound, Poisson limit")
methods <- c("fine99", "coarse99", "fine49", "coarse49")
scope <- "balanced_prefix"

bin_row <- function(s, m, event_name) {
  row <- binary[binary$scope == scope & binary$scenario == s & binary$method == m &
    binary$event == event_name, ]
  stopifnot(nrow(row) == 1L)
  row
}
set_row <- function(s, m) {
  row <- set_summary[set_summary$scope == scope & set_summary$scenario == s &
    set_summary$method == m, ]
  stopifnot(nrow(row) == 1L)
  row
}
rate <- function(row) sprintf("%d/%d", row$successes, row$denominator)
ci <- function(row) if (row$denominator) sprintf("%.3f to %.3f", row$lower95, row$upper95) else "NA"
dec <- function(x) if (is.finite(x)) sprintf("%.3f", x) else "NA"
line <- function(...) paste0("| ", paste(..., sep = " | "), " |")

grid <- datasets[datasets$panel == "grid", ]
full_grid <- grid[grid$present, ]
balanced <- grid[grid$prefix, ]
direct <- points[points$panel == "grid" & points$prefix & points$draws == 99L, ]
balanced_sets <- sets[sets$prefix & sets$method == "fine99", ]
stopifnot(nrow(grid) == 48L, !anyDuplicated(grid$id), nrow(direct) == nrow(balanced),
  nrow(balanced_sets) == nrow(balanced),
  all(full_grid$status %in% c("initial_fit_failed", "point_fit_failed",
    "point_inference_failed", "point_pair_failed", "point_complete")))

coverage <- unlist(lapply(seq_along(scenario), function(i) {
  s <- scenario[i]
  row <- bin_row(s, "fine99", "conditional_coverage")
  op <- bin_row(s, "fine99", "unconditional_coverage")
  av <- bin_row(s, "fine99", "available")
  truth <- point_binary[point_binary$scope == scope & point_binary$panel == "grid" &
    point_binary$scenario == s & point_binary$draws == 99L &
    point_binary$event == "conditional_acceptance", ]
  line(labels[i], rate(av), paste(rate(op), paste0("(", ci(op), ")")),
    paste(rate(row), paste0("(", ci(row), ")")),
    paste(rate(truth), paste0("(", ci(truth), ")")))
}))
shape <- unlist(lapply(seq_along(scenario), function(i) {
  s <- scenario[i]
  z <- set_row(s, "fine99")
  line(labels[i], rate(bin_row(s, "fine99", "full_range")),
    rate(bin_row(s, "fine99", "any_bound")),
    rate(bin_row(s, "fine99", "disconnected")),
    rate(bin_row(s, "fine99", "singleton_component")),
    rate(bin_row(s, "fine99", "zero_log_length")), dec(z$median_log_span), dec(z$median_profile_jaccard), dec(z$max_log_bracket))
}))
variants <- unlist(lapply(seq_along(scenario), function(i) {
  s <- scenario[i]
  z <- set_row(s, "fine99")
  rates <- vapply(methods, function(m) rate(bin_row(s, m, "unconditional_coverage")), "")
  line(labels[i], rates[1L], rates[2L], rates[3L], rates[4L], z$any_ambiguity,
    paste(dec(z$certain_coverage_min), "to", dec(z$possible_coverage_max)))
}))
contrasts <- c("coarse99_vs_fine99", "fine49_vs_fine99", "coarse49_vs_coarse99",
  "fine99_vs_direct99")
paired_rows <- unlist(lapply(contrasts, function(key) {
  z <- paired[paired$scope == scope & paired$comparison == key, ]
  values <- vapply(scenario, function(s) {
    one <- z[z$scenario == s, ]
    row <- one[one$event == "coverage_changed", ]
    direction <- one$successes[match(c("first_only", "second_only"), one$event)]
    sprintf("%s (%s), %d / %d", rate(row), ci(row), direction[1L], direction[2L])
  }, "")
  line(key, values[1L], values[2L], values[3L])
}))

n_different <- sum(vapply(scenario, function(s) {
  row <- paired[paired$scope == scope & paired$scenario == s &
    paired$comparison == "fine99_vs_direct99" & paired$event == "coverage_changed", ]
  row$successes
}, 0))
node99 <- nodes[nodes$id %in% balanced$id, ]
zero <- sets[sets$prefix & sets$available & sets$log_span == 0, c("id", "method")]
zero_detail <- if (!nrow(zero)) "All available sets had positive total log length." else {
  by_method <- split(zero$id, factor(zero$method, levels = methods), drop = TRUE)
  c("Available zero-length sets by method", "",
    vapply(names(by_method), function(m)
      paste0("- `", m, "` had ", paste(by_method[[m]], collapse = ", "), "."), ""))
}
failed <- sum(full_grid$grid_failed_pairs)
point_failed <- sum(full_grid$point_failed_pairs)
worked <- read("returned_components")
worked <- worked[worked$id == "near_bound-02" & worked$method == "fine99", ]
worked_text <- if (nrow(worked) == 1L) c("", sprintf(paste(
  "For `near_bound-02`, the accepted component was [%g, %g] m.",
  "The adjacent rejected node was %g m.",
  "This singleton contains only %g m, so the 730 m generating scale was outside it."),
  worked$lower, worked$upper, worked$lower_outside, worked$lower)) else character()
interior <- balanced_sets[balanced_sets$scenario == "interior" & balanced_sets$available, ]
interior_nodes <- node99[node99$id %in% interior$id, ]
interior_component <- read("returned_components")
interior_component <- interior_component[interior_component$id %in% interior$id &
  interior_component$method == "fine99", ]
geometry <- if (nrow(interior)) {
  count_accepted <- table(vapply(split(interior_nodes, interior_nodes$id),
    function(x) sum(x$bootstrap_accept), integer(1L)))
  neighbors <- vapply(split(interior_nodes, interior_nodes$id), function(x)
    sum(x$bootstrap_accept & x$sigma %in% c(200, 280)), integer(1L))
  misses <- interior_component[interior_component$id %in%
    interior$id[!interior$covers], ]
  bounds <- unique(interior_component$id[interior_component$lower <= 240 &
    interior_component$upper >= 240])
  stopifnot(setequal(interior$id[interior$covers], bounds))
  miss_below <- sum(vapply(split(misses, misses$id), function(x) all(x$upper < 240),
    logical(1L)))
  miss_above <- sum(vapply(split(misses, misses$id), function(x) all(x$lower > 240),
    logical(1L)))
  c("", sprintf(paste("Interior sets had these accepted-node counts by dataset:",
    "%s. %d of %d sets had exactly one component."), paste(paste0(
      names(count_accepted), " nodes in ", as.integer(count_accepted), " datasets"),
      collapse = ", "), sum(interior$components == 1L), nrow(interior)),
    sprintf(paste("Both fixed nodes flanking 240 m (200 and 280 m) were accepted in",
      "%d covering and %d noncovering datasets. Among %d noncovering sets, %d lay",
      "entirely below 240 m and %d entirely above it."), sum(neighbors[bounds] == 2L),
      sum(neighbors[setdiff(interior$id, bounds)] == 2L),
      sum(!interior$covers), miss_below, miss_above))
} else character()
text <- c("# NB2 bootstrap finite-grid confidence-set results",
  "",
  sprintf("The grid panel completed %d of %d planned datasets. The balanced prefix contains %d datasets.",
    nrow(full_grid), nrow(grid), nrow(balanced)),
  sprintf(paste("It attempted %s grid refit pairs (%d failed, %d warned) and %s",
    "separately seeded truth-test pairs (%d failed, %d warned)."),
    format(sum(full_grid$grid_pairs), big.mark = ","), failed, sum(full_grid$grid_warned_pairs),
    format(sum(full_grid$point_pairs), big.mark = ","), point_failed,
    sum(full_grid$point_warned_pairs)),
  sprintf("The balanced prefix contains %d ambiguous 99-draw node decisions among %d tested nodes.",
    sum(node99$mc_ambiguous, na.rm = TRUE), sum(!is.na(node99$bootstrap_accept))),
  "",
  "## Coverage and availability",
  "",
  "Coverage means membership of the generating scale in the literal union of returned closed intervals.",
  "The unconditional rate counts an unavailable set as no returned covering set.",
  "The conditional rate uses only available sets.",
  "Parentheses give exact 95% binomial intervals across datasets.",
  "The direct test column gives pointwise inclusion among datasets with an available direct decision.",
  "",
  "| Scenario | Available | Unconditional coverage | Conditional coverage | Direct truth inclusion |",
  "|:--|--:|--:|--:|--:|", coverage,
  "",
  sprintf("The fine grid and direct truth test disagreed in %d balanced-prefix datasets with both decisions available.",
    n_different),
  "The generating scale was off the fixed grid. Its direct test also used a separately seeded bootstrap.",
  "The comparison reflects both differences in the target and Monte Carlo uncertainty.",
  "",
  "## Localization",
  "",
  "Lengths are measured on the log-scale domain. Span is a fraction of the full search length.",
  "Jaccard is descriptive overlap with the root-refined chi-squared profile set.",
  "The latter is a likelihood comparator. It supplies no calibration truth.",
  "",
  "| Scenario | Full range | Bound contact | Disconnected | Singleton component | Zero total length | Median span | Median profile overlap | Maximum log bracket width |",
  "|:--|--:|--:|--:|--:|--:|--:|--:|--:|", shape,
  "", zero_detail, worked_text, geometry,
  "",
  "## Grid and Monte Carlo sensitivity",
  "",
  "The coarse projection removes nodes without refitting. Its decisions at shared nodes are identical to those of the fine run.",
  "The 49-draw result uses the first 49 bootstrap draws at each node. It is paired with the 99-draw result.",
  "The last column bounds observed coverage when ambiguous node decisions and unavailable sets can vary.",
  "It is a sensitivity range for these datasets.",
  "",
  "| Scenario | Fine 99 | Coarse 99 | Fine 49 | Coarse 49 | Fine 99 cases with ambiguity | Fine 99 sensitivity range |",
  "|:--|--:|--:|--:|--:|--:|--:|", variants,
  "",
  "The paired table reports the fraction of datasets whose coverage event changes. Parentheses give exact 95% intervals.",
  "The final two counts in each cell give datasets covered only by the first and only by the second method, respectively.",
  "For the direct test, the paired event is truth inclusion.",
  "",
  "| Comparison | Interior | Weak effect | Near bound |",
  "|:--|--:|--:|--:|", paired_rows,
  "",
  "Calibration remains experimental. The evidence is specific to the declared grids, a single Gaussian kernel predictor, and a finite bootstrap budget.",
  "See the [separate pointwise panel](POINTWISE_RESULTS.md) for the 96-dataset truth-test results.")
writeLines(text, args[2L])
