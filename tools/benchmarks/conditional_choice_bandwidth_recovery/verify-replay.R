args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1L] else "tools/benchmarks/conditional_choice_bandwidth_recovery/results"
read <- function(parts, name) do.call(rbind, lapply(file.path(path, parts, paste0(name, ".csv")),
  read.csv, stringsAsFactors = FALSE))
prior_dirs <- c("initial_grid_part1", "initial_grid_part2")
final_dirs <- c("part1", "part2")
sort_by <- function(x, cols) {
  x <- x[do.call(order, x[cols]), ]
  rownames(x) <- NULL
  x
}
for (nm in c("design", "assessment", "starts")) {
  key <- switch(nm, design = "id", assessment = c("id", "choice"), starts = c("id", "start"))
  a <- sort_by(read(prior_dirs, nm), key)
  b <- sort_by(read(final_dirs, nm), key)
  if (nm == "starts") {
    a$seconds <- b$seconds <- NULL
  }
  stopifnot(identical(a, b))
}
a <- sort_by(read(prior_dirs, "scenarios"), "id")
b <- sort_by(read(final_dirs, "scenarios"), "id")
fields <- c("id", grep("^fitted_|^training_|^eligible_starts$", names(a), value = TRUE))
stopifnot(identical(a[fields], b[fields]))
out <- data.frame(id = a$id, initial_bandwidth = a$independent_bandwidth,
  refined_bandwidth = b$independent_bandwidth, initial_nll = a$independent_nll,
  refined_nll = b$independent_nll, initial_constrained_score = a$reference_constrained_score,
  refined_constrained_score = b$reference_constrained_score)
write.csv(out, file.path(path, "independent-reference-refinement.csv"), row.names = FALSE)
cat("Public outputs, start diagnostics except timing, and held-out probabilities reproduced exactly.\n")
cat("Initial projected score max", max(abs(out$initial_constrained_score)),
  "; refined max", max(abs(out$refined_constrained_score)), "\n")
