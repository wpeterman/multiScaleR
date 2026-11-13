# Independently verify the compact fit and paired assessment records.
root <- "tools/benchmarks/bayesian_scales/uncorrected_laplace_20261112"
archive <- file.path(root, "results")
paths <- list.files(archive, pattern = "[.]rds$", full.names = TRUE)
configs <- lapply(paths[grepl("configuration-", paths, fixed = TRUE)], readRDS)
stopifnot(length(configs) == 2L, identical(configs[[1]]$hash, configs[[2]]$hash),
  identical(configs[[1]]$plan, configs[[2]]$plan))
expected <- configs[[1]]$hash
stopifnot(all(vapply(names(expected), function(f)
  digest::digest(file = f, algo = "sha256"), character(1)) == expected))
records <- lapply(paths[!grepl("configuration-", paths, fixed = TRUE)], readRDS)
keys <- vapply(records, function(x) paste(x$case, x$method), character(1))
stopifnot(length(keys) == 36L, !anyDuplicated(keys),
  all(vapply(records, function(x) identical(x$hash, expected), logical(1))))
fits <- read.csv(file.path(archive, "fits.csv"))
stopifnot(nrow(fits) == 36L, setequal(paste(fits$case, fits$method), keys),
  all(fits$status %in% c("eligible", "rejected", "failed")))
for (r in records) {
  row <- fits[fits$case == r$case & fits$method == r$method, , drop = FALSE]
  stopifnot(identical(r$status, row$status))
  if (is.null(r$intervals)) next
  stopifnot(all(r$intervals$lower <= r$intervals$median),
    all(r$intervals$median <= r$intervals$upper),
    all(r$intervals$covered == with(r$intervals, lower <= truth & truth <= upper)),
    all(r$intervals$cdf_at_truth >= 0 & r$intervals$cdf_at_truth <= 1),
    isTRUE(all.equal(r$joint_log_score, row$joint_log_score, tolerance = 1e-12)))
}
p <- read.csv(file.path(archive, "paired_probabilities.csv"))
q <- read.csv(file.path(archive, "paired_intervals.csv"))
stopifnot(all(abs(abs(p$estimate - p$estimate_nuts) - p$absolute_difference) < 1e-12),
  identical(p$material, p$absolute_difference > pmax(.05, 5 * p$combined_mcse)),
  all(abs(q$width_ratio - (q$upper - q$lower) / (q$upper_nuts - q$lower_nuts)) < 1e-12))
cat(sprintf("Verified %d fit records, %d paired scale probabilities, and %d paired intervals.\n",
  nrow(fits), nrow(p), nrow(q)))
