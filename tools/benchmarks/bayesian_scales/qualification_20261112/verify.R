# Verify completed qualifications from the archived case-level summaries.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
path <- args[[1L]]
read <- function(name) utils::read.csv(file.path(path, paste0(name, ".csv")))
fits <- read("fits")
ranks <- read("ranks")
eligibility <- read("eligibility")
calibration <- read("rank-calibration")
probability <- read("probabilities")
paired <- read("paired-assessment")
stopifnot(nrow(fits) == 1152L, all(eligibility$planned == 48L),
  !anyDuplicated(fits[c("case", "method")]),
  all(fits$status %in% c("eligible", "rejected", "failed", "unrun")))
key <- function(x) paste(x$family, x$design, x$method)
for (i in seq_len(nrow(eligibility))) {
  x <- eligibility[i, ]
  observed <- fits[key(fits) == key(x), , drop = FALSE]
  counts <- tabulate(match(observed$status, c("eligible", "rejected", "failed", "unrun")), 4L)
  stopifnot(identical(counts, as.integer(x[c("eligible", "rejected", "failed", "unrun")])))
  lower <- if (x$eligible) stats::qbeta(.05, x$eligible, x$planned - x$eligible + 1L) else 0
  stopifnot(isTRUE(all.equal(lower, x$lower_95)),
    x$availability_qualified == (x$unrun == 0L && lower >= .90))
}
for (i in seq_len(nrow(calibration))) {
  x <- calibration[i, ]
  observed <- ranks[key(ranks) == key(x) & ranks$quantity == x$quantity &
    ranks$status == "eligible", "cdf_at_truth"]
  covered <- sum(observed >= .05 & observed <= .95)
  stopifnot(x$n == length(observed), isTRUE(all.equal(x$coverage_90, covered / x$n)),
    isTRUE(all.equal(c(x$lower_95, x$upper_95),
      as.numeric(stats::binom.test(covered, x$n)$conf.int))))
}
stopifnot(all(is.finite(fits$joint_log_score[fits$status == "eligible"])))
for (i in seq_len(nrow(paired))) {
  x <- paired[i, ]
  prefix <- paste(x$family, x$design, "", sep = "-")
  p <- probability[startsWith(probability$case, prefix) & probability$method == x$method &
    probability$reference_status == "eligible" & probability$importance_status == "eligible", ]
  usable <- is.finite(p$combined_mcse)
  gap <- abs(p$estimate_nuts[usable] - p$estimate_importance[usable])
  material <- sum(gap > pmax(.05, 5 * p$combined_mcse[usable]))
  stopifnot(x$common_datasets == length(unique(p$case)),
    x$usable_events == sum(usable), x$unavailable_events == sum(!usable),
    x$material_differences == material)
}
cat(nrow(fits), "planned records,", nrow(calibration), "calibration rows,",
  nrow(probability), "paired probability records verified.\n")

