# Review completed and unavailable draws without rerunning or dropping fits.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
output <- args[1L]
read <- function(name) utils::read.csv(file.path(output, paste0(name, ".csv")))
fits <- read("fits")
eligibility <- read("eligibility")
stopifnot(!anyDuplicated(fits[c("case", "method")]))
print(eligibility[eligibility$planned > eligibility$unrun, ], row.names = FALSE)
if (!file.exists(file.path(output, "ranks.csv"))) quit(save = "no")
ranks <- read("ranks")
stopifnot(all(ranks$cdf_at_truth >= -1e-12 & ranks$cdf_at_truth <= 1 + 1e-12))
ranks$cdf_at_truth <- pmin(1, pmax(0, ranks$cdf_at_truth))
split_columns <- function(data, fields) split(data,
  do.call(interaction, c(unname(as.list(data[fields])), list(drop = TRUE))))
valid <- ranks[ranks$status == "eligible", , drop = FALSE]
rank_cells <- split_columns(valid, c("family", "design", "method", "quantity"))
calibration <- do.call(rbind, lapply(rank_cells, function(x) {
  p <- sort(x$cdf_at_truth)
  n <- length(p)
  D <- max(seq_len(n) / n - p, p - (seq_len(n) - 1L) / n)
  # Simultaneous band for three monitored CDFs within a design-method cell.
  band <- sqrt(log(2 * 3 / .05) / (2 * n))
  coverage <- binom.test(sum(p >= .05 & p <= .95), n)$conf.int
  cbind(x[1L, c("family", "design", "method", "quantity")],
    n = n, coverage_90 = mean(p >= .05 & p <= .95), lower_95 = coverage[1L],
    upper_95 = coverage[2L], ecdf_deviation = D, ecdf_band = band,
    ecdf_flag = D > band)
}))
if (nrow(as.data.frame(calibration)))
  write.csv(calibration, file.path(output, "rank-calibration.csv"), row.names = FALSE)
if (!file.exists(file.path(output, "probabilities.csv"))) quit(save = "no")
probability <- read("probabilities")
p <- probability[probability$reference_status == "eligible" &
  probability$importance_status == "eligible", , drop = FALSE]
p$family <- fits$family[match(p$case, fits$case)]
p$design <- fits$design[match(p$case, fits$case)]
p_cells <- split_columns(p, c("family", "design", "method"))
paired <- do.call(rbind, lapply(p_cells, function(x) {
  usable <- is.finite(x$combined_mcse)
  key <- x[1L, c("family", "design", "method")]
  data <- fits[fits$case %in% unique(x$case) & fits$status == "eligible", ]
  ref <- data[data$method == "nuts", ]
  alt <- data[data$method == key$method, ]
  loss <- alt$joint_log_score - ref$joint_log_score[match(alt$case, ref$case)]
  cbind(key, common_datasets = length(unique(x$case)), usable_events = sum(usable),
    unavailable_events = sum(!usable), material_differences = sum(x$material_difference[usable]),
    median_absolute_probability_difference = median(x$absolute_difference[usable]),
    max_absolute_probability_difference = if (any(usable)) max(x$absolute_difference[usable]) else NA,
    median_joint_log_loss_difference = median(loss),
    max_abs_joint_log_loss_difference = max(abs(loss)))
}))
if (nrow(as.data.frame(paired)))
  write.csv(paired, file.path(output, "paired-assessment.csv"), row.names = FALSE)
print(paired, row.names = FALSE)

# Plot only complete qualification panels. The study unit is the dataset.
available <- eligibility[eligibility$availability_qualified, , drop = FALSE]
if (nrow(available)) {
  panels <- split_columns(available, c("family", "design"))
  colors <- c(nuts = "#222222", joint_importance = "#1b7494", conditional_beta = "#b45c20")
  labels <- c(nuts = "NUTS", joint_importance = "Joint proposal",
    conditional_beta = "Conditional coefficients")
  quantities <- c(beta = "Habitat coefficient", log_distance = "First-year log distance",
    scale_sd = "Between-year scale SD")
  for (panel in panels) {
    family <- panel$family[1L]
    design <- panel$design[1L]
    keep <- valid$family == family & valid$design == design & valid$method %in% panel$method
    figure <- file.path(output, paste0("ranks-", family, "-", design, ".png"))
    grDevices::png(figure, width = 1550, height = 570, res = 150)
    op <- graphics::par(mfrow = c(1, 3), mar = c(4.2, 4.3, 2.5, 0.8))
    for (quantity in names(quantities)) {
      sel <- valid[keep & valid$quantity == quantity, ]
      band <- max(calibration$ecdf_band[calibration$family == family &
        calibration$design == design & calibration$quantity == quantity &
        calibration$method %in% panel$method])
      graphics::plot(NA_real_, xlim = c(0, 1), ylim = c(0, 1), asp = 1,
        main = quantities[[quantity]], xlab = "Posterior CDF at truth", ylab = "Empirical CDF")
      x <- seq(0, 1, length.out = 101)
      graphics::polygon(c(x, rev(x)), c(pmin(1, x + band), rev(pmax(0, x - band))),
        col = "#e9edf0", border = NA)
      graphics::abline(0, 1, col = "#767676", lty = 2)
      for (method in panel$method) {
        p <- sel$cdf_at_truth[sel$method == method]
        graphics::lines(stats::ecdf(p), col = colors[[method]], lwd = 2, do.points = FALSE)
      }
      if (quantity == "beta")
        graphics::legend("bottomright", unname(labels[panel$method]),
          col = unname(colors[panel$method]), lwd = 2, bty = "n", cex = .72)
    }
    graphics::par(op)
    grDevices::dev.off()
  }
}
