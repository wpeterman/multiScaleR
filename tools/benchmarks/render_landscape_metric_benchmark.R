args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    paste(
      "Usage: Rscript render_landscape_metric_benchmark.R",
      "<raw_csv> <summary_csv> <output_dir>"
    ),
    call. = FALSE
  )
}

raw_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
summary_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- args[[3]]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw <- utils::read.csv(raw_path, stringsAsFactors = FALSE)
summary_df <- utils::read.csv(summary_path, stringsAsFactors = FALSE)

composition_metrics <- c(
  "shdi", "shei", "sidi", "siei", "msidi", "msiei", "pr", "prd", "rpr", "ta"
)
edge_metrics <- c("ed", "te", "lsi")
adjacency_metrics <- c("pladj", "contag")

metric_family <- function(metric) {
  if (metric %in% composition_metrics) {
    return("Composition")
  }
  if (metric %in% edge_metrics) {
    return("Edge")
  }
  if (metric %in% adjacency_metrics) {
    return("Adjacency")
  }
  "Other"
}

summary_df$family <- vapply(summary_df$metric, metric_family, character(1))

family_summary <- unique(summary_df[, c(
  "size", "metric", "family", "method", "mean_sec",
  "median_sec", "sd_sec", "min_sec", "max_sec",
  "moving_mean_sec", "fft_mean_sec", "speedup_x"
)])

family_runtime <- stats::aggregate(
  cbind(mean_sec, median_sec) ~ size + family + method,
  data = family_summary,
  FUN = mean
)

family_speedup <- unique(family_summary[, c(
  "size", "family", "moving_mean_sec", "fft_mean_sec", "speedup_x"
)])
family_speedup <- stats::aggregate(
  cbind(moving_mean_sec, fft_mean_sec, speedup_x) ~ size + family,
  data = family_speedup,
  FUN = mean
)
family_speedup <- family_speedup[order(family_speedup$size, family_speedup$family), ]

selected_metrics <- c("shdi", "ed", "te", "pladj", "contag")
metric_500 <- unique(summary_df[
  summary_df$size == 500 & summary_df$metric %in% selected_metrics,
  c("metric", "method", "mean_sec", "moving_mean_sec", "fft_mean_sec", "speedup_x")
])
metric_500 <- metric_500[metric_500$method == "moving_window", ]
metric_500 <- metric_500[order(match(metric_500$metric, selected_metrics)), ]

fmt_num <- function(x, digits = 3) {
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

md_table <- function(df, digits = 3) {
  cols <- names(df)
  header <- paste0("| ", paste(cols, collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", length(cols)), collapse = " | "), " |")
  rows <- apply(df, 1, function(row) {
    paste0("| ", paste(row, collapse = " | "), " |")
  })
  c(header, sep, rows)
}

family_table <- data.frame(
  Size = family_speedup$size,
  Family = family_speedup$family,
  `Moving Mean (s)` = fmt_num(family_speedup$moving_mean_sec),
  `FFT Mean (s)` = fmt_num(family_speedup$fft_mean_sec),
  `Speedup (x)` = fmt_num(family_speedup$speedup_x, digits = 1),
  check.names = FALSE
)

metric_table <- data.frame(
  Metric = metric_500$metric,
  `Moving Mean (s)` = fmt_num(metric_500$moving_mean_sec),
  `FFT Mean (s)` = fmt_num(metric_500$fft_mean_sec),
  `Speedup (x)` = fmt_num(metric_500$speedup_x, digits = 1),
  check.names = FALSE
)

runtime_plot <- family_runtime
runtime_plot$family <- factor(runtime_plot$family,
                              levels = c("Composition", "Edge", "Adjacency"))
runtime_plot$method <- factor(runtime_plot$method,
                              levels = c("moving_window", "fft"),
                              labels = c("Moving window", "FFT"))

speedup_plot <- family_speedup
speedup_plot$family <- factor(speedup_plot$family,
                              levels = c("Composition", "Edge", "Adjacency"))

png_path <- file.path(output_dir, "landscape_metric_benchmark_runtime.png")

grDevices::png(filename = png_path, width = 1800, height = 900, res = 150)
op <- par(no.readonly = TRUE)
on.exit({
  par(op)
  grDevices::dev.off()
}, add = TRUE)

layout(matrix(c(1, 2), nrow = 1), widths = c(2.2, 1.1))

par(mar = c(5, 5, 3, 1))
plot(
  NA,
  xlim = range(runtime_plot$size),
  ylim = range(log10(runtime_plot$mean_sec)),
  xaxt = "n",
  yaxt = "n",
  xlab = "Raster size (cells per side)",
  ylab = "log10 mean runtime (seconds)",
  bty = "n"
)
axis(1, at = sort(unique(runtime_plot$size)))
y_ticks <- pretty(runtime_plot$mean_sec)
axis(2, at = log10(y_ticks[y_ticks > 0]), labels = fmt_num(y_ticks[y_ticks > 0], 2), las = 1)
title("Landscape metric benchmark: moving window vs FFT")
grid(nx = NA, ny = NULL, col = "grey90", lty = "solid")

family_cols <- c(
  "Composition" = "#2E86AB",
  "Edge" = "#F18F01",
  "Adjacency" = "#C73E1D"
)
method_lty <- c("Moving window" = 1, "FFT" = 2)
method_pch <- c("Moving window" = 16, "FFT" = 17)

for (fam in levels(runtime_plot$family)) {
  for (meth in levels(runtime_plot$method)) {
    d <- runtime_plot[runtime_plot$family == fam & runtime_plot$method == meth, ]
    lines(d$size, log10(d$mean_sec),
          col = family_cols[[fam]],
          lty = method_lty[[meth]],
          lwd = 2)
    points(d$size, log10(d$mean_sec),
           col = family_cols[[fam]],
           pch = method_pch[[meth]],
           cex = 1.1)
  }
}

legend("topleft",
       legend = levels(runtime_plot$family),
       col = unname(family_cols),
       lty = 1,
       lwd = 3,
       bty = "n",
       inset = 0.01,
       title = "Metric family")
legend("bottomright",
       legend = levels(runtime_plot$method),
       lty = unname(method_lty),
       pch = unname(method_pch),
       lwd = 2,
       bty = "n",
       inset = 0.01,
       title = "Method")

par(mar = c(5, 4, 3, 1))
bar_heights <- tapply(speedup_plot$speedup_x,
                      list(speedup_plot$family, speedup_plot$size),
                      identity)
barplot(
  bar_heights,
  beside = TRUE,
  col = unname(family_cols[rownames(bar_heights)]),
  ylab = "Mean speedup (x)",
  xlab = "Raster size",
  names.arg = colnames(bar_heights),
  las = 1,
  border = NA,
  main = "FFT speedup by family"
)

grDevices::dev.off()
on.exit(NULL, add = FALSE)

md_lines <- c(
  "# Landscape Metric Benchmark",
  "",
  "Benchmark configuration:",
  "- 10 independent categorical rasters per size",
  "- Raster sizes: 100x100, 250x250, 500x500",
  "- 4 land-cover classes",
  "- Circular moving window radius: 3 cells",
  "- Methods compared: conventional moving window (`terra::focal`) vs multiScaleR FFT whole-raster implementation",
  "",
  "## Family-level summary",
  ""
)
md_lines <- c(md_lines, md_table(family_table), "", "## Selected 500x500 metrics", "")
md_lines <- c(md_lines, md_table(metric_table), "", "## Figure", "")
md_lines <- c(md_lines,
              paste0("![Landscape metric benchmark runtime](./",
                     basename(png_path), ")"),
              "")

md_path <- file.path(output_dir, "landscape_metric_benchmark_summary.md")
writeLines(md_lines, con = md_path, useBytes = TRUE)

utils::write.csv(raw,
                 file = file.path(output_dir, "landscape_metric_benchmark_raw.csv"),
                 row.names = FALSE)
utils::write.csv(summary_df,
                 file = file.path(output_dir, "landscape_metric_benchmark_summary.csv"),
                 row.names = FALSE)

cat("WROTE_MD=", normalizePath(md_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("WROTE_PNG=", normalizePath(png_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
