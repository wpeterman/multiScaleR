family_path <- file.path("tools", "benchmarks", "landscape_metric_benchmark_family_summary.csv")
metric_path <- file.path("tools", "benchmarks", "landscape_metric_benchmark_selected_500_metrics.csv")
md_path <- file.path("tools", "benchmarks", "landscape_metric_benchmark_summary.md")
png_path <- file.path("tools", "benchmarks", "landscape_metric_benchmark_runtime.png")

family_df <- utils::read.csv(family_path, stringsAsFactors = FALSE)
metric_df <- utils::read.csv(metric_path, stringsAsFactors = FALSE)

family_df$family <- factor(family_df$family,
                           levels = c("Composition", "Edge", "Adjacency"))
family_df$method <- factor(family_df$method,
                           levels = c("Moving window", "FFT"))

speedup_df <- reshape(
  family_df[, c("size", "family", "method", "mean_sec")],
  idvar = c("size", "family"),
  timevar = "method",
  direction = "wide"
)
names(speedup_df) <- sub("^mean_sec\\.", "", names(speedup_df))
speedup_df$speedup_x <- speedup_df$`Moving window` / speedup_df$FFT
speedup_df <- speedup_df[order(speedup_df$size, speedup_df$family), ]

fmt_num <- function(x, digits = 3) {
  formatC(x, format = "f", digits = digits)
}

md_table <- function(df) {
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  c(header, sep, rows)
}

family_table <- data.frame(
  Size = speedup_df$size,
  Family = as.character(speedup_df$family),
  `Moving Mean (s)` = fmt_num(speedup_df$`Moving window`, 3),
  `FFT Mean (s)` = fmt_num(speedup_df$FFT, 3),
  `Speedup (x)` = fmt_num(speedup_df$speedup_x, 1),
  check.names = FALSE
)

metric_table <- data.frame(
  Metric = metric_df$metric,
  `Moving Mean (s)` = fmt_num(metric_df$moving_mean_sec, 3),
  `FFT Mean (s)` = fmt_num(metric_df$fft_mean_sec, 3),
  `Speedup (x)` = fmt_num(metric_df$speedup_x, 1),
  check.names = FALSE
)

grDevices::png(filename = png_path, width = 1800, height = 900, res = 150)
op <- par(no.readonly = TRUE)
layout(matrix(c(1, 2), nrow = 1), widths = c(2.1, 1.1))

family_cols <- c(
  "Composition" = "#2E86AB",
  "Edge" = "#F18F01",
  "Adjacency" = "#C73E1D"
)
method_lty <- c("Moving window" = 1, "FFT" = 2)
method_pch <- c("Moving window" = 16, "FFT" = 17)

par(mar = c(5, 5, 3, 1))
plot(
  NA,
  xlim = range(family_df$size),
  ylim = range(log10(family_df$mean_sec)),
  xaxt = "n",
  yaxt = "n",
  xlab = "Raster size (cells per side)",
  ylab = "log10 mean runtime (seconds)",
  bty = "n"
)
axis(1, at = sort(unique(family_df$size)))
y_ticks <- c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100)
axis(2, at = log10(y_ticks), labels = fmt_num(y_ticks, 2), las = 1)
title("Landscape metric benchmark: moving window vs FFT")
grid(nx = NA, ny = NULL, col = "grey90", lty = "solid")

for (fam in levels(family_df$family)) {
  for (meth in levels(family_df$method)) {
    d <- family_df[family_df$family == fam & family_df$method == meth, ]
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
       legend = levels(family_df$family),
       col = unname(family_cols),
       lty = 1,
       lwd = 3,
       bty = "n",
       inset = 0.01,
       title = "Metric family")
legend("bottomright",
       legend = levels(family_df$method),
       lty = unname(method_lty),
       pch = unname(method_pch),
       lwd = 2,
       bty = "n",
       inset = 0.01,
       title = "Method")

par(mar = c(5, 4, 3, 1))
bar_heights <- tapply(speedup_df$speedup_x,
                      list(speedup_df$family, speedup_df$size),
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

par(op)
grDevices::dev.off()

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

writeLines(md_lines, con = md_path, useBytes = TRUE)

cat("WROTE_MD=", normalizePath(md_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("WROTE_PNG=", normalizePath(png_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
