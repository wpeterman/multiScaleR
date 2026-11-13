# Convert a rendered vignette to plain text and extract its figures for inspection.
# Usage: Rscript html_inspect.R <html-file> <out-prefix>
args <- commandArgs(trailingOnly = TRUE)
html <- args[1]; prefix <- args[2]
txt <- paste0(prefix, ".txt")
rmarkdown::pandoc_convert(html, to = "plain", from = "html", output = txt,
                          options = c("--wrap=none"))
x <- paste(readLines(html, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
m <- regmatches(x, gregexpr("data:image/png;base64,[A-Za-z0-9+/=]+", x))[[1]]
cat("figures:", length(m), "\n")
for (i in seq_along(m)) {
  raw <- jsonlite::base64_dec(sub("^data:image/png;base64,", "", m[i]))
  writeBin(raw, sprintf("%s_fig%02d.png", prefix, i))
}
cat("text lines:", length(readLines(txt, warn = FALSE)), "\n")
