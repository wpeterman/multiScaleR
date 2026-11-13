# Strip YAML and code chunks from an Rmd so the style audit sees prose only.
# Usage: Rscript prose_only.R <in.Rmd> <out.md>
args <- commandArgs(trailingOnly = TRUE)
x <- readLines(args[1], warn = FALSE, encoding = "UTF-8")
yaml <- which(x == "---")
if (length(yaml) >= 2) x <- x[-(yaml[1]:yaml[2])]
fence <- cumsum(grepl("^```", x))
x <- x[fence %% 2 == 0 & !grepl("^```", x)]
x <- x[!grepl("^\\s*\\|", x)]          # drop table rows
x <- gsub("`r [^`]*`", "N", x)         # inline R becomes a number placeholder
writeLines(x, args[2], useBytes = TRUE)
cat("em-dashes:", sum(grepl(intToUtf8(0x2014), readLines(args[1], warn = FALSE, encoding = "UTF-8"), fixed = TRUE)), "\n")
cat("prose words:", sum(lengths(strsplit(trimws(x[nzchar(trimws(x))]), "\\s+"))), "\n")
