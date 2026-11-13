# Render one vignette against the private scratch build.
# Usage: Rscript render_one.R <vignette-file-name> <output-dir>
args <- commandArgs(trailingOnly = TRUE)
scr <- "C:/Users/peterman.73/AppData/Local/Temp/claude/9e3b7d3898bd820b/459cb932-8333-43df-b0da-7d9a2422aa59/scratchpad"
.libPaths(c(file.path(scr, "rlib"), .libPaths()))
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
res <- rmarkdown::render(file.path("vignettes", args[1]), output_dir = out_dir,
                         intermediates_dir = tempdir(), quiet = TRUE,
                         envir = new.env())
cat("rendered:", res, "\n")
cat("seconds:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "\n")
