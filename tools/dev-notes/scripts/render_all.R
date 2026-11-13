# Render every vignette against the private build, one fresh R process each, and record time/status.
scr <- "C:/Users/peterman.73/AppData/Local/Temp/claude/9e3b7d3898bd820b/459cb932-8333-43df-b0da-7d9a2422aa59/scratchpad"
rs <- file.path(R.home("bin"), "Rscript.exe")
out <- file.path(scr, "render_final")
for (f in sort(list.files("vignettes", pattern = "[.]Rmd$"))) {
  t0 <- Sys.time()
  log <- suppressWarnings(system2(rs, c("--vanilla", shQuote(file.path(scr, "render_one.R")), f, shQuote(out)),
                 stdout = TRUE, stderr = TRUE))
  st <- attr(log, "status"); if (is.null(st)) st <- 0
  cat(sprintf("%-36s exit %d  %5.1fs  warnings:%d\n", f, st,
      as.numeric(difftime(Sys.time(), t0, units = "secs")), sum(grepl("Warning", log))))
  if (st != 0) cat(tail(log, 8), sep = "\n")
}
