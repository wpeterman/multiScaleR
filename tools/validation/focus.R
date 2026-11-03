input <- jsonlite::fromJSON(commandArgs(TRUE)[[1L]])
.libPaths(c(input$dependency_overlay, .libPaths()))
invisible(lapply(list.files(file.path(input$package, "R"), pattern = "[.]R$", full.names = TRUE), parse))
invisible(lapply(list.files(file.path(input$package, "tests/testthat"), pattern = "[.]R$", full.names = TRUE), parse))
# Keep requested metadata and all test tasks; bound actual validation sockets.
trace("makePSOCKcluster", where = asNamespace("parallel"), print = FALSE,
  tracer = quote({ if (is.numeric(names) && length(names) == 1L) names <- min(names, 2L)
                  else if (length(names) > 2L) names <- names[seq_len(2L)] }))
result <- testthat::test_local(input$package, filter = input$filter,
  reporter = "summary", stop_on_failure = FALSE)
d <- as.data.frame(result)
counts <- list(passed = sum(d$passed), failed = sum(d$failed), errors = sum(d$error),
  warnings = sum(d$warning), skips = sum(d$skipped))
value <- list(counts = counts, by_file = lapply(split(seq_len(nrow(d)), basename(d$file)), function(i)
  list(passed = sum(d$passed[i]), failed = sum(d$failed[i]), errors = sum(d$error[i]),
    warnings = sum(d$warning[i]), skips = sum(d$skipped[i]))))
stopifnot(!file.exists(input$result), !file.exists(paste0(input$result, ".tmp")))
jsonlite::write_json(value, paste0(input$result, ".tmp"), auto_unbox = TRUE, pretty = TRUE)
stopifnot(file.rename(paste0(input$result, ".tmp"), input$result))
if (counts$failed + counts$errors) stop("Retain the focused failures; do not replay this phase.")
