input <- jsonlite::fromJSON(commandArgs(TRUE)[[1L]])
.libPaths(c(dirname(input$installed_library), input$dependency_overlay, .libPaths()))
stopifnot(identical(normalizePath(find.package("MASS"), winslash = "/"),
  normalizePath(file.path(input$dependency_overlay, "MASS"), winslash = "/")))
stopifnot(identical(normalizePath(find.package("multiScaleR"), winslash = "/"),
  normalizePath(input$installed_library, winslash = "/")))
library(multiScaleR)
cache <- get(".msr_hw_cache", asNamespace("multiScaleR"))
cache$physical_cores <- 4L; cache$logical_cores <- 4L
# Preserve every test/task and requested worker metadata while bounding the
# actual socket pool to two workers beside this single validation driver.
trace("makePSOCKcluster", where = asNamespace("parallel"), print = FALSE,
  tracer = quote({ if (is.numeric(names) && length(names) == 1L) names <- min(names, 2L)
                  else if (length(names) > 2L) names <- names[seq_len(2L)] }))
setwd(input$exact_archive_tests)
result <- source("testthat.R", echo = FALSE)
d <- as.data.frame(result$value)
value <- lapply(split(seq_len(nrow(d)), basename(d$file)), function(i)
  list(passed = sum(d$passed[i]), failed = sum(d$failed[i]), errors = sum(d$error[i]),
    warnings = sum(d$warning[i]), skips = sum(d$skipped[i])))
stopifnot(!file.exists(input$result), !file.exists(paste0(input$result, ".tmp")))
jsonlite::write_json(value, paste0(input$result, ".tmp"), auto_unbox = TRUE, pretty = TRUE)
stopifnot(file.rename(paste0(input$result, ".tmp"), input$result))
