# Focused feedback from a genuine private installation. This avoids source-only
# package-environment serialization artifacts without suppressing warnings.
setup_started <- proc.time()[["elapsed"]]
input <- jsonlite::fromJSON(commandArgs(TRUE)[1L])
audit_root <- input$audit_root
if (is.null(audit_root)) audit_root <- file.path(dirname(input$result), "a")
dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(MSR_GDR_TEST_AUDIT_ROOT = normalizePath(audit_root, winslash = "/", mustWork = TRUE))
if (!is.character(input$filter) || length(input$filter) != 1L || !nzchar(input$filter))
  stop("Installed focused feedback requires an explicit nonempty filter.")
.libPaths(c(dirname(input$installed_library), input$dependency_overlay, .libPaths()))
stopifnot(identical(normalizePath(find.package("multiScaleR"), winslash = "/"),
  normalizePath(input$installed_library, winslash = "/")))
library(multiScaleR)
cache <- get(".msr_hw_cache", asNamespace("multiScaleR"))
cache$physical_cores <- 4L; cache$logical_cores <- 4L
trace("makePSOCKcluster", where = asNamespace("parallel"), print = FALSE,
  tracer = quote({ if (is.numeric(names) && length(names) == 1L) names <- min(names, 2L)
                  else if (length(names) > 2L) names <- names[seq_len(2L)] }))
setwd(input$tests_directory)
test_started <- proc.time()[["elapsed"]]
result <- testthat::test_check("multiScaleR", filter = input$filter,
  reporter = "summary", stop_on_failure = FALSE)
elapsed <- proc.time()[["elapsed"]] - test_started
d <- as.data.frame(result)
counts <- list(passed = sum(d$passed), failed = sum(d$failed), errors = sum(d$error),
  warnings = sum(d$warning), skips = sum(d$skipped))
value <- list(counts = counts, by_file = lapply(split(seq_len(nrow(d)), basename(d$file)), function(i)
  list(passed = sum(d$passed[i]), failed = sum(d$failed[i]), errors = sum(d$error[i]),
    warnings = sum(d$warning[i]), skips = sum(d$skipped[i]))),
  mode = "installed_focused", release_eligible = FALSE)
stopifnot(!file.exists(input$result), !file.exists(paste0(input$result, ".tmp")))
jsonlite::write_json(value, paste0(input$result, ".tmp"), pretty = TRUE, auto_unbox = TRUE)
stopifnot(file.rename(paste0(input$result, ".tmp"), input$result))
source(file.path(Sys.getenv("MSR_VALIDATION_RUNTIME"), "test_metrics.R"), local = TRUE)
msr_write_test_timings(d, input$result, elapsed, test_started - setup_started)
if (!setequal(names(value$by_file), input$expected_files)) stop("Focused file selection differed from its plan.")
if (counts$failed + counts$errors) stop("Retain the focused failures; declare a new continuation.")
