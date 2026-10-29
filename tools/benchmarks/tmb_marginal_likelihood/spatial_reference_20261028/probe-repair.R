study <- normalizePath("tools/benchmarks/tmb_marginal_likelihood/spatial_reference_20261028")
preceding <- "tools/benchmarks/tmb_marginal_likelihood/functional_precision_20261028"
qualified <- jsonlite::read_json(file.path(preceding, "check01.json"), simplifyVector = TRUE)
snapshot <- jsonlite::read_json(file.path(preceding, "snapshot01-input.json"), simplifyVector = TRUE)
libraries <- jsonlite::read_json(file.path(preceding, "compatibility01.json"), simplifyVector = TRUE)$libraries
.libPaths(c(file.path(qualified$workspace, "multiScaleR.Rcheck"), libraries, .libPaths()))
library(multiScaleR)
source(file.path(snapshot$package, "tests/testthat/helper-bayesian-scales.R"))
source(file.path(study, "adapter.R"))
stopifnot(length(parse(file.path(study, "run-fit-v2.R"))) > 0L)
case <- sr_case("null")
original <- multiScaleR:::.msr_spatial_scale_adfun(case$object)
original_initial <- lapply(seq_len(4L), function(i) unname(original$par) + .05 * sin(seq_along(original$par) + i))
invisible(original$fn(original_initial[[1L]]))
lost_names <- is.null(names(original$env$last.par.best))
objective <- multiScaleR:::.msr_spatial_scale_adfun(case$object)
named_initial <- lapply(seq_len(4L), function(i) objective$par + .05 * sin(seq_along(objective$par) + i))
initial <- lapply(named_initial, unname)
values <- vapply(named_initial, objective$fn, 0.0)
gradients <- lapply(named_initial, objective$gr)
retained_names <- identical(names(objective$env$last.par.best), names(objective$par))
# Test the exact tmbstan init-sanitizer indexing operation without sampling.
sanitized <- list()
sanitized[[names(objective$env$last.par.best)[1L]]] <- numeric()
sanitized$y <- initial[[1L]]
record <- list(passed = lost_names && retained_names && identical(initial, original_initial) &&
    all(is.finite(values)) && all(vapply(gradients, function(x) all(is.finite(x)), TRUE)),
  original_unnamed_preflight_loses_RTMB_names = lost_names, named_preflight_retains_names = retained_names,
  all_four_numeric_starts_identical = identical(initial, original_initial),
  values = values, gradient_max = vapply(gradients, function(x) max(abs(x)), 0.0),
  tmbstan_initial_mapping_valid = identical(sanitized$y, original_initial[[1L]]),
  driver_parses = TRUE, note = "No sampler called; no target, prior, numeric start or gate changed.",
  repaired_driver_md5 = unname(tools::md5sum(file.path(study, "run-fit-v2.R"))))
sr_atomic(record, file.path(study, "repair-unit01.json"))
print(record)
stopifnot(record$passed)
