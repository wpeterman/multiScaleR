study <- normalizePath("tools/benchmarks/tmb_marginal_likelihood/spatial_reference_20261028")
preceding <- "tools/benchmarks/tmb_marginal_likelihood/functional_precision_20261028"
qualified <- jsonlite::read_json(file.path(preceding, "check01.json"), simplifyVector = TRUE)
snapshot <- jsonlite::read_json(file.path(preceding, "snapshot01-input.json"), simplifyVector = TRUE)
libraries <- jsonlite::read_json(file.path(preceding, "compatibility01.json"), simplifyVector = TRUE)$libraries
.libPaths(c(file.path(qualified$workspace, "multiScaleR.Rcheck"), libraries, .libPaths()))
library(multiScaleR)
source(file.path(snapshot$package, "tests/testthat/helper-bayesian-scales.R"))
source(file.path(study, "adapter.R"))
stage <- "case"; error_calls <- character()
record <- tryCatch(withCallingHandlers({
  case <- sr_case("null")
  stage <- "objective"
  objective <- multiScaleR:::.msr_spatial_scale_adfun(case$object)
  stage <- "initial"
  initial <- lapply(seq_len(4L), function(i) unname(objective$par) + .05 * sin(seq_along(objective$par) + i))
  stage <- "value"
  values <- vapply(initial, objective$fn, 0.0)
  stage <- "gradient"
  gradients <- lapply(initial, objective$gr)
  list(passed = all(is.finite(values)) && all(vapply(gradients, function(x) all(is.finite(x)), TRUE)),
    stage = stage, values = values, gradient_max = vapply(gradients, function(x) max(abs(x)), 0.0),
    parameter_names = names(objective$par), objective_class = class(objective),
    note = "No sampler called.")
}, error = function(e) { error_calls <<- vapply(sys.calls(), function(x) paste(deparse(x), collapse = " "), "") }),
  error = function(e) list(passed = FALSE, stage = stage, error = conditionMessage(e), calls = error_calls,
    RTMB = as.character(packageVersion("RTMB")), TMB = as.character(packageVersion("TMB")),
    note = "No sampler called."))
sr_atomic(record, file.path(study, "probe01.json"))
print(record)
