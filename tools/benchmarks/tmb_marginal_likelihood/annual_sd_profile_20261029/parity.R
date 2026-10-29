.libPaths(c("tools/benchmarks/joint_scale_assessment/phase31/library",
  "C:/Users/peterman.73/.codex/scratch/multiScaleR-bayes-library", .libPaths()))
study <- "tools/benchmarks/tmb_marginal_likelihood/annual_sd_profile_20261029"
mode <- commandArgs(trailingOnly = TRUE)[1L]
path <- file.path(study, paste0(mode, ".json")); stopifnot(!file.exists(path))
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
for (helper in list.files("tests/testthat", "^helper.*[.]R$", full.names = TRUE)) source(helper)
fixture <- new.env(parent = globalenv())
eval(parse(file = "tests/testthat/test-factorized-sd-profile.R")[[1L]], fixture)
old <- new.env(parent = asNamespace("multiScaleR"))
sys.source(file.path(study, "predecessor/R/factorized_scale_profile.R"), old)
answer <- tryCatch({
  f <- fixture$response_sd_profile_fixture(); rng <- .Random.seed
  functions <- c(".msr_response_profile_problem", ".msr_response_profile_attempt",
    ".msr_response_profile_limits", ".msr_factorized_response_sd_profile")
  formals_same <- vapply(functions, function(nm)
    identical(formals(get(nm, old)), formals(get(nm, asNamespace("multiScaleR")))), TRUE)
  stopifnot(all(formals_same))
  rows <- list()
  for (sd in c(0, .2, .5)) {
    a <- old$.msr_response_profile_problem(f$spec, sd, f$parameters)
    b <- multiScaleR:::.msr_response_profile_problem(f$spec, sd, f$parameters)
    stopifnot(identical(a$par, b$par), identical(a$free, b$free), identical(a$fixed, b$fixed))
    for (offset in c(0, .01)) {
      x <- a$par + offset * sin(seq_along(a$par))
      same <- c(value = identical(a$fn(x), b$fn(x)), gradient = identical(a$gr(x), b$gr(x)),
        hessian = identical(a$he(x), b$he(x)), parameters = identical(a$parameters(x), b$parameters(x)))
      stopifnot(all(same))
      rows[[length(rows) + 1L]] <- list(sd = sd, offset = offset, exact = as.list(same))
    }
  }
  stopifnot(identical(.Random.seed, rng))
  list(passed = TRUE, fixed_vector_cases = rows, original_formals_unchanged = as.list(formals_same), RNG_unchanged = TRUE,
    scope = "Function-reference values/gradients/Hessians/extraction only; no closed producer or optimizer rerun.")
}, error = function(e) list(passed = FALSE, error = conditionMessage(e)))
writeLines(jsonlite::toJSON(answer, auto_unbox = TRUE), paste0(path, ".tmp"))
stopifnot(file.rename(paste0(path, ".tmp"), path), answer$passed)
