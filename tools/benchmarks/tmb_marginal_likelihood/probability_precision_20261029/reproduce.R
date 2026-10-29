study <- "tools/benchmarks/tmb_marginal_likelihood/probability_precision_20261029"
prior <- "tools/benchmarks/tmb_marginal_likelihood/functional_precision_20261028"
qualified <- jsonlite::read_json(file.path(prior, "check01.json"), simplifyVector = TRUE)
snapshot <- jsonlite::read_json(file.path(prior, "snapshot01-input.json"), simplifyVector = TRUE)
libraries <- jsonlite::read_json(file.path(prior, "compatibility01.json"), simplifyVector = TRUE)$libraries
.libPaths(c(file.path(qualified$workspace, "multiScaleR.Rcheck"), libraries, .libPaths()))
library(multiScaleR)
source(file.path(snapshot$package, "tests/testthat/helper-bayesian-scales.R"))
stopifnot(!file.exists(file.path(study, "reproduction01.json")))
namespace <- asNamespace("posterior")
functions <- c("rhat", "rhat.default", ".rhat", ".rhat_basic", ".split_chains", ".z_scale")
functions <- functions[vapply(functions, exists, TRUE, envir = namespace, inherits = FALSE)]
source_text <- unlist(lapply(functions, function(nm)
  c(paste0("# posterior:::", nm), deparse(get(nm, envir = namespace, inherits = FALSE)), "")))
writeLines(source_text, file.path(study, "posterior-source.txt"))
reproduce <- function() {
  f <- bayesian_scale_fixture(n_site = 3L, n_year = 2L)
  fit <- multiScale_bayes(count ~ habitat, f$inputs, f$data, row_id = "key", method = "prepare")
  fit$specification$unit_conv <- 1
  fit$sampler <- "deterministic chain fixture only"
  fit$diagnostics <- list(passed = FALSE, marker = "must remain failed")
  n <- 400L * 4L
  columns <- unlist(lapply(names(fit$parameters), function(nm) {
    q <- length(fit$parameters[[nm]])
    if (!q) character() else if (q == 1L) nm else paste0(nm, "[", seq_len(q), "]")
  }))
  values <- matrix(rep(unlist(fit$parameters), each = n), nrow = n)
  set.seed(917613)
  values[, match("scale_beta", columns)] <- sample(rep(c(-1, 1), each = n / 2L))
  raw <- array(values, c(400L, 4L, ncol(values)), dimnames = list(NULL, NULL, columns))
  indicator <- 1 * (raw[, , "scale_beta"] <= 0)
  expected <- c(estimate = mean(indicator), mcse = posterior::mcse_mean(indicator),
    ess_mean = posterior::ess_mean(indicator), rhat = posterior::rhat(indicator))
  testthat::local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) raw, .package = "multiScaleR")
  before <- serialize(fit, NULL); rng <- .Random.seed
  report <- diagnose_bayes_probabilities(fit, list(habitat = 0), scale = "log", mcse_target = .1)
  list(expected = as.list(expected), actual = report$probabilities,
    expected_mean_precision_available = all(is.finite(expected[c("mcse", "ess_mean")])) &&
      all(expected[c("mcse", "ess_mean")] > 0),
    fit_unchanged = identical(before, serialize(fit, NULL)), RNG_unchanged = identical(rng, .Random.seed),
    raw_failure_retained = identical(report$raw_diagnostics, fit$diagnostics),
    overall_passed_flag_absent = !"passed" %in% names(report))
}
result <- reproduce()
result$regression_reproduced <- isTRUE(result$expected_mean_precision_available) &&
  !isTRUE(result$actual$precision_available)
result$posterior_version <- as.character(packageVersion("posterior"))
result$R <- R.version.string
result$installed_package <- find.package("multiScaleR")
result$source_md5 <- as.list(tools::md5sum(c("R/bayes_scale_probabilities.R", file.path(study, "posterior-source.txt"))))
result$scope <- "Deterministic original-chain fixture; no sampler, simulation panel or scientific dataset. The expected old reporting failure is retained."
temporary <- file.path(study, "reproduction01.json.tmp")
writeLines(jsonlite::toJSON(result, auto_unbox = TRUE, null = "null", digits = 16), temporary)
stopifnot(file.rename(temporary, file.path(study, "reproduction01.json")))
print(result[c("expected", "actual", "regression_reproduced", "posterior_version")])
stopifnot(identical(isTRUE(result$actual$precision_available), isTRUE(result$expected_mean_precision_available)))
