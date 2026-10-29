.libPaths(c("tools/benchmarks/joint_scale_assessment/phase31/library",
  "C:/Users/peterman.73/.codex/scratch/multiScaleR-bayes-library", .libPaths()))
study <- "tools/benchmarks/tmb_marginal_likelihood/probability_precision_20261029"
stopifnot(!file.exists(file.path(study, "parity01.json")))
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
source("tests/testthat/helper-bayesian-scales.R")
old_environment <- new.env(parent = asNamespace("multiScaleR"))
sys.source(file.path(study, "predecessor/R/bayes_scale_probabilities.R"), envir = old_environment)
old_function <- old_environment$diagnose_bayes_probabilities
stopifnot(identical(formals(old_function), formals(diagnose_bayes_probabilities)))
check_case <- function(random, kernel, balanced = FALSE) {
  f <- bayesian_scale_fixture(n_site = 3L, n_year = 2L, kernel = kernel)
  fit <- if (random) prepare_bayesian_fixture(f) else
    multiScale_bayes(count ~ habitat, f$inputs, f$data, row_id = "key", method = "prepare")
  fit$sampler <- "parity fixture only"
  fit$diagnostics <- list(passed = FALSE, marker = "preserved raw failure")
  columns <- unlist(lapply(names(fit$parameters), function(nm) {
    n <- length(fit$parameters[[nm]])
    if (!n) character() else if (n == 1L) nm else paste0(nm, "[", seq_len(n), "]")
  }))
  n <- 600L * 4L
  set.seed(89031)
  values <- matrix(rep(unlist(fit$parameters), each = n), nrow = n)
  values <- values + matrix(rnorm(length(values), sd = .3), nrow = n)
  if (balanced) {
    fit$specification$unit_conv <- 1
    values[, match("scale_beta", columns)] <- sample(rep(c(-1, 1), each = n / 2L))
  }
  raw <- array(values, c(600L, 4L, ncol(values)), dimnames = list(NULL, NULL, columns))
  testthat::local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) raw, .package = "multiScaleR")
  args <- list(object = fit, thresholds = list(habitat = if (balanced) 0 else c(100, 250, 400)),
    scale = if (balanced) "log" else "distance", mcse_target = .1)
  before <- serialize(fit, NULL); rng <- .Random.seed
  old <- do.call(old_function, args)
  stopifnot(identical(.Random.seed, rng), identical(serialize(fit, NULL), before))
  current <- do.call(diagnose_bayes_probabilities, args)
  stopifnot(identical(.Random.seed, rng), identical(serialize(fit, NULL), before))
  fixed <- setdiff(names(old$probabilities), c("precision_available", "precision_target_met", "reason"))
  stopifnot(identical(old$probabilities[fixed], current$probabilities[fixed]),
    identical(old$chain_probabilities, current$chain_probabilities))
  root <- setdiff(names(old), "probabilities")
  stopifnot(identical(old[root], current[root]))
  p <- old$probabilities
  changed <- !is.finite(p$rhat) & !p$constant & is.finite(p$mcse) & p$mcse > 0 &
    is.finite(p$ess_mean) & p$ess_mean > 0
  stopifnot(identical(old$probabilities[!changed, , drop = FALSE],
    current$probabilities[!changed, names(old$probabilities), drop = FALSE]))
  if (any(changed)) stopifnot(all(current$probabilities$precision_available[changed]),
    !any(current$probabilities$rhat_available[changed]),
    all(current$probabilities$precision_target_met[changed] == (p$mcse[changed] <= .1)),
    all(current$probabilities$reason[changed] == ""))
  list(random = random, kernel = kernel, balanced = balanced,
    numeric_and_chain_outputs_identical = TRUE, root_output_and_raw_diagnostics_identical = TRUE,
    expected_availability_corrections = sum(changed), fit_RNG_unchanged = TRUE)
}
results <- list()
for (kernel in c("gaussian", "exponential")) for (random in c(FALSE, TRUE))
  results[[length(results) + 1L]] <- check_case(random, kernel)
results[[length(results) + 1L]] <- check_case(FALSE, "gaussian", TRUE)
record <- list(passed = TRUE, public_formals_identical = TRUE, results = results,
  source = as.list(tools::md5sum(c("R/bayes_scale_probabilities.R",
    file.path(study, "predecessor/R/bayes_scale_probabilities.R")))))
temporary <- file.path(study, "parity01.json.tmp")
writeLines(jsonlite::toJSON(record, auto_unbox = TRUE), temporary)
stopifnot(file.rename(temporary, file.path(study, "parity01.json")))
