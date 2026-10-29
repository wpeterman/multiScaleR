.libPaths(c("tools/benchmarks/joint_scale_assessment/phase31/library",
  "C:/Users/peterman.73/.codex/scratch/multiScaleR-bayes-library", .libPaths()))
study <- "tools/benchmarks/tmb_marginal_likelihood/prediction_contracts_20261028"
stopifnot(!file.exists(file.path(study, "parity01.json")))
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
for (helper in list.files("tests/testthat", "^helper.*[.]R$", full.names = TRUE)) source(helper)
legacy <- new.env(parent = asNamespace("multiScaleR"))
for (name in c("bayes_scale_methods.R", "bayes_scale_predict_scales.R", "bayes_scale_simulate.R",
               "prepare_prediction.R", "prepare_native_prediction.R", "factorized_scale_posterior.R"))
  sys.source(file.path(study, "predecessor/R", name), legacy)
functions <- c("predict.multiScaleR_bayes", "predict_scales", "simulate.multiScaleR_bayes",
               "prepare_prediction", "predict.multiScaleR_prediction", "print.multiScaleR_prediction")
stopifnot(all(vapply(functions, function(nm)
  identical(formals(get(nm, legacy)), formals(get(nm, asNamespace("multiScaleR")))), TRUE)))
strip <- function(x) { attr(x, "prediction_contract") <- NULL; x }
records <- list()
for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3, kernel = kernel)
  if (family == "binomial") f$data$count <- as.numeric(f$data$count > 2)
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f, family), 9)
  fit$draws$beta[,] <- matrix(seq_len(length(fit$draws$beta)) / 100, 9)
  nd <- f$data[c(4, 2, 4), ]; nd$year <- c("future", "y1", "future")
  nd$site <- c("new", "s2", "new")
  initial <- fit
  for (data in list(NULL, nd)) {
    args <- list(object = fit, newdata = data, kernel_inputs = if (is.null(data)) NULL else f$inputs,
      new_scale_levels = "sample", new_response_levels = "sample", seed = 407)
    for (summary in c(TRUE, FALSE)) for (type in c("response", "link")) {
      current <- do.call(predict.multiScaleR_bayes, c(args, list(summary = summary, type = type)))
      old <- do.call(legacy$predict.multiScaleR_bayes, c(args, list(summary = summary, type = type)))
      stopifnot(identical(strip(current), old))
    }
    set.seed(102); state <- .Random.seed
    current <- do.call(simulate.multiScaleR_bayes, c(args, list(nsim = 3)))
    stopifnot(identical(.Random.seed, state))
    old <- do.call(legacy$simulate.multiScaleR_bayes, c(args, list(nsim = 3)))
    stopifnot(identical(strip(current), old), identical(.Random.seed, state))
  }
  for (summary in c(TRUE, FALSE)) for (scale in c("distance", "log")) {
    args <- list(object = fit, newdata = nd, new_levels = "sample", seed = 408,
                 summary = summary, scale = scale)
    stopifnot(identical(strip(do.call(predict_scales, args)), do.call(legacy$predict_scales, args)))
  }
  stopifnot(identical(initial, fit))
  records[[length(records) + 1L]] <- list(kernel = kernel, family = family,
    current_and_predecessor_identical_except_contract = TRUE, RNG_preserved = TRUE)
}
# The candidate's full numerical parameter-map description is structural only.
# Fixed-vector likelihood reconstruction must retain the same conditional values.
f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp, f$inputs, f$data,
  list(habitat = ~1), "gaussian", "key", 200)
p <- s$parameters; p$scale_z <- numeric()
current <- multiScaleR:::.msr_factorized_reconstruct(s, p)
old <- legacy$.msr_factorized_reconstruct(s, p)
current$prediction_contract <- NULL
stopifnot(identical(current, old))
out <- list(passed = TRUE, public_formals_unchanged = functions, posterior_cases = records,
            conditional_reconstruction_identical = TRUE)
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE), file.path(study, "parity01.json.tmp"))
stopifnot(file.rename(file.path(study, "parity01.json.tmp"), file.path(study, "parity01.json")))
