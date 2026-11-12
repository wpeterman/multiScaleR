# Stored results for vignettes/bayesian-computation.Rmd.
#
# Runs a short pilot and one long slice-sampling fit of the bundled five-year
# example and keeps only their summaries (diagnostics, precision reports, run
# time, object size). The fits themselves are not stored. Slice sampling needs no optional fitting
# packages; chain diagnostics need 'posterior' and parallel chains 'processx'.
#
# Run from the package root:  Rscript data-raw/cache_bayesian-computation.R
# Output: inst/extdata/bayesian-computation-cache.rds (xz compressed)

library(multiScaleR)

example <- readRDS(file.path("inst", "extdata", "bayesian-annual-example.rds"))
dat      <- example$data
inputs   <- example$inputs
response <- count ~ habitat + factor(year) + (1 | site)
scales   <- list(habitat = ~1 + (1 | year))

priors   <- list(scale_median = 250, scale_group_sd = 0.5)

# Keep the summaries of one slice fit that the vignette prints or cites
keep_summaries <- function(fit, seconds) {
  list(
    seconds     = seconds,
    sampling    = fit$sampling,                       # chains, iter, warmup, seed, cores
    fit_bytes   = as.numeric(object.size(fit)),
    chain_dim   = dim(fit$chains),                    # iterations, chains, parameters
    diagnostics = fit$diagnostics,                    # stored sampler checks
    precision   = diagnose_bayes_scales(fit)[c("available", "parameters",
                    "quantiles", "chains", "iterations", "derived_passed",
                    "passed", "criteria")],
    scales      = posterior_scales(fit),              # annual scales, meters
    group_sd    = summary(fit)$group_sd
  )
}

# 1. Short pilot run on one worker (object name in the vignette: pilot)
pilot_seconds <- system.time(
  pilot <- multiScale_bayes(
    response, inputs, dat, scales = scales, row_id = "site_year",
    family = "nbinom2", priors = priors,
    method = "slice", chains = 2, iter = 1000, seed = 5212
  )
)[["elapsed"]]

# 2. The same pilot on two workers, to record the cost of parallel startup
pilot2_seconds <- system.time(
  pilot2 <- multiScale_bayes(
    response, inputs, dat, scales = scales, row_id = "site_year",
    family = "nbinom2", priors = priors,
    method = "slice", chains = 2, iter = 1000, seed = 5212, cores = 2
  )
)[["elapsed"]]

# 3. Long run on four workers (object name in the vignette: slice_long)
long_seconds <- system.time(
  slice_long <- multiScale_bayes(
    response, inputs, dat, scales = scales, row_id = "site_year",
    family = "nbinom2", priors = priors,
    method = "slice", chains = 4, iter = 8000, seed = 5213, cores = 4
  )
)[["elapsed"]]

cache <- list(
  pilot = keep_summaries(pilot, pilot_seconds),
  pilot_two_cores = list(seconds = pilot2_seconds,
                         same_draws = identical(pilot$chains, pilot2$chains)),
  long  = keep_summaries(slice_long, long_seconds),
  # Precision of P(sigma <= threshold) for each year, from the long run
  probabilities = diagnose_bayes_probabilities(
    slice_long, thresholds = list(habitat = c(100, 250)), mcse_target = 0.005
  )[c("available", "probabilities", "mcse_target")],
  session = list(multiScaleR = as.character(packageVersion("multiScaleR")),
                 R = R.version.string, date = Sys.Date())
)

saveRDS(cache, file.path("inst", "extdata", "bayesian-computation-cache.rds"),
        compress = "xz")
