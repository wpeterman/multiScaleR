# Stored results for vignettes/bayesian-scales.Rmd.
#
# The vignette shows one slice-sampling fit to the running example
# (landscape_counts, surv_pts, landscape.tif). The fit takes roughly 20 to 40
# seconds, so it is stored here and the vignette shows the call unevaluated.
# The stored fit keeps its ordered chains, so diagnose_bayes_scales() and
# chain plots work on it.
#
# Run from the package root with the current build of multiScaleR installed
# and the 'posterior' package available:
#   Rscript data-raw/cache_bayesian-scales.R
# Output: inst/extdata/bayesian-scales-cache.rds (xz compressed).

library(multiScaleR)
stopifnot(requireNamespace("posterior", quietly = TRUE))

data("landscape_counts")
data("surv_pts")
land_rast <- terra::rast(system.file("extdata", "landscape.tif",
                                     package = "multiScaleR"))

# Same preparation as the vignette: binned summaries within 1700 m
inputs <- kernel_prep(pts = terra::vect(surv_pts), raster_stack = land_rast,
                      max_D = 1700, kernel = "gaussian", bin = TRUE,
                      verbose = FALSE)

# Same call as the vignette chunk `fit-shared`
t0 <- Sys.time()
fit <- multiScale_bayes(counts ~ site + land1 + land2, inputs,
                        landscape_counts, family = "poisson",
                        priors = list(scale_median = c(land1 = 300,
                                                       land2 = 300)),
                        method = "slice", seed = 5211)
seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Chain diagnostics for the fitted log distances (needs 'posterior')
precision <- diagnose_bayes_scales(fit)
stopifnot(isTRUE(precision$available))

# The same model with the reference scales moved to 250 m and 500 m. The
# vignette quotes its coefficients to show that a Bayesian coefficient refers
# to one standard deviation of the covariate at scale_median. Only the two
# summary tables are stored.
fit_ref <- multiScale_bayes(counts ~ site + land1 + land2, inputs,
                            landscape_counts, family = "poisson",
                            priors = list(scale_median = c(land1 = 250,
                                                           land2 = 500)),
                            method = "slice", seed = 5211)
reference_check <- list(coefficients = summary(fit_ref)$coefficients,
                        scales = posterior_scales(fit_ref),
                        passed = summary(fit_ref)$diagnostics$passed)

cache <- list(fit = fit, precision = precision, seconds = seconds,
              reference_check = reference_check,
              package_version = as.character(packageVersion("multiScaleR")),
              created = Sys.Date())
saveRDS(cache, file.path("inst", "extdata", "bayesian-scales-cache.rds"),
        compress = "xz")
cat("seconds:", round(seconds, 1), "\n")
cat("bytes:", file.size(file.path("inst", "extdata",
                                  "bayesian-scales-cache.rds")), "\n")
