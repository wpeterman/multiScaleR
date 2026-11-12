# Stored results for vignettes/spatial-response-fields.Rmd.
#
# The vignette fits a Poisson count model with a glmmTMB spatial field through
# multiScale_optim(). Each candidate scale needs a glmmTMB refit, so the search
# takes most of a minute. The fit is stored here and the vignette shows the
# call unevaluated. Everything else in the vignette runs live.
#
# The stored fit drops its copy of the kernel inputs (about 44 MB in memory).
# The vignette attaches the inputs it prepares itself, which are identical.
#
# Run from the package root with the current build of multiScaleR installed
# and the 'glmmTMB' package available:
#   Rscript data-raw/cache_spatial-response-fields.R
# Output: inst/extdata/spatial-response-fields-cache.rds (xz compressed).

library(multiScaleR)
library(terra)
stopifnot(requireNamespace("glmmTMB", quietly = TRUE))

data("landscape_counts")
data("surv_pts")
pts <- vect(surv_pts)
xy  <- crds(pts)
land_rast <- rast(system.file("extdata", "landscape.tif",
                              package = "multiScaleR"))

# ---- Simulated counts: same code as the vignette chunk `simulate` ----------
truth  <- kernel_scale.raster(land_rast[[c("land1", "land2")]],
                              sigma = c(250, 500), kernel = "gaussian",
                              verbose = FALSE)
x_true <- scale(as.matrix(extract(truth, pts, ID = FALSE)))

set.seed(1)
d_mat <- as.matrix(dist(xy))
field <- drop(t(chol(0.6^2 * exp(-d_mat / 800))) %*% rnorm(nrow(xy)))
site  <- rnorm(nrow(xy))
log_mu <- 1 + 0.3 * site - 0.5 * x_true[, "land1"] +
  0.7 * x_true[, "land2"] + field
dat <- data.frame(counts = rpois(nrow(xy), exp(log_mu)), site = site)

# ---- Same preparation and starting model as the vignette -------------------
inputs <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
                      kernel = "gaussian", bin = TRUE, verbose = FALSE)
df <- data.frame(dat, inputs$kernel_dat)
df$pos <- glmmTMB::numFactor(xy[, 1], xy[, 2])
df$grp <- factor(rep("all", nrow(df)))

start_field <- glmmTMB::glmmTMB(
  counts ~ site + land1 + land2 + exp(pos + 0 | grp),
  family = poisson(), data = df,
  start  = list(theta = c(log(0.5), log(500)))
)

# Same call as the vignette chunk `fit-field`
t0 <- Sys.time()
fit_field <- multiScale_optim(start_field, inputs, n_cores = 1,
                              verbose = FALSE)
fit_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Drop the large duplicated pieces. The vignette restores kernel_inputs.
fit_field$kernel_inputs       <- NULL
fit_field$opt_context         <- NULL
fit_field$fitted_mod_original <- NULL

# ---- The bundled counts, which were simulated without a spatial field ------
# The vignette reports these numbers in one sentence to show what a field
# estimates when the data hold no residual spatial pattern.
dfr <- data.frame(landscape_counts, inputs$kernel_dat)
dfr$pos <- df$pos
dfr$grp <- df$grp
run_glm <- multiScale_optim(
  glm(counts ~ site + land1 + land2, family = poisson(), data = dfr),
  inputs, n_cores = 1, verbose = FALSE)
run_field <- multiScale_optim(
  glmmTMB::glmmTMB(counts ~ site + land1 + land2 + exp(pos + 0 | grp),
                   family = poisson(), data = dfr,
                   start = list(theta = c(log(0.5), log(500)))),
  inputs, n_cores = 1, verbose = FALSE)
bundled <- list(
  sigma_glm   = run_glm$scale_est,
  sigma_field = run_field$scale_est,
  field_par   = exp(glmmTMB::getME(run_field$opt_mod, "theta"))  # SD, range (m)
)

cache <- list(
  fit_field   = fit_field,
  dat         = dat,
  bundled     = bundled,
  fit_seconds = fit_seconds,
  versions    = list(multiScaleR = as.character(packageVersion("multiScaleR")),
                     glmmTMB = as.character(packageVersion("glmmTMB")),
                     TMB = as.character(packageVersion("TMB"))),
  created     = Sys.time()
)

saveRDS(cache, file.path("inst", "extdata",
                         "spatial-response-fields-cache.rds"),
        compress = "xz")
