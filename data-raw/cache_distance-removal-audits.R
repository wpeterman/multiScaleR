# Stored results for vignettes/distance-removal-audits.Rmd
#
# The vignette's spatial abundance field section shows its calls with
# eval = FALSE, because fit_gdr_response_field() runs only with unmarked 1.5.2
# and TMB 1.9.25. This script runs those calls and stores the small table the
# vignette prints. It takes about 15 seconds.
#
# Run from the package root with the development build installed or loaded:
#   Rscript data-raw/cache_distance-removal-audits.R
# Output: inst/extdata/distance-removal-audits-cache.rds (xz compressed)

library(multiScaleR)
library(terra)

stopifnot(as.character(utils::packageVersion("unmarked")) == "1.5.2",
          as.character(utils::packageVersion("TMB")) == "1.9.25")

# ---- Landscape and survey points (same inputs as the vignette) -------------
data("surv_pts")
pts <- vect(surv_pts)
land_rast <- rast(system.file("extdata", "landscape.tif",
                              package = "multiScaleR"))

# Kernel-weighted habitat at a 250 m scale, standardized, one row per point
inputs_250 <- kernel_prep(pts = pts, raster_stack = land_rast[["land1"]],
                          max_D = 1700, kernel = "gaussian", sigma = 250,
                          verbose = FALSE)

# ---- Sixteen square blocks over the extent of the survey points ------------
xy <- crds(pts)
n_side <- 4
width <- diff(range(xy[, 1])) / n_side
height <- diff(range(xy[, 2])) / n_side
col_id <- pmin(n_side - 1, floor((xy[, 1] - min(xy[, 1])) / width))
row_id <- pmin(n_side - 1, floor((xy[, 2] - min(xy[, 2])) / height))
block <- factor(paste0("b", col_id, row_id))

# Every survey in a block carries the coordinates of the block center
block_xy <- cbind(x = min(xy[, 1]) + (col_id + 0.5) * width,
                  y = min(xy[, 2]) + (row_id + 0.5) * height)

# ---- Block effects with an exponential correlation range of 1500 m ---------
field_range <- 1500
field_sd <- 0.5
centers <- block_xy[!duplicated(block), , drop = FALSE]
corr <- exp(-as.matrix(dist(centers)) / field_range)
set.seed(3)
effect <- drop(t(chol(corr)) %*% rnorm(nrow(centers))) * field_sd
block_effect <- effect[match(block, block[!duplicated(block)])]

# ---- Simulated surveys: density depends on habitat and the block effect ----
sim_blocks <- sim_dat_gdr(
  density          = exp(log(400) + 0.5 * inputs_250$kernel_dat$land1 +
                           block_effect),
  distance_scale   = 40,
  removal          = 0.2,
  dist_breaks      = c(0, 25, 50, 75),
  interval_lengths = c(2, 3, 5),
  interval_units   = "minutes",
  site_covs        = data.frame(inputs_250$kernel_dat, block = block),
  seed             = 102
)
umf_blocks <- sim_blocks$observed
rownames(block_xy) <- rownames(unmarked::siteCovs(umf_blocks))

# ---- The calls shown in the vignette ----------------------------------------
# 1. Template: independent random intercept for each block (TMB engine)
template <- unmarked::gdistremoval(
  lambdaformula   = ~ land1 + (1 | block),
  removalformula  = ~ 1,
  distanceformula = ~ 1,
  data = umf_blocks, keyfun = "halfnorm", mixture = "P",
  output = "density", unitsOut = "kmsq", K = 80, engine = "TMB",
  starts = c(log(400), 0, log(40), qlogis(0.2), log(0.5))
)

# 2. A second, independent set of starting values, named by parameter block
cold <- c("beta_lambda[1]"   = log(300),    # log density at average habitat
          "beta_lambda[2]"   = 0,           # habitat effect
          "beta_dist[1]"     = log(30),     # log distance scale (m)
          "beta_rem[1]"      = qlogis(0.4), # logit removal probability
          "lsigma_lambda[1]" = log(0.3))    # log SD of the block effects
settings <- list(method = "BFGS", se = TRUE, threads = 1L,
                 control = list(maxit = 2500L, reltol = 1e-10))

# 3. Replace the independent block effects with a spatially correlated field
elapsed <- system.time(
  field_fit <- fit_gdr_response_field(
    model = template, coords = block_xy, crs = 32617, range = field_range,
    cold_starts = list(cold = cold),
    audit_dir = tempfile("gdr-field-audit-"),
    tmb_settings = settings
  )
)[["elapsed"]]

# ---- Small tables for the vignette ------------------------------------------
block_sd <- function(fit) {
  # SD of the block effects, on the log density scale
  as.numeric(unmarked::sigma(fit)$sigma[1])
}
row_for <- function(fit, label) {
  est <- unmarked::coef(fit)
  data.frame(
    model        = label,
    density      = unname(exp(est["lambda(Int)"])),
    land1_effect = unname(est["lambda(land1)"]),
    block_sd     = block_sd(fit),
    logLik       = as.numeric(unmarked::logLik(fit)),
    AIC          = fit@AIC
  )
}
compare <- rbind(row_for(template, "Independent blocks"),
                 row_for(field_fit, "Correlated field"))

cache <- list(
  field = list(
    compare     = compare,
    n_blocks    = nlevels(block),
    range       = field_range,
    class       = class(field_fit)[1],
    seconds     = round(elapsed),
    versions    = c(unmarked = as.character(utils::packageVersion("unmarked")),
                    TMB = as.character(utils::packageVersion("TMB")),
                    multiScaleR = as.character(utils::packageVersion("multiScaleR")))
  )
)

saveRDS(cache, file.path("inst", "extdata", "distance-removal-audits-cache.rds"),
        compress = "xz")
print(cache$field$compare)
cat("bytes:", file.size(file.path("inst", "extdata",
                                  "distance-removal-audits-cache.rds")), "\n")
