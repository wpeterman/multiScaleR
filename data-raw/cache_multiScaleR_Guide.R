# Stored fits for vignettes/multiScaleR_Guide.Rmd.
#
# The zero-inflated and occupancy examples each take several seconds to
# simulate, prepare, and fit. Together they add about 20 s to the vignette
# build, so their fitted objects are stored and the vignette shows the calls
# with eval = FALSE.
#
# The code between the "must match" markers is copied from the vignette. Keep
# the two in step: same object names, arguments, and seeds.
#
# Run from the package root with the current build of the package installed:
#   Rscript data-raw/cache_multiScaleR_Guide.R
# Output: inst/extdata/multiScaleR_Guide-cache.rds (xz compressed)

library(multiScaleR)
library(terra)
library(unmarked)

# Drop what summary(), coef(), and plot() do not need: the optimization
# context, the starting model, and the per-cell and binned landscape data.
strip_fit <- function(x) {
  x$opt_context         <- NULL
  x$fitted_mod_original <- NULL
  x$kernel_inputs$binned  <- NULL
  x$kernel_inputs$d_list  <- NULL
  x$kernel_inputs$raw_cov <- NULL
  x
}

## ---- must match the vignette: zero-inflated counts ----
set.seed(12345)
zi_rast <- sim_rast(user_seed = 999, dim = 500, resolution = 30)
zi_rast <- subset(zi_rast, c("bin1", "cont1"))

interior <- buffer(vect(ext(zi_rast)), -1000)
zi_pts   <- spatSample(zi_rast[[1]], size = 400, method = "random",
                       ext = ext(interior), as.points = TRUE)

truth <- kernel_prep(pts = zi_pts, raster_stack = zi_rast, max_D = 1500,
                     kernel = "gaussian", sigma = c(400, 200),
                     verbose = FALSE)$kernel_dat

excess_zero <- rbinom(400, size = 1, prob = plogis(-0.25 - 1 * truth$cont1))
abundance   <- rnbinom(400, mu = exp(0.5 + 0.75 * truth$bin1), size = 2)
zi_dat      <- data.frame(cnt = ifelse(excess_zero == 1, 0, abundance))

zi_inputs <- kernel_prep(pts = zi_pts, raster_stack = zi_rast, max_D = 1500,
                         kernel = "gaussian", verbose = FALSE)
zi_df <- data.frame(zi_dat, zi_inputs$kernel_dat)

zinb_mod <- pscl::zeroinfl(cnt ~ bin1 | cont1, dist = "negbin", data = zi_df)
zinb_opt <- multiScale_optim(zinb_mod, zi_inputs, verbose = FALSE)
## ---- end zero-inflated counts ----

## ---- must match the vignette: occupancy ----
um_rast <- subset(sim_rast(user_seed = 999, dim = 250), c("bin1", "cont2"))

occ_sim <- sim_dat_unmarked(alpha = 0.25, beta = c(-0.75, 1),
                            sigma = c(225, 75),
                            kernel = "gaussian", type = "occ",
                            n_points = 250, n_surv = 5, det = 0.5,
                            raster_stack = um_rast, max_D = 800,
                            user_seed = 555)

occ_inputs <- kernel_prep(pts = occ_sim$pts, raster_stack = um_rast,
                          max_D = 800, kernel = "gaussian", verbose = FALSE)
occ_umf <- unmarkedFrameOccu(y = occ_sim$y, siteCovs = occ_inputs$kernel_dat)

occ_mod <- occu(~ 1 ~ bin1 + cont2, data = occ_umf)
opt_occ <- multiScale_optim(occ_mod, occ_inputs, verbose = FALSE)
## ---- end occupancy ----

cache <- list(
  zinb_opt = strip_fit(zinb_opt),
  opt_occ  = strip_fit(opt_occ),
  built_with = as.character(packageVersion("multiScaleR"))
)

out <- file.path("inst", "extdata", "multiScaleR_Guide-cache.rds")
saveRDS(cache, out, compress = "xz")
message("Saved ", out, " (", format(file.size(out), big.mark = ","), " bytes)")
