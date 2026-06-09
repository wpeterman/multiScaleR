# Regenerate the precomputed `unmarked::pcount` object used by the
# "Optimization with unmarked" section of the multiScaleR_Guide vignette.
#
# The pcount optimization refits a K = 100 N-mixture model at every optimizer
# step, so it is the one computationally intensive fit in that vignette
# (~4-5 minutes). To keep the vignette build fast, its result is precomputed
# here, stripped to the elements needed for display + projection, and saved to
# inst/extdata/opt_umf_p.rds (loaded in the vignette via readRDS()).
#
# Re-run this whenever the simulation spec or the modeling code changes, so the
# vignette's printed pcount results stay current. Every OTHER model in the
# vignette is computed live at build time and does not need regeneration.
#
# Usage (from the package root):
#   source("data-raw/generate_pcount_vignette_object.R")

library(multiScaleR)
library(terra)
library(unmarked)

## --- This block must match the vignette's pcount example exactly ---
rs <- sim_rast(user_seed = 999, dim = 250)
rs <- terra::subset(rs, c(1, 4))
s_dat <- sim_dat_unmarked(alpha = 1, beta = c(0.75, -0.75), kernel = "gaussian",
                          sigma = c(75, 150), n_points = 75, n_surv = 5, det = 0.5,
                          type = "count", raster_stack = rs, max_D = 550, user_seed = 555)
kernel_inputs <- kernel_prep(pts = s_dat$pts, raster_stack = rs, max_D = 550,
                             kernel = "gaus", verbose = FALSE)
umf <- unmarkedFramePCount(y = s_dat$y, siteCovs = kernel_inputs$kernel_dat)
mod0_umf.p <- pcount(~1 ~bin1 + cont2, data = umf, K = 100)

opt_umf.p <- multiScale_optim(fitted_mod = mod0_umf.p, kernel_inputs = kernel_inputs)

## --- Strip to the elements needed for summary()/plot()/projection ---
## (drops the per-cell data, binned summaries, and the optimization context /
## original model, which are only needed to re-optimize or profile.)
opt_umf.p$opt_context         <- NULL
opt_umf.p$fitted_mod_original <- NULL
opt_umf.p$kernel_inputs$binned  <- NULL
opt_umf.p$kernel_inputs$d_list  <- NULL
opt_umf.p$kernel_inputs$raw_cov <- NULL

saveRDS(opt_umf.p, "inst/extdata/opt_umf_p.rds", version = 2)
message("Saved inst/extdata/opt_umf_p.rds (",
        format(file.info("inst/extdata/opt_umf_p.rds")$size, big.mark = ","), " bytes)")
