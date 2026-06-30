# Generate cached objects used by CRAN-friendly vignettes.
#
# Run from the package root after code changes that affect vignette output:
#   Rscript data-raw/generate_vignette_cache.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package `pkgload` is required.", call. = FALSE)
}

pkgload::load_all(".", helpers = FALSE, quiet = TRUE)

library(terra)

pkg_extdata <- function(...) {
  file.path("inst", "extdata", ...)
}

pack_plot_raster <- function(x, fact = 4) {
  terra::wrap(terra::aggregate(x, fact = fact, fun = mean, na.rm = TRUE))
}

strip_prep_payload <- function(x) {
  x$kernel_inputs <- NULL
  x
}

data("landscape_counts", package = "multiScaleR")
data("surv_pts", package = "multiScaleR")

dat <- landscape_counts
pts <- vect(surv_pts)
land_rast <- rast(pkg_extdata("landscape.tif"))

kernel_inputs <- kernel_prep(
  pts = pts,
  raster_stack = land_rast,
  max_D = 1700,
  kernel = "gaussian",
  store_cell_data = FALSE,
  verbose = FALSE
)

df <- data.frame(dat, kernel_inputs$kernel_dat)
mod0 <- glm(counts ~ site + land1 + land2 + land3,
            family = poisson(),
            data = df)
opt1 <- multiScale_optim(fitted_mod = mod0,
                         kernel_inputs = kernel_inputs,
                         verbose = FALSE)

mod0_2 <- glm(counts ~ site + land1 + land2,
              family = poisson(),
              data = df)
opt2 <- multiScale_optim(fitted_mod = mod0_2,
                         kernel_inputs = kernel_inputs,
                         verbose = FALSE)

prof2 <- profile_sigma(opt2, n_pts = 15, verbose = FALSE)
rast_opt <- kernel_scale.raster(raster_stack = land_rast,
                                multiScaleR = opt2,
                                verbose = FALSE)
r_scaled <- kernel_scale.raster(raster_stack = land_rast,
                                multiScaleR = opt2,
                                scale_center = TRUE,
                                clamp = TRUE,
                                verbose = FALSE)
pred <- terra::predict(r_scaled, opt2$opt_mod, type = "response")

landcover <- classify(
  land_rast$land2,
  rcl = matrix(
    c(-Inf, -0.5, 1,
      -0.5,  0.5, 2,
      0.5,  Inf, 3),
    ncol = 3,
    byrow = TRUE
  )
)
names(landcover) <- "landcover"
metric_rasters <- c(land_rast$land1, landcover)
landscape_vars <- msr_vars(
  land1_prop = kernel_var("land1"),
  land1_ed_500 = landscape_var("land1", metric = "ed", radius = 500),
  landcover_shdi_500 = landscape_var("landcover", metric = "shdi", radius = 500)
)
landscape_inputs <- kernel_prep(
  pts = pts,
  raster_stack = metric_rasters,
  max_D = 1200,
  kernel = "gaussian",
  scale_vars = landscape_vars,
  store_cell_data = FALSE,
  verbose = FALSE
)
landscape_df <- data.frame(dat, landscape_inputs$kernel_dat)
landscape_mod <- glm(
  counts ~ site + land1_prop + land1_ed_500 + landcover_shdi_500,
  family = poisson(),
  data = landscape_df
)
landscape_opt <- multiScale_optim(
  fitted_mod = landscape_mod,
  kernel_inputs = landscape_inputs,
  verbose = FALSE
)
landscape_projected <- kernel_scale.raster(
  raster_stack = metric_rasters,
  multiScaleR = landscape_opt,
  scale_center = TRUE,
  clamp = TRUE,
  verbose = FALSE
)

surface_vars <- msr_vars(
  land2_mean = kernel_var("land2"),
  land3_sq_500 = surface_var("land3", metric = "sq", radius = 500),
  land3_sa_500 = surface_var("land3", metric = "sa", radius = 500)
)
surface_inputs <- kernel_prep(
  pts = pts,
  raster_stack = land_rast,
  max_D = 1500,
  kernel = "gaussian",
  scale_vars = surface_vars,
  store_cell_data = FALSE,
  verbose = FALSE
)
surface_df <- data.frame(dat, surface_inputs$kernel_dat)
surface_mod <- glm(
  counts ~ site + land2_mean + land3_sq_500 + land3_sa_500,
  family = poisson(),
  data = surface_df
)
surface_opt <- multiScale_optim(
  fitted_mod = surface_mod,
  kernel_inputs = surface_inputs,
  verbose = FALSE
)
surface_projected <- kernel_scale.raster(
  raster_stack = land_rast,
  multiScaleR = surface_opt,
  scale_center = TRUE,
  clamp = TRUE,
  verbose = FALSE
)

cache <- list(
  guide = list(
    opt1 = opt1,
    opt2 = opt2,
    prof2 = prof2,
    rast_opt = pack_plot_raster(rast_opt)
  ),
  quickstart = list(
    opt = opt2,
    prof = prof2,
    r_scaled = pack_plot_raster(r_scaled),
    pred = pack_plot_raster(pred)
  ),
  landscape = list(
    opt = strip_prep_payload(landscape_opt),
    projected = pack_plot_raster(landscape_projected)
  ),
  surface = list(
    opt = strip_prep_payload(surface_opt),
    projected = pack_plot_raster(surface_projected)
  )
)

saveRDS(cache, pkg_extdata("vignette_cache.rds"), compress = "xz")
