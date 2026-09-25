# Year-specific landscapes in one multiScaleR model
# Run with library(multiScaleR), or devtools::load_all() for a source checkout.
# Each point is extracted only from its own year's map.
library(multiScaleR)
library(terra)

# 1. Create two aligned, projected habitat maps (1 = habitat, 0 = other).
# With real data, replace these maps with your annual SpatRaster objects.
# Do not average numeric land-cover class codes. Use binary class layers for
# composition, or explicit landscape_var() specifications for landscape metrics.
set.seed(47)
habitat_2024 <- rast(nrows = 100, ncols = 100, resolution = 20,
                     xmin = 0, xmax = 2000, ymin = 0, ymax = 2000,
                     crs = "EPSG:26915")
xy <- xyFromCell(habitat_2024, seq_len(ncell(habitat_2024)))
values(habitat_2024) <- as.integer(
  sin(xy[, 1] / 140) + cos(xy[, 2] / 190) > 0
)
habitat_2025 <- habitat_2024
managed <- xy[, 1] > 700 & xy[, 1] < 1300 &
  xy[, 2] > 700 & xy[, 2] < 1300
habitat_2025[which(managed)] <- 0
names(habitat_2024) <- names(habitat_2025) <- "habitat"
stopifnot(compareGeom(habitat_2024, habitat_2025),
          identical(names(habitat_2024), names(habitat_2025)))

# One row per observation, with unique IDs even for repeated locations.
# Keep every buffer inside the raster extent and use consistent map units.
obs <- data.frame(
  year = factor(rep(c("2024", "2025"), 100)),
  easting = runif(200, 420, 1580),
  northing = runif(200, 420, 1580)
)
rownames(obs) <- paste0("nest_", seq_len(nrow(obs)))
pts <- sf::st_as_sf(obs, coords = c("easting", "northing"), crs = 26915)

# 2. Prepare one pooled object. The function applies each year's raster to
# every nest from that year, including nests outside management polygons whose
# landscape buffers intersect changed cells. Bins use all observations.
pooled <- kernel_prep_by_group(
  pts, raster_stacks = list("2024" = habitat_2024,
                            "2025" = habitat_2025),
  group = obs$year, max_D = 400, sigma = 100,
  kernel = "gaussian", bin = TRUE, nbins = 256L,
  store_cell_data = FALSE, verbose = FALSE
)
stopifnot(identical(rownames(pooled$kernel_dat), rownames(obs)),
          identical(unname(pooled$raster_group), as.character(obs$year)))

# 3. Create an illustrative binary response and fit one pooled model.
# This response represents known survival over the SAME observation interval.
# Actual nest/juvenile analyses need a likelihood appropriate to exposure time,
# censoring, detection, and repeated observations. This GLM is a wiring example.
# Replace this simulation with your observed response and suitable model.
obs$survived <- rbinom(nrow(obs), 1,
  plogis(-0.3 + 1.2 * pooled$kernel_dat$habitat +
           0.4 * (obs$year == "2025")))
model_data <- cbind(obs, pooled$kernel_dat)
initial_model <- glm(survived ~ habitat + year, data = model_data,
                     family = binomial(), na.action = na.fail)
fit <- multiScale_optim(initial_model, kernel_inputs = pooled,
                        n_cores = 1, verbose = FALSE)

# One habitat scale and coefficient are estimated from both years together.
# The year term allows baseline survival to differ between years. Omitting a
# habitat-by-year interaction assumes a shared slope; it does not establish it.
print(fit$scale_est)          # Gaussian sigma in map units (here, meters).
print(coef(fit$opt_mod))      # Log-odds coefficients; habitat is standardized.
print(fit$diagnostics)       # Inspect convergence and scale-boundary diagnostics.

# Keep both maps on the same grid with the same names and class coding.
# Missing cells may differ across years; the bin summaries track valid counts.
# Configuration metrics require their own scale_vars definitions in BOTH calls;
# the binary habitat mean above measures composition, not configuration.
# For landscape/surface metrics, set store_cell_data = TRUE.
