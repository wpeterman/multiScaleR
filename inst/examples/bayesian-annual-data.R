# Reproduce the artificial data used in vignette("bayesian-scales").
# The seed controls maps, site effects and NB2 responses. This script does not
# fit a model or write files. Distances are in projected map units (meters).
make_bayesian_annual_data <- function() {
  set.seed(521)
  r <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 6000,
    ymin = 0, ymax = 6000, crs = "EPSG:26915")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  years <- as.character(2021:2025)
  maps <- lapply(seq_along(years), function(j) {
    z <- r
    terra::values(z) <- sin(xy[, 1] / 400 + j / 4) + cos(xy[, 2] / 550) +
      rnorm(nrow(xy), sd = .2)
    names(z) <- "habitat";z
  })
  names(maps) <- years
  locations <- expand.grid(x = seq(1400, 4600, length.out = 5),
    y = seq(1400, 4600, length.out = 5))
  dat <- expand.grid(site = paste0("s", seq_len(nrow(locations))),
    year = years, stringsAsFactors = FALSE)
  site_index <- match(dat$site, paste0("s", seq_len(nrow(locations))))
  # The third and fifth years have fewer surveyed sites.
  keep <- (dat$year != "2023" | site_index %in% seq(1, 25, by = 3)) &
    (dat$year != "2025" | site_index <= 16)
  dat <- dat[keep, ];site_index <- site_index[keep]
  dat$site <- factor(dat$site);dat$year <- factor(dat$year, levels = years)
  dat$site_year <- paste(dat$site, dat$year, sep = "_")
  rownames(dat) <- NULL
  points <- sf::st_as_sf(locations[site_index, ], coords = c("x", "y"), crs = 26915)
  rownames(points) <- dat$site_year
  inputs <- multiScaleR::kernel_prep_by_group(points, maps,
    group = as.character(dat$year), max_D = 1200, sigma = 250, nbins = 48,
    store_cell_data = FALSE, verbose = FALSE)
  # Compute the generating predictor from the same finite-cell bin summaries.
  # Holding its reference transformation fixed matches the Bayesian contract.
  b <- inputs$binned;row <- match(dat$site_year, b$point_ids)
  values <- b$vsum$habitat[row, , drop = FALSE]
  counts <- b$csum$habitat[row, , drop = FALSE]
  kernel_mean <- function(distance) {
    w <- exp(-.5 * outer(1 / (distance / inputs$unit_conv), b$bin_rep)^2)
    rowSums(w * values) / rowSums(w * counts)
  }
  true_scale <- c(180, 250, 350, 260, 210)
  reference <- kernel_mean(rep(250, nrow(dat)))
  habitat <- (kernel_mean(true_scale[as.integer(dat$year)]) - mean(reference)) / sd(reference)
  site_effect <- rnorm(nrow(locations), sd = .25)
  year_effect <- c(-.1, .1, 0, -.1, .1)
  eta <- log(15) + year_effect[as.integer(dat$year)] + habitat + site_effect[site_index]
  dat$count <- rnbinom(nrow(dat), mu = exp(eta), size = 8)
  list(data = dat, inputs = inputs,
    truth = data.frame(year = years, scale = true_scale,
      observations = as.integer(table(dat$year))),
    generating = list(seed = 521L, population_scale = 250, scale_sd = .35,
      scale_z = log(true_scale / 250) / .35, response_sd = .25,
      response_intercept = log(15), habitat_coefficient = 1, dispersion = 8))
}
