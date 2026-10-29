# Reproduce the small artificial future landscape scenario used in
# vignette("bayesian-prediction"). No fitting, random sampling or file writing.
# This is a scenario derived from existing bin means, not a landscape forecast.
make_bayesian_prediction_scenario <- function(annual = NULL) {
  if (is.null(annual)) annual <- readRDS(system.file("extdata",
    "bayesian-annual-example.rds", package = "multiScaleR"))
  dat <- expand.grid(site = c("s1", "s2", "new_site_A", "new_site_B"),
    year = c("2026", "2027"), stringsAsFactors = FALSE)
  dat$site_year <- paste(dat$site, dat$year, sep = "_")
  # The last two positions are relabeled as unsurveyed model units for this
  # artificial demonstration. An actual forecast needs the new sites' maps.
  keys <- paste0(rep(c("s1", "s2", "s3", "s4"), 2), "_2025")
  rows <- match(keys, annual$inputs$binned$point_ids)
  stopifnot(!anyNA(rows))
  inputs <- annual$inputs; bins <- inputs$binned
  bins$vsum <- lapply(bins$vsum, function(x) x[rows, , drop = FALSE])
  bins$csum <- lapply(bins$csum, function(x) x[rows, , drop = FALSE])
  # Add a transparent smooth shift to each ring's mean, preserving cell counts.
  shift <- outer(rep(1:2, each = 4), sin(6 * bins$bin_rep),
    function(year, distance) .15 * year * distance)
  bins$vsum$habitat <- bins$vsum$habitat + shift * bins$csum$habitat
  bins$point_ids <- dat$site_year; bins$n_points <- nrow(dat)
  inputs$binned <- bins[c("bin_rep", "vsum", "csum", "point_ids", "n_points")]
  inputs <- inputs[intersect(names(inputs), c("binned", "min_D", "max_D",
    "unit_conv", "kernel", "scale_vars"))]
  list(newdata = dat, kernel_inputs = inputs)
}
