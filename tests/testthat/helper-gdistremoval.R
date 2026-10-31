gdr_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(715)
    raster <- terra::rast(nrows = 24, ncols = 24, xmin = 0, xmax = 480,
      ymin = 0, ymax = 480, crs = "EPSG:26917")
    xy <- terra::xyFromCell(raster, seq_len(terra::ncell(raster)))
    terra::values(raster) <- sin(xy[, 1] / 31) + cos(xy[, 2] / 43)
    names(raster) <- "hab"
    pts <- expand.grid(x = seq(90, 390, length.out = 6), y = seq(90, 390, length.out = 4))
    points <- sf::st_as_sf(pts, coords = c("x", "y"), crs = 26917)
    inputs <- kernel_prep(points, raster, max_D = 80, sigma = 30,
      kernel = "gaussian", bin = FALSE, verbose = FALSE)
    m <- nrow(pts); tt <- 2L; j <- 3L
    sc <- data.frame(site = seq_len(m), hab = inputs$kernel_dat$hab)
    annual <- data.frame(visit = rep(c(-0.5, 0.5), m))
    obs <- data.frame(tod = rep(c(-0.5, 0, 0.5), m * tt))
    abundance <- rpois(m, exp(2 + 0.4 * sc$hab))
    yd <- yr <- matrix(0, m, tt * j)
    for (site in seq_len(m)) for (t in seq_len(tt)) {
      count <- rbinom(1, abundance[site], 0.45)
      cols <- (t - 1L) * j + seq_len(j)
      yd[site, cols] <- as.vector(rmultinom(1, count, c(0.25, 0.45, 0.3)))
      yr[site, cols] <- as.vector(rmultinom(1, count, c(0.5, 0.3, 0.2)))
    }
    frame <- unmarked::unmarkedFrameGDR(yd, yr, numPrimary = tt, siteCovs = sc,
      yearlySiteCovs = annual, obsCovs = obs, dist.breaks = c(0, 25, 60, 100),
      unitsIn = "m", period.lengths = c(1, 2, 2))
    model <- unmarked::gdistremoval(~hab, ~1, ~tod, ~visit,
      data = frame, K = 35, se = FALSE)
    context <- build_opt_context(model, inputs$raw_cov, join_by = sc["site"],
      unit_conv = inputs$unit_conv)
    evaluated <- kernel_scale_fn(30 / inputs$unit_conv, inputs$d_list, inputs$raw_cov,
      inputs$kernel, model, opt_context = context, mod_return = TRUE)
    object <- structure(list(opt_mod = evaluated$mod, kernel_inputs = inputs,
      scl_params = evaluated$scl_params, join_by = sc["site"],
      scale_est = data.frame(Mean = 30, SE = NA_real_, row.names = "hab"),
      min_D = inputs$min_D, max_D = inputs$max_D), class = "multiScaleR")
    cache <<- list(model = model, frame = frame, inputs = inputs,
      object = object, raster = raster, points = points)
    cache
  }
})

# Known coefficients supply an independent analytical oracle. This modified
# object is used only for prediction/simulation, never for likelihood ranking.
gdr_oracle <- function(mixture = "P", output = "abund", keyfun = "halfnorm") {
  model <- gdr_fixture()$model
  values <- list(lambda = c(log(8), 0), phi = qlogis(0.7),
    rem = c(qlogis(0.3), 0), dist = c(log(50), log(2)))
  model@data@yearlySiteCovs$visit <- rep(c(0, 1), nrow(model@data@yDistance))
  for (nm in names(values)) model@estimates@estimates[[nm]]@estimates[] <- values[[nm]]
  model@mixture <- mixture
  model@output <- output
  model@keyfun <- keyfun
  if (mixture != "P") {
    alpha <- model@estimates@estimates$phi
    alpha@name <- "Dispersion"
    alpha@short.name <- "alpha"
    alpha@estimates[] <- if (mixture == "NB") log(2) else qlogis(0.4)
    alpha@invlink <- "exp"
    model@estimates@estimates$alpha <- alpha
  }
  if (keyfun == "hazard") {
    shape <- model@estimates@estimates$phi
    shape@name <- "Hazard-rate (scale)"
    shape@short.name <- "scale"
    shape@estimates[] <- log(2)
    shape@invlink <- "exp"
    model@estimates@estimates$scale <- shape
  }
  model
}
