# A small simulated step-selection data set for conditional_logit tests.
# Choices respond to forest smoothed with a Gaussian kernel, sigma = 80 m.
# Returns the step table, its points, the forest raster, kernel inputs at
# max_D = 300 m, and the step table joined to the kernel covariate.
ssf_fixture <- function() {
  set.seed(42)
  r <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 1200,
                   ymin = 0, ymax = 1200, crs = "EPSG:32617")
  terra::values(r) <- stats::rnorm(terra::ncell(r))
  forest <- terra::focal(r, w = terra::focalMat(r, 60, "Gauss"), fun = "sum",
                         na.rm = TRUE) > 0
  forest <- forest * 1
  names(forest) <- "forest"
  perceived <- terra::focal(forest, w = terra::focalMat(forest, 80, "Gauss"),
                            fun = "sum", na.rm = TRUE)
  n_sets <- 150; n_alt <- 4
  sx <- rep(stats::runif(n_sets, 400, 800), each = n_alt)
  sy <- rep(stats::runif(n_sets, 400, 800), each = n_alt)
  len <- stats::rgamma(n_sets * n_alt, shape = 2, scale = 50)
  ang <- stats::runif(n_sets * n_alt, 0, 2 * pi)
  steps <- data.frame(step_id = rep(seq_len(n_sets), each = n_alt),
                      x = sx + len * cos(ang), y = sy + len * sin(ang),
                      sl_km = len / 1000)
  f <- terra::extract(perceived, as.matrix(steps[, c("x", "y")]))[, 1]
  z <- (f - mean(f)) / stats::sd(f)
  score <- exp(1.2 * z - 2 * steps$sl_km)
  steps$case_ <- 0
  for (s in seq_len(n_sets)) {
    rows <- which(steps$step_id == s)
    steps$case_[sample(rows, 1, prob = score[rows])] <- 1
  }
  pts <- terra::vect(steps, geom = c("x", "y"), crs = "EPSG:32617")
  inputs <- kernel_prep(pts = pts, raster_stack = forest, max_D = 300,
                        kernel = "gaussian", verbose = FALSE)
  list(steps = steps, pts = pts, forest = forest, inputs = inputs,
       df = data.frame(steps, inputs$kernel_dat))
}

# Native clogit() log-likelihood on kernel values prepared at sigma `s`,
# an independent check on the backend's internal kernel engine.
ssf_native_loglik <- function(fx, formula, s) {
  coxph <- survival::coxph
  k <- kernel_prep(pts = fx$pts, raster_stack = fx$forest, max_D = 300,
                   kernel = "gaussian", sigma = s, verbose = FALSE)
  as.numeric(stats::logLik(survival::clogit(formula,
    data = data.frame(fx$steps, k$kernel_dat))))
}
