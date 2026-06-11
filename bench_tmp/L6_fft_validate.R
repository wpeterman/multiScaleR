## Validate the FFT raster projection of the new metrics against landscapemetrics
## computed on the circular buffer (NA outside) at several interior cells.
.libPaths(c("C:/Temp/msrlib", .libPaths()))
suppressPackageStartupMessages({ library(multiScaleR); library(terra); library(landscapemetrics) })
options(width = 130)
adj_fft  <- multiScaleR:::.landscape_adjacency_raster_fft
comp_fft <- multiScaleR:::.landscape_composition_raster_fft

dim <- 40; radius <- 8
set.seed(3)
r <- terra::rast(nrows = dim, ncols = dim, xmin = 0, xmax = dim, ymin = 0, ymax = dim, crs = "EPSG:3857")
terra::values(r) <- sample(1:4, dim * dim, replace = TRUE)
M <- terra::as.matrix(r, wide = TRUE)
mkwin <- function(rr, cc) {  # masked circular window as a standalone raster
  ri <- (rr - radius):(rr + radius); ci <- (cc - radius):(cc + radius)
  sub <- M[ri, ci]; dmat <- sqrt(outer((ri - rr)^2, (ci - cc)^2, "+")); sub[dmat > radius] <- NA
  storage.mode(sub) <- "double"
  w <- terra::rast(sub); terra::ext(w) <- c(0, ncol(sub), 0, nrow(sub)); terra::crs(w) <- "EPSG:3857"; w
}
pts <- list(c(20, 20), c(15, 25), c(28, 12), c(25, 30), c(12, 18))

## landscape-level metrics: project once, compare extracted cell value to lm-on-window
cat("=== FFT projection vs landscapemetrics on the circular window ===\n")
lvl <- list(iji = lsm_l_iji, ent = lsm_l_ent, condent = lsm_l_condent,
            joinent = lsm_l_joinent, mutinf = lsm_l_mutinf, relmutinf = lsm_l_relmutinf)
cat(sprintf("%-10s %12s\n", "metric", "max|fft-lm|"))
for (nm in names(lvl)) {
  pr <- terra::as.matrix(adj_fft(r, radius = radius, base = 2, metric = nm), wide = TRUE)
  md <- 0
  for (p in pts) {
    fftv <- pr[p[1], p[2]]
    lmv <- suppressWarnings(lvl[[nm]](mkwin(p[1], p[2]))$value)
    md <- max(md, abs(fftv - lmv), na.rm = TRUE)
  }
  cat(sprintf("%-10s %12.3e\n", nm, md))
}

## class-level: clumpy / pland / ca for each present class at each point
cat("\n=== class-level FFT projection vs landscapemetrics on the window ===\n")
cl <- pl <- ca <- 0
for (p in pts) {
  w <- mkwin(p[1], p[2]); lmc <- suppressWarnings(lsm_c_clumpy(w))
  lmp <- suppressWarnings(lsm_c_pland(w)); lmca <- suppressWarnings(lsm_c_ca(w))
  for (k in 1:4) {
    cfft <- terra::as.matrix(adj_fft(r, radius = radius, metric = "clumpy", focal_class = k), wide = TRUE)[p[1], p[2]]
    pfft <- terra::as.matrix(comp_fft(r, radius = radius, metric = "pland", focal_class = k), wide = TRUE)[p[1], p[2]]
    afft <- terra::as.matrix(comp_fft(r, radius = radius, metric = "ca", focal_class = k), wide = TRUE)[p[1], p[2]]
    lc <- lmc$value[lmc$class == k]; if (length(lc) == 0) lc <- NA
    lp <- lmp$value[lmp$class == k]; if (length(lp) == 0) lp <- 0
    la <- lmca$value[lmca$class == k]; if (length(la) == 0) la <- 0
    cl <- max(cl, abs(cfft - lc), na.rm = TRUE); pl <- max(pl, abs(pfft - lp), na.rm = TRUE)
    ca <- max(ca, abs(afft - la), na.rm = TRUE)
  }
}
cat(sprintf("clumpy max|fft-lm| = %.3e\npland  max|fft-lm| = %.3e\nca     max|fft-lm| = %.3e\n", cl, pl, ca))
cat("L6 DONE\n")
