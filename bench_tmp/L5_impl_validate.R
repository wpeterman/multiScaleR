## Validate the IMPLEMENTED metrics (production C++ direct path + FFT raster path)
## against landscapemetrics. The C++ adjacency/composition kernels, given a buffer
## large enough to cover an entire small raster, return the whole-landscape value.
.libPaths(c("C:/Temp/msrlib", .libPaths()))
suppressPackageStartupMessages({ library(multiScaleR); library(terra); library(landscapemetrics) })
options(width = 130)
.landscape_adjacency_by_buffer <- multiScaleR:::.landscape_adjacency_by_buffer
.landscape_composition_by_buffer <- multiScaleR:::.landscape_composition_by_buffer

mk <- function(seed, dim = 14, nclass = 4) { set.seed(seed)
  r <- terra::rast(nrows = dim, ncols = dim, xmin = 0, xmax = dim, ymin = 0, ymax = dim, crs = "EPSG:3857")
  terra::values(r) <- sample(seq_len(nclass), terra::ncell(r), replace = TRUE); r }

## ---- (1) direct C++ path vs landscapemetrics on the full raster ----
cat("=== (1) direct C++ path vs landscapemetrics (whole raster; base=2 for info theory) ===\n")
land_lvl <- list(
  iji       = lsm_l_iji,    ent     = lsm_l_ent,     condent = lsm_l_condent,
  joinent   = lsm_l_joinent, mutinf  = lsm_l_mutinf,  relmutinf = lsm_l_relmutinf,
  contag    = lsm_l_contag, ai = lsm_l_ai, pladj = lsm_l_pladj
)
adj_cpp <- function(r, metric, base = 2, focal = NULL) {
  v <- terra::values(r)[, 1]; n <- length(v)
  .landscape_adjacency_by_buffer(d = rep(0, n), r_stack.df = data.frame(x = v),
    cells = seq_len(n), radius = 1e6, n_cols = terra::ncol(r),
    metric = metric, base = base, focal_class = focal, validate = FALSE)[[1]]
}
comp_cpp <- function(r, metric, focal = NULL) {
  v <- terra::values(r)[, 1]; n <- length(v)
  .landscape_composition_by_buffer(d = rep(0, n), r_stack.df = data.frame(x = v),
    radius = 1e6, metric = metric, resolution = terra::res(r)[1],
    focal_class = focal, validate = FALSE)[[1]]
}
maxdiff <- list()
for (seed in c(1, 4, 9)) {
  r <- mk(seed)
  for (nm in names(land_lvl)) {
    got <- adj_cpp(r, nm); ref <- land_lvl[[nm]](r)$value
    maxdiff[[nm]] <- max(maxdiff[[nm]], abs(got - ref), na.rm = TRUE)
  }
}
for (nm in names(land_lvl))
  cat(sprintf("  %-10s max|cpp - lm| = %.3e\n", nm, maxdiff[[nm]]))

## class-level: clumpy / pland / ca per class
cat("\n  class-level (clumpy / pland / ca), max over classes & seeds:\n")
cd <- pd <- cad <- 0
for (seed in c(1, 4, 9)) {
  r <- mk(seed); res <- terra::res(r)[1]
  lmc <- lsm_c_clumpy(r); lmp <- lsm_c_pland(r); lmca <- lsm_c_ca(r)
  for (k in sort(unique(terra::values(r)[, 1]))) {
    cd  <- max(cd,  abs(adj_cpp(r, "clumpy", focal = k) - lmc$value[lmc$class == k]))
    pd  <- max(pd,  abs(comp_cpp(r, "pland", focal = k) - lmp$value[lmp$class == k]))
    cad <- max(cad, abs(comp_cpp(r, "ca",    focal = k) - lmca$value[lmca$class == k]))
  }
}
cat(sprintf("  %-10s max|cpp - lm| = %.3e\n", "clumpy", cd))
cat(sprintf("  %-10s max|cpp - lm| = %.3e\n", "pland", pd))
cat(sprintf("  %-10s max|cpp - lm| = %.3e\n", "ca", cad))
cat("L5 DONE\n")
