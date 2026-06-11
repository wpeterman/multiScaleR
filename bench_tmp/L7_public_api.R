## End-to-end smoke test through the public API + the raster projection path.
.libPaths(c("C:/Temp/msrlib", .libPaths()))
suppressPackageStartupMessages({ library(multiScaleR); library(terra) })
options(width = 130)

set.seed(11)
r <- terra::rast(nrows = 30, ncols = 30, xmin = 0, xmax = 30, ymin = 0, ymax = 30, crs = "EPSG:3857")
terra::values(r) <- sample(1:4, 900, replace = TRUE)
names(r) <- "lc"

## (1) build specs spanning all new families; check the data frame + print
vars <- msr_vars(
  juxta   = landscape_var("lc", metric = "iji",       radius = 6),
  entropy = landscape_var("lc", metric = "ent",       radius = 6, base = 2),
  relinfo = landscape_var("lc", metric = "relmutinf", radius = 6),
  clump2  = landscape_var("lc", metric = "clumpy",    radius = 6, class = 2),
  cover3  = landscape_var("lc", metric = "pland",     radius = 6, class = 3),
  area1   = landscape_var("lc", metric = "ca",        radius = 6, class = 1)
)
cat("=== (1) msr_vars data frame (note `class` column) ===\n")
print(vars)
stopifnot("class" %in% names(vars), vars$class[vars$covariate == "clump2"] == 2)

## (2) argument validation
cat("\n=== (2) class-argument validation ===\n")
chk <- function(label, expr) {
  res <- tryCatch({ expr; "NO ERROR" }, error = function(e) conditionMessage(e))
  cat(sprintf("  %-28s -> %s\n", label, res))
}
chk("clumpy without class", landscape_var("lc", "clumpy"))
chk("pland without class", landscape_var("lc", "pland"))
chk("iji with class (forbidden)", landscape_var("lc", "iji", class = 2))
chk("clumpy with class=2 (ok)", landscape_var("lc", "clumpy", class = 2))

## (3) project all six through the kernel_scale.raster code path
cat("\n=== (3) raster projection via .msr_scale_vars_raster (fixed radius) ===\n")
proj <- multiScaleR:::.msr_scale_vars_raster(
  raster_stack = r, scale_vars = vars, sigma = numeric(0), shape = NULL,
  kernel = "gaussian", pct_wt = 0.95, fft = TRUE, na.rm = TRUE, verbose = FALSE
)
print(terra::global(proj, fun = "mean", na.rm = TRUE))
cat("\nlayer names:", paste(names(proj), collapse = ", "), "\n")
cat("L7 DONE\n")
