# Functional check: fold helpers into cross_validate(), and spline formulas
# through fit / predict / effect curves / CV. Scratch library only.
scr <- "C:/Users/peterman.73/AppData/Local/Temp/claude/9e3b7d3898bd820b/459cb932-8333-43df-b0da-7d9a2422aa59/scratchpad"
.libPaths(c(file.path(scr, "rlib"), .libPaths()))
suppressPackageStartupMessages({library(multiScaleR); library(terra); library(splines)})
try_ <- function(label, expr) {
  t0 <- Sys.time()
  out <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  [warning]", conditionMessage(w), "\n"); invokeRestart("muffleWarning")}),
    error = function(e) { cat("  [ERROR]", conditionMessage(e), "\n"); NULL })
  cat(sprintf("== %s: %s (%.1fs)\n", label, if (is.null(out)) "FAILED" else "ok",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(out)
}
data("landscape_counts"); data("surv_pts")
land <- terra::rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))
pts <- surv_pts
inputs <- kernel_prep(pts = pts, raster_stack = land, max_D = 1700,
                      kernel = "gaussian", bin = TRUE, verbose = FALSE)
dat <- data.frame(landscape_counts, inputs$kernel_dat)
ids <- rownames(dat); cat("n =", nrow(dat), " id head:", head(ids, 3), "\n")

cat("\n######## FOLD HELPERS\n")
sf_folds <- try_("spatial_cv_folds(sf points)", spatial_cv_folds(as.data.frame(sf::st_coordinates(pts)), block_size = 1000, k = 3, crs = 32617))
if (!is.null(sf_folds)) print(sf_folds)
tv_folds <- try_("spatial_cv_folds(SpatVector)", spatial_cv_folds(terra::vect(sf::st_set_crs(pts, 32617)), block_size = 1000, k = 3))
buf_folds <- try_("spatial_cv_folds(buffer = 200)", spatial_cv_folds(sf::st_set_crs(pts, 32617), block_size = 1000, k = 3, buffer = 200))
if (!is.null(buf_folds)) print(lengths(attr(buf_folds, "fold_info")$excluded))
grp <- try_("grouped_cv_folds", grouped_cv_folds(rep(seq_len(20), length.out = nrow(dat)), k = 4, ids = ids))
cv1 <- try_("cross_validate(profile_glm, spatial folds)",
  cross_validate(counts ~ site + land1 + land2, dat, inputs, sf_folds, family = "poisson"))
if (!is.null(cv1)) { cat("complete:", cv1$complete, "\n"); print(cv1$fold_results) }
cv0 <- try_("cross_validate_glm null", cross_validate_glm(counts ~ 1, dat, sf_folds, family = "poisson"))
if (!is.null(cv1) && !is.null(cv0)) print(try_("cv_tab", cv_tab(list(spatial = cv1, null = cv0))))
cv2 <- try_("cross_validate(refit, grouped folds)",
  cross_validate(counts ~ site + land1 + land2, dat, inputs, grp, family = "poisson", backend = "refit"))
if (!is.null(cv2)) cat("complete:", cv2$complete, "\n")

cat("\n######## SPLINES\n")
k1 <- 0; bk <- c(-3, 3)   # knots on the standardized kernel covariate
f_fixed <- counts ~ site + ns(land1, knots = 0, Boundary.knots = c(-3, 3)) + land2
f_df    <- counts ~ site + ns(land1, df = 3) + land2
m_fixed <- glm(f_fixed, family = poisson(), data = dat)
m_df    <- glm(f_df, family = poisson(), data = dat)
s1 <- try_("refit backend, fixed knots", multiScale_optim(m_fixed, inputs, n_cores = 1, verbose = FALSE))
s2 <- try_("refit backend, adaptive df", multiScale_optim(m_df, inputs, n_cores = 1, verbose = FALSE))
s3 <- try_("profile_glm, fixed knots", multiScale_optim(m_fixed, inputs, backend = "profile_glm", n_cores = 1, verbose = FALSE))
s4 <- try_("profile_glm, adaptive df", multiScale_optim(m_df, inputs, backend = "profile_glm", n_cores = 1, verbose = FALSE))
for (nm in c("s1", "s2", "s3", "s4")) { o <- get(nm); if (!is.null(o)) { cat(nm, "\n"); print(o$scale_est) } }
grDevices::pdf(NULL)
for (nm in c("s1", "s2", "s3")) {
  o <- get(nm); if (is.null(o)) next
  try_(paste(nm, "summary"), utils::capture.output(print(summary(o))))
  pr <- try_(paste(nm, "predict()"), predict(o))
  p <- try_(paste(nm, "plot_marginal_effects"), plot_marginal_effects(o))
  if (!is.null(p)) cat("  panels:", names(p), " caption null:", is.null(p[[1]]$labels$caption), "\n")
  try_(paste(nm, "plot()"), plot(o))
  try_(paste(nm, "kernel_scale.raster"), kernel_scale.raster(land, multiScaleR = o, verbose = FALSE))
  try_(paste(nm, "aic_tab"), aic_tab(list(o)))
}
if (!is.null(s3)) {
  try_("s3 predict joint", predict(s3, interval = "confidence", uncertainty = "joint"))
  try_("s3 profile summary", utils::capture.output(print(summary(s3, profile = TRUE))))
}
cvs <- try_("cross_validate spline fixed (profile_glm)",
  cross_validate(f_fixed, dat, inputs, sf_folds, family = "poisson"))
if (!is.null(cvs)) { cat("complete:", cvs$complete, "\n"); print(cvs$fold_results); if (!cvs$complete) print(cvs$errors) }
cvd <- try_("cross_validate spline df (refit)",
  cross_validate(f_df, dat, inputs, sf_folds, family = "poisson", backend = "refit"))
if (!is.null(cvd)) { cat("complete:", cvd$complete, "\n"); if (!cvd$complete) print(cvd$errors) }
cvg <- try_("cross_validate_glm spline on site", cross_validate_glm(counts ~ ns(site, df = 3), dat, sf_folds, family = "poisson"))
if (!is.null(cvg)) cat("complete:", cvg$complete, "\n")
