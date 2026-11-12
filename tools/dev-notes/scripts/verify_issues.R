# Independent reproduction of code issues reported during the vignette rollout.
# Run from the worktree root. Uses the private scratch build only.
scr <- "C:/Users/peterman.73/AppData/Local/Temp/claude/9e3b7d3898bd820b/459cb932-8333-43df-b0da-7d9a2422aa59/scratchpad"
.libPaths(c(file.path(scr, "rlib"), .libPaths()))
suppressPackageStartupMessages({library(multiScaleR); library(terra); library(splines); library(survival)})
grDevices::pdf(NULL)
chk <- function(label, expr) {
  t0 <- Sys.time(); warns <- character()
  val <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")}),
    error = function(e) structure(conditionMessage(e), class = "err"))
  cat(sprintf("\n### %s  [%.1fs]\n", label, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (inherits(val, "err")) cat("  ERROR:", substr(unclass(val), 1, 300), "\n")
  if (length(warns)) cat("  WARNINGS:", paste(substr(unique(warns), 1, 160), collapse = " | "), "\n")
  invisible(val)
}
ok <- function(x) !inherits(x, "err")

data("landscape_counts"); data("surv_pts")
land <- rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))
inputs <- kernel_prep(pts = surv_pts, raster_stack = land, max_D = 1700, kernel = "gaussian", bin = TRUE, verbose = FALSE)
dat <- data.frame(landscape_counts, inputs$kernel_dat)
f0 <- counts ~ site + land1 + land2
mod <- glm(f0, family = poisson(), data = dat)

fit_refit <- chk("A. default refit fit: is there a `backend` element?", multiScale_optim(mod, inputs, n_cores = 1, verbose = FALSE))
cat("  'backend' %in% names(fit):", "backend" %in% names(fit_refit), " value:", deparse(fit_refit$backend), "\n")

in_fix <- kernel_prep(pts = surv_pts, raster_stack = land, max_D = 1700, kernel = "fixed", verbose = FALSE)
d_fix <- data.frame(landscape_counts, in_fix$kernel_dat); m_fix <- glm(f0, family = poisson(), data = d_fix)
fx1 <- chk("B1. fixed kernel, default start", multiScale_optim(m_fix, in_fix, n_cores = 1, verbose = FALSE))
if (ok(fx1)) { print(fx1$scale_est); cat("  logLik:", as.numeric(logLik(fx1$opt_mod)), "\n") }
fx2 <- chk("B2. fixed kernel, start_strategy = 'screen'", multiScale_optim(m_fix, in_fix, n_cores = 1, verbose = FALSE, start_strategy = "screen"))
if (ok(fx2)) { print(fx2$scale_est); cat("  logLik:", as.numeric(logLik(fx2$opt_mod)), "\n") }

fit_pg <- chk("C0. profile_glm fit", multiScale_optim(mod, inputs, backend = "profile_glm", n_cores = 1, verbose = FALSE))
p1 <- chk("C1. profile_sigma(profile_glm fit) with defaults", profile_sigma(fit_pg, verbose = FALSE))
p2 <- chk("C2. profile_sigma(profile_glm fit, sigma_range = c(100, 600))", profile_sigma(fit_pg, sigma_range = c(100, 600), verbose = FALSE))
p3 <- chk("C3. profile_sigma(default refit fit) with defaults", profile_sigma(fit_refit, verbose = FALSE))

dat$logc <- log1p(dat$counts)
m_lm <- lm(logc ~ site + land1 + land2, data = dat)
m_gg <- glm(logc ~ site + land1 + land2, family = gaussian(), data = dat)
chk("D1. profile_glm with lm()", multiScale_optim(m_lm, inputs, backend = "profile_glm", n_cores = 1, verbose = FALSE))
chk("D2. profile_glm with glm(gaussian)", multiScale_optim(m_gg, inputs, backend = "profile_glm", n_cores = 1, verbose = FALSE))

xy <- as.data.frame(terra::crds(terra::vect(surv_pts)))
folds <- spatial_cv_folds(xy, block_size = 1000, k = 3, crs = 32617)
cv_p  <- chk("E0a. Poisson CV, default metrics", cross_validate(f0, dat, inputs, folds, family = "poisson", keep_selection = TRUE))
cv_nb <- chk("E0b. NB2 CV, default metrics", cross_validate(f0, dat, inputs, folds, family = "nbinom2", keep_selection = TRUE))
if (ok(cv_p) && ok(cv_nb)) {
  cat("  Poisson metrics:", cv_p$settings$metrics, " NB2 metrics:", cv_nb$settings$metrics, "\n")
  chk("E1. select_cv(mixed Poisson + NB2, defaults)", select_cv(list(pois = cv_p, nb = cv_nb), criterion = "AICc"))
  cv_p2 <- cross_validate(f0, dat, inputs, folds, family = "poisson", keep_selection = TRUE, metrics = c("rmse", "mae", "neg_log_score"))
  chk("E2. select_cv after matching Poisson metrics", select_cv(list(pois = cv_p2, nb = cv_nb), criterion = "AICc"))
}
lean <- kernel_prep(pts = surv_pts, raster_stack = land, max_D = 1700, kernel = "gaussian", verbose = FALSE)
cat("\n  default kernel_prep store_cell_data formal:", deparse(formals(kernel_prep)$store_cell_data), "\n")
cv_l <- chk("E3. cross_validate with default kernel_prep() inputs", cross_validate(f0, data.frame(landscape_counts, lean$kernel_dat), lean, folds, family = "poisson"))
if (ok(cv_l)) cat("  complete:", cv_l$complete, "\n")
lean2 <- chk("E4. kernel_prep(store_cell_data = FALSE)", kernel_prep(pts = surv_pts, raster_stack = land, max_D = 1700, kernel = "gaussian", store_cell_data = FALSE, verbose = FALSE))
if (ok(lean2)) { cv_l2 <- chk("E5. cross_validate with store_cell_data = FALSE inputs", cross_validate(f0, data.frame(landscape_counts, lean2$kernel_dat), lean2, folds, family = "poisson")); if (ok(cv_l2)) cat("  complete:", cv_l2$complete, "\n") }

cat("\n### F. msr_vars: base on msidi; rpr without classes_max\n")
cover <- classify(land$land3, rcl = matrix(c(-Inf, .45, 1, .45, .55, 2, .55, Inf, 3), ncol = 3, byrow = TRUE)); names(cover) <- "cover"
mk <- function(...) kernel_prep(pts = surv_pts, raster_stack = cover, max_D = 800, scale_vars = msr_vars(...), verbose = FALSE)
a <- chk("F1. msidi default base", mk(m = landscape_var("cover", metric = "msidi", radius = 500)))
b <- chk("F2. msidi base = 2", mk(m = landscape_var("cover", metric = "msidi", radius = 500, base = 2)))
if (ok(a) && ok(b)) cat("  identical raw values:", isTRUE(all.equal(attr(a$kernel_dat, "scaled:center"), attr(b$kernel_dat, "scaled:center"))), " centers:", unlist(a$scl_params$mean), unlist(b$scl_params$mean), "\n")
s1 <- chk("F3. shdi default base", mk(m = landscape_var("cover", metric = "shdi", radius = 500)))
s2 <- chk("F4. shdi base = 2", mk(m = landscape_var("cover", metric = "shdi", radius = 500, base = 2)))
if (ok(s1) && ok(s2)) cat("  shdi centers:", unlist(s1$scl_params$mean), unlist(s2$scl_params$mean), "\n")
chk("F5. rpr without classes_max", mk(m = landscape_var("cover", metric = "rpr", radius = 500)))
chk("F6. rpr with classes_max = 3", mk(m = landscape_var("cover", metric = "rpr", radius = 500, classes_max = 3)))

cat("\n### G. plot_kernel dashed-line colour by kernel\n")
for (k in c("gaussian", "exp", "expow", "fixed")) {
  p <- tryCatch(if (k == "expow") plot_kernel(sigma = 300, kernel = k, beta = 5) else plot_kernel(sigma = 300, kernel = k), error = function(e) conditionMessage(e))
  if (inherits(p, "ggplot")) { gb <- ggplot2::ggplot_build(p); cols <- unlist(lapply(gb$data, function(d) unique(d$colour))); cat(" ", k, ": layer colours =", paste(cols, collapse = ", "), "\n") } else cat(" ", k, ":", p, "\n")
}

cat("\n### H. user-name paths in shipped rds files\n")
for (f in list.files(system.file("extdata", package = "multiScaleR"), pattern = "[.]rds$", full.names = TRUE)) {
  txt <- utils::capture.output(str(readRDS(f), max.level = 6, vec.len = 4, nchar.max = 300, list.len = 200))
  hits <- grep("Users/|Users\\\\", txt, value = TRUE)
  cat(" ", basename(f), ": lines with a Users path =", length(hits), "\n"); if (length(hits)) cat("    e.g.", gsub("Users[/\\\\]+[^/\\\\ \"]+", "Users/<name>", substr(trimws(hits[1]), 1, 150)), "\n")
}
cat("\n### I. land2 value range:", paste(round(unlist(terra::global(land$land2, "range", na.rm = TRUE)), 3), collapse = " to "), "\n")

cat("\n### J. conditional_logit\n")
setwd("vignettes")
src <- readLines(file.path(scr, "purl_step-selection.R")); eval(parse(text = src[1:90]), envir = globalenv())
setwd("..")
cat("  sigma SE:", opt$scale_est["forest", "SE"], "\n")
chk("J1. profile_sigma(conditional_logit fit)", profile_sigma(opt, verbose = FALSE))
chk("J2. summary(conditional_logit fit, profile = TRUE)", summary(opt, profile = TRUE))
m_sp <- clogit(case_ ~ ns(forest, knots = 0, Boundary.knots = c(-3, 3)) + sl_km + strata(step_id), data = df, model = TRUE)
chk("J3. conditional_logit with fixed-knot ns()", multiScale_optim(fitted_mod = m_sp, kernel_inputs = inputs, backend = "conditional_logit", n_cores = 1, verbose = FALSE))
m_sp2 <- clogit(case_ ~ splines::bs(forest, knots = 0, Boundary.knots = c(-4, 4)) + sl_km + strata(step_id), data = df, model = TRUE)
chk("J4. conditional_logit with fixed-knot splines::bs()", multiScale_optim(fitted_mod = m_sp2, kernel_inputs = inputs, backend = "conditional_logit", n_cores = 1, verbose = FALSE))

cat("\n### K. glmmTMB spatial field\n")
if (suppressWarnings(requireNamespace("glmmTMB", quietly = TRUE))) {
  setwd("vignettes")
  src <- readLines(file.path(scr, "purl_spatial-response-fields.R")); eval(parse(text = src[1:76]), envir = globalenv())
  setwd("..")
  cat("  'backend' in cached field fit:", "backend" %in% names(fit_field), "\n")
  chk("K1. response_field_geometry(real field fit)", response_field_geometry(fit_field, coords = xy, crs = 32617))
  ff <- fit_field; ff$backend <- "refit"
  chk("K2. same, after setting fit$backend <- 'refit' by hand", response_field_geometry(ff, coords = xy, crs = 32617))
  g_xy <- cbind(x = c(1010, 2030), y = c(1500, 2500))
  g_in <- kernel_prep(pts = vect(g_xy, type = "points"), raster_stack = land_rast, max_D = 1700, kernel = "gaussian", bin = TRUE, verbose = FALSE)
  g_dat <- data.frame(site = 0, pos = glmmTMB::numFactor(g_xy[, 1], g_xy[, 2]), grp = factor("all"))
  chk("K3. predict at new positions, no allow.new.levels", predict(fit_field, g_dat, g_in))
  chk("K4. predict at new positions, allow.new.levels = TRUE", predict(fit_field, g_dat, g_in, allow.new.levels = TRUE))
  s0 <- chk("K5. glmmTMB field starting model with DEFAULT starts", suppressMessages(glmmTMB::glmmTMB(counts ~ site + land1 + land2 + exp(pos + 0 | grp), family = poisson(), data = df)))
  if (ok(s0)) {
    cat("  starting model pdHess:", s0$sdr$pdHess, " theta:", round(glmmTMB::getME(s0, "theta"), 3), "\n")
    f0f <- chk("K6. multiScale_optim on that default-start model", multiScale_optim(s0, inputs, n_cores = 1, verbose = FALSE))
    if (ok(f0f)) print(f0f$scale_est)
  }
} else cat("  glmmTMB not installed\n")
cat("\nDONE\n")
