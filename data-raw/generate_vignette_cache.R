# Regenerate every stored result that the vignettes load in place of a slow
# fit. Each vignette shows the real call with eval = FALSE and loads the stored
# result in a hidden chunk. The code in each section must match its vignette:
# the same object names, arguments, and seeds. Keep the two in step.
#
# Run from the package root with the current development build installed. The
# vignettes are built against the installed package, and the Bayesian
# parallel runs start workers that load it.
#   Rscript data-raw/generate_vignette_cache.R                  # every section
#   Rscript data-raw/generate_vignette_cache.R quickstart nmix  # named sections
#
# Section               Vignette                   Output in inst/extdata/          Time
# quickstart            v01-quickstart             cache-quickstart.rds             ~10 s
# user-guide            v02-user-guide             cache-user-guide.rds             ~20 s
# nmix                  v02-user-guide             cache-user-guide-nmix.rds        ~5 min
# spatial-fields        v11-spatial-fields         cache-spatial-fields.rds         ~2 min
# distance-removal      v12-distance-removal       cache-distance-removal.rds       ~15 s
# bayesian-scales       v13-bayesian-scales        cache-bayesian-scales.rds        ~1 min
# bayesian-computation  v15-bayesian-computation   cache-bayesian-computation.rds   several min
#
# Every file is xz compressed. Check the built tarball size after a rebuild,
# because CRAN asks that packages stay under 5 MB.

library(multiScaleR)
library(terra)

# Stored results must come from the code in this checkout.
local({
  source_version <- read.dcf("DESCRIPTION", fields = "Version")[1, 1]
  installed <- as.character(utils::packageVersion("multiScaleR"))
  if (package_version(installed) != package_version(source_version))
    stop("Install the current build first: the installed multiScaleR is ",
         installed, " and DESCRIPTION says ", source_version, ".", call. = FALSE)
})

save_cache <- function(x, file) {
  path <- file.path("inst", "extdata", file)
  saveRDS(x, path, compress = "xz")
  message("Saved ", path, " (", format(file.size(path), big.mark = ","), " bytes)")
  invisible(path)
}

# Drop what summary(), coef(), and plot() do not need: the optimization
# context, the starting model, and the per-cell and binned landscape data.
strip_fit <- function(x) {
  x$opt_context           <- NULL
  x$fitted_mod_original   <- NULL
  x$kernel_inputs$binned  <- NULL
  x$kernel_inputs$d_list  <- NULL
  x$kernel_inputs$raw_cov <- NULL
  x
}

# Aggregated, wrapped rasters are small and survive saveRDS().
pack_plot_raster <- function(x, fact = 4) {
  terra::wrap(terra::aggregate(x, fact = fact, fun = mean, na.rm = TRUE))
}

# Each section is a quoted block evaluated in the global environment, as a
# top-level script would run. Model formulas then keep the global environment,
# which saveRDS() does not serialize. Inside a function they would capture the
# function frame and carry every raster and data object into the stored file.
sections <- list()

# ---- quickstart: the running example (v01-quickstart) ----------------------
# The Quick Start shows these calls with eval = FALSE. The object names match
# the vignette, so the stored call prints exactly as written there.
sections$quickstart <- quote({
  data("landscape_counts", package = "multiScaleR")
  data("surv_pts", package = "multiScaleR")
  pts <- vect(surv_pts)
  land_rast <- rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))
  # Lean inputs keep the stored fit small. The estimates match default inputs.
  kernel_inputs <- kernel_prep(pts = pts, raster_stack = land_rast,
                               max_D = 1700, kernel = "gaussian",
                               store_cell_data = FALSE, verbose = FALSE)
  df <- data.frame(landscape_counts, kernel_inputs$kernel_dat)
  mod <- glm(counts ~ site + land1 + land2, family = poisson(), data = df)
  invisible(utils::capture.output(
    opt <- multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs)))
  prof <- profile_sigma(opt, n_pts = 15, verbose = FALSE)
  r_scaled <- kernel_scale.raster(raster_stack = land_rast, multiScaleR = opt,
                                  scale_center = TRUE, clamp = TRUE,
                                  verbose = FALSE)
  pred <- terra::predict(r_scaled, opt$opt_mod, type = "response")
  save_cache(list(opt = opt, prof = prof,
                  r_scaled = pack_plot_raster(r_scaled),
                  pred = pack_plot_raster(pred),
                  built_with = as.character(packageVersion("multiScaleR"))),
             "cache-quickstart.rds")
})

# ---- user-guide: zero-inflated and occupancy fits (v02-user-guide) ---------
# Each takes several seconds to simulate, prepare, and fit.
sections[["user-guide"]] <- quote({
  library(unmarked)
  ## ---- must match the vignette: zero-inflated counts ----
  set.seed(12345)
  zi_rast <- sim_rast(user_seed = 999, dim = 500, resolution = 30)
  zi_rast <- subset(zi_rast, c("bin1", "cont1"))

  interior <- buffer(vect(ext(zi_rast)), -1000)
  zi_pts   <- spatSample(zi_rast[[1]], size = 400, method = "random",
                         ext = ext(interior), as.points = TRUE)

  truth <- kernel_prep(pts = zi_pts, raster_stack = zi_rast, max_D = 1500,
                       kernel = "gaussian", sigma = c(400, 200),
                       verbose = FALSE)$kernel_dat

  excess_zero <- rbinom(400, size = 1, prob = plogis(-0.25 - 1 * truth$cont1))
  abundance   <- rnbinom(400, mu = exp(0.5 + 0.75 * truth$bin1), size = 2)
  zi_dat      <- data.frame(cnt = ifelse(excess_zero == 1, 0, abundance))

  zi_inputs <- kernel_prep(pts = zi_pts, raster_stack = zi_rast, max_D = 1500,
                           kernel = "gaussian", verbose = FALSE)
  zi_df <- data.frame(zi_dat, zi_inputs$kernel_dat)

  zinb_mod <- pscl::zeroinfl(cnt ~ bin1 | cont1, dist = "negbin", data = zi_df)
  zinb_opt <- multiScale_optim(zinb_mod, zi_inputs, verbose = FALSE)
  ## ---- end zero-inflated counts ----

  ## ---- must match the vignette: occupancy ----
  um_rast <- subset(sim_rast(user_seed = 999, dim = 250), c("bin1", "cont2"))

  occ_sim <- sim_dat_unmarked(alpha = 0.25, beta = c(-0.75, 1),
                              sigma = c(225, 75),
                              kernel = "gaussian", type = "occ",
                              n_points = 250, n_surv = 5, det = 0.5,
                              raster_stack = um_rast, max_D = 800,
                              user_seed = 555)

  occ_inputs <- kernel_prep(pts = occ_sim$pts, raster_stack = um_rast,
                            max_D = 800, kernel = "gaussian", verbose = FALSE)
  occ_umf <- unmarkedFrameOccu(y = occ_sim$y, siteCovs = occ_inputs$kernel_dat)

  occ_mod <- occu(~ 1 ~ bin1 + cont2, data = occ_umf)
  opt_occ <- multiScale_optim(occ_mod, occ_inputs, verbose = FALSE)
  ## ---- end occupancy ----

  save_cache(list(zinb_opt = strip_fit(zinb_opt), opt_occ = strip_fit(opt_occ),
                  built_with = as.character(packageVersion("multiScaleR"))),
             "cache-user-guide.rds")
})

# ---- nmix: the N-mixture fit (v02-user-guide) -------------------------------
# The pcount optimization refits a K = 100 N-mixture model at every optimizer
# step, so it is the slowest fit in the Guide.
sections$nmix <- quote({
  library(unmarked)
  ## --- This block must match the vignette's pcount example exactly ---
  rs <- sim_rast(user_seed = 999, dim = 250)
  rs <- terra::subset(rs, c(1, 4))
  s_dat <- sim_dat_unmarked(alpha = 1, beta = c(0.75, -0.75), kernel = "gaussian",
                            sigma = c(75, 150), n_points = 75, n_surv = 5, det = 0.5,
                            type = "count", raster_stack = rs, max_D = 550, user_seed = 555)
  kernel_inputs <- kernel_prep(pts = s_dat$pts, raster_stack = rs, max_D = 550,
                               kernel = "gaus", verbose = FALSE)
  umf <- unmarkedFramePCount(y = s_dat$y, siteCovs = kernel_inputs$kernel_dat)
  mod0_umf.p <- pcount(~1 ~bin1 + cont2, data = umf, K = 100)

  opt_umf.p <- multiScale_optim(fitted_mod = mod0_umf.p, kernel_inputs = kernel_inputs)

  ## --- Strip to the elements needed for summary()/plot()/projection ---
  ## (drops the per-cell data, binned summaries, and the optimization context /
  ## original model, which are only needed to re-optimize or profile.)
  opt_umf.p$opt_context         <- NULL
  opt_umf.p$fitted_mod_original <- NULL
  opt_umf.p$kernel_inputs$binned  <- NULL
  opt_umf.p$kernel_inputs$d_list  <- NULL
  opt_umf.p$kernel_inputs$raw_cov <- NULL

  save_cache(opt_umf.p, "cache-user-guide-nmix.rds")
})

# ---- spatial-fields: glmmTMB response field (v11-spatial-fields) ------------
# Each candidate scale needs a glmmTMB refit. The stored fit drops its copy of
# the kernel inputs, which the vignette attaches from its own identical call.
sections[["spatial-fields"]] <- quote({
  stopifnot(requireNamespace("glmmTMB", quietly = TRUE))
  data("landscape_counts")
  data("surv_pts")
  pts <- vect(surv_pts)
  xy  <- crds(pts)
  land_rast <- rast(system.file("extdata", "landscape.tif",
                                package = "multiScaleR"))

  # ---- Simulated counts: same code as the vignette chunk `simulate` ----------
  truth  <- kernel_scale.raster(land_rast[[c("land1", "land2")]],
                                sigma = c(250, 500), kernel = "gaussian",
                                verbose = FALSE)
  x_true <- scale(as.matrix(extract(truth, pts, ID = FALSE)))

  set.seed(1)
  d_mat <- as.matrix(dist(xy))
  field <- drop(t(chol(0.6^2 * exp(-d_mat / 800))) %*% rnorm(nrow(xy)))
  site  <- rnorm(nrow(xy))
  log_mu <- 1 + 0.3 * site - 0.5 * x_true[, "land1"] +
    0.7 * x_true[, "land2"] + field
  dat <- data.frame(counts = rpois(nrow(xy), exp(log_mu)), site = site)

  # ---- Same preparation and starting model as the vignette -------------------
  inputs <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
                        kernel = "gaussian", bin = TRUE, verbose = FALSE)
  df <- data.frame(dat, inputs$kernel_dat)
  df$pos <- glmmTMB::numFactor(xy[, 1], xy[, 2])
  df$grp <- factor(rep("all", nrow(df)))

  start_field <- glmmTMB::glmmTMB(
    counts ~ site + land1 + land2 + exp(pos + 0 | grp),
    family = poisson(), data = df,
    start  = list(theta = c(log(0.5), log(500)))
  )

  # Same call as the vignette chunk `fit-field`
  t0 <- Sys.time()
  fit_field <- multiScale_optim(start_field, inputs, n_cores = 1,
                                verbose = FALSE)
  fit_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Drop the large duplicated pieces. The vignette restores kernel_inputs.
  fit_field$kernel_inputs       <- NULL
  fit_field$opt_context         <- NULL
  fit_field$fitted_mod_original <- NULL

  # ---- The bundled counts, which were simulated without a spatial field ------
  # The vignette reports these numbers in one sentence to show what a field
  # estimates when the data hold no residual spatial pattern.
  dfr <- data.frame(landscape_counts, inputs$kernel_dat)
  dfr$pos <- df$pos
  dfr$grp <- df$grp
  run_glm <- multiScale_optim(
    glm(counts ~ site + land1 + land2, family = poisson(), data = dfr),
    inputs, n_cores = 1, verbose = FALSE)
  run_field <- multiScale_optim(
    glmmTMB::glmmTMB(counts ~ site + land1 + land2 + exp(pos + 0 | grp),
                     family = poisson(), data = dfr,
                     start = list(theta = c(log(0.5), log(500)))),
    inputs, n_cores = 1, verbose = FALSE)
  bundled <- list(
    sigma_glm   = run_glm$scale_est,
    sigma_field = run_field$scale_est,
    field_par   = exp(glmmTMB::getME(run_field$opt_mod, "theta"))  # SD, range (m)
  )

  cache <- list(
    fit_field   = fit_field,
    dat         = dat,
    bundled     = bundled,
    fit_seconds = fit_seconds,
    versions    = list(multiScaleR = as.character(packageVersion("multiScaleR")),
                       glmmTMB = as.character(packageVersion("glmmTMB")),
                       TMB = as.character(packageVersion("TMB"))),
    created     = Sys.time()
  )

  save_cache(cache, "cache-spatial-fields.rds")
})

# ---- distance-removal: GDR abundance field table (v12-distance-removal) ----
# fit_gdr_response_field() runs only with unmarked 1.5.2 and TMB 1.9.25, so the
# vignette shows these calls with eval = FALSE and prints the stored table.
sections[["distance-removal"]] <- quote({
  stopifnot(as.character(utils::packageVersion("unmarked")) == "1.5.2",
            as.character(utils::packageVersion("TMB")) == "1.9.25")

  # ---- Landscape and survey points (same inputs as the vignette) -------------
  data("surv_pts")
  pts <- vect(surv_pts)
  land_rast <- rast(system.file("extdata", "landscape.tif",
                                package = "multiScaleR"))

  # Kernel-weighted habitat at a 250 m scale, standardized, one row per point
  inputs_250 <- kernel_prep(pts = pts, raster_stack = land_rast[["land1"]],
                            max_D = 1700, kernel = "gaussian", sigma = 250,
                            verbose = FALSE)

  # ---- Sixteen square blocks over the extent of the survey points ------------
  xy <- crds(pts)
  n_side <- 4
  width <- diff(range(xy[, 1])) / n_side
  height <- diff(range(xy[, 2])) / n_side
  col_id <- pmin(n_side - 1, floor((xy[, 1] - min(xy[, 1])) / width))
  row_id <- pmin(n_side - 1, floor((xy[, 2] - min(xy[, 2])) / height))
  block <- factor(paste0("b", col_id, row_id))

  # Every survey in a block carries the coordinates of the block center
  block_xy <- cbind(x = min(xy[, 1]) + (col_id + 0.5) * width,
                    y = min(xy[, 2]) + (row_id + 0.5) * height)

  # ---- Block effects with an exponential correlation range of 1500 m ---------
  field_range <- 1500
  field_sd <- 0.5
  centers <- block_xy[!duplicated(block), , drop = FALSE]
  corr <- exp(-as.matrix(dist(centers)) / field_range)
  set.seed(3)
  effect <- drop(t(chol(corr)) %*% rnorm(nrow(centers))) * field_sd
  block_effect <- effect[match(block, block[!duplicated(block)])]

  # ---- Simulated surveys: density depends on habitat and the block effect ----
  sim_blocks <- sim_dat_gdr(
    density          = exp(log(400) + 0.5 * inputs_250$kernel_dat$land1 +
                             block_effect),
    distance_scale   = 40,
    removal          = 0.2,
    dist_breaks      = c(0, 25, 50, 75),
    interval_lengths = c(2, 3, 5),
    interval_units   = "minutes",
    site_covs        = data.frame(inputs_250$kernel_dat, block = block),
    seed             = 102
  )
  umf_blocks <- sim_blocks$observed
  rownames(block_xy) <- rownames(unmarked::siteCovs(umf_blocks))

  # ---- The calls shown in the vignette ----------------------------------------
  # 1. Template: independent random intercept for each block (TMB engine)
  template <- unmarked::gdistremoval(
    lambdaformula   = ~ land1 + (1 | block),
    removalformula  = ~ 1,
    distanceformula = ~ 1,
    data = umf_blocks, keyfun = "halfnorm", mixture = "P",
    output = "density", unitsOut = "kmsq", K = 80, engine = "TMB",
    starts = c(log(400), 0, log(40), qlogis(0.2), log(0.5))
  )

  # 2. A second, independent set of starting values, named by parameter block
  cold <- c("beta_lambda[1]"   = log(300),    # log density at average habitat
            "beta_lambda[2]"   = 0,           # habitat effect
            "beta_dist[1]"     = log(30),     # log distance scale (m)
            "beta_rem[1]"      = qlogis(0.4), # logit removal probability
            "lsigma_lambda[1]" = log(0.3))    # log SD of the block effects
  settings <- list(method = "BFGS", se = TRUE, threads = 1L,
                   control = list(maxit = 2500L, reltol = 1e-10))

  # 3. Replace the independent block effects with a spatially correlated field
  elapsed <- system.time(
    field_fit <- fit_gdr_response_field(
      model = template, coords = block_xy, crs = 32617, range = field_range,
      cold_starts = list(cold = cold),
      audit_dir = tempfile("gdr-field-audit-"),
      tmb_settings = settings
    )
  )[["elapsed"]]

  # ---- Small tables for the vignette ------------------------------------------
  block_sd <- function(fit) {
    # SD of the block effects, on the log density scale
    as.numeric(unmarked::sigma(fit)$sigma[1])
  }
  row_for <- function(fit, label) {
    est <- unmarked::coef(fit)
    data.frame(
      model        = label,
      density      = unname(exp(est["lambda(Int)"])),
      land1_effect = unname(est["lambda(land1)"]),
      block_sd     = block_sd(fit),
      logLik       = as.numeric(unmarked::logLik(fit)),
      AIC          = fit@AIC
    )
  }
  compare <- rbind(row_for(template, "Independent blocks"),
                   row_for(field_fit, "Correlated field"))

  cache <- list(
    field = list(
      compare     = compare,
      n_blocks    = nlevels(block),
      range       = field_range,
      class       = class(field_fit)[1],
      seconds     = round(elapsed),
      versions    = c(unmarked = as.character(utils::packageVersion("unmarked")),
                      TMB = as.character(utils::packageVersion("TMB")),
                      multiScaleR = as.character(utils::packageVersion("multiScaleR")))
    )
  )

  save_cache(cache, "cache-distance-removal.rds")
  print(cache$field$compare)
})

# ---- bayesian-scales: one slice-sampling fit (v13-bayesian-scales) ----------
# The stored fit keeps its ordered chains, so diagnose_bayes_scales() and
# chain plots work on it. Needs the 'posterior' package.
sections[["bayesian-scales"]] <- quote({
  stopifnot(requireNamespace("posterior", quietly = TRUE))
  data("landscape_counts")
  data("surv_pts")
  land_rast <- terra::rast(system.file("extdata", "landscape.tif",
                                       package = "multiScaleR"))

  # Same preparation as the vignette: binned summaries within 1700 m
  inputs <- kernel_prep(pts = terra::vect(surv_pts), raster_stack = land_rast,
                        max_D = 1700, kernel = "gaussian", bin = TRUE,
                        verbose = FALSE)

  # Same call as the vignette chunk `fit-shared`
  t0 <- Sys.time()
  fit <- multiScale_bayes(counts ~ site + land1 + land2, inputs,
                          landscape_counts, family = "poisson",
                          priors = list(scale_median = c(land1 = 300,
                                                         land2 = 300)),
                          method = "slice", seed = 5211)
  seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Chain diagnostics for the fitted log distances (needs 'posterior')
  precision <- diagnose_bayes_scales(fit)
  stopifnot(isTRUE(precision$available))

  # The same model with the reference scales moved to 250 m and 500 m. The
  # vignette quotes its coefficients to show that a Bayesian coefficient refers
  # to one standard deviation of the covariate at scale_median. Only the two
  # summary tables are stored.
  fit_ref <- multiScale_bayes(counts ~ site + land1 + land2, inputs,
                              landscape_counts, family = "poisson",
                              priors = list(scale_median = c(land1 = 250,
                                                             land2 = 500)),
                              method = "slice", seed = 5211)
  reference_check <- list(coefficients = summary(fit_ref)$coefficients,
                          scales = posterior_scales(fit_ref),
                          passed = summary(fit_ref)$diagnostics$passed)

  cache <- list(fit = fit, precision = precision, seconds = seconds,
                reference_check = reference_check,
                package_version = as.character(packageVersion("multiScaleR")),
                created = Sys.Date())

  save_cache(cache, "cache-bayesian-scales.rds")
})

# ---- bayesian-computation: run summaries (v15-bayesian-computation) ---------
# A short pilot and one long slice-sampling run of the bundled five-year
# example. Only their summaries are stored. Parallel chains need 'processx'.
sections[["bayesian-computation"]] <- quote({
  example <- readRDS(file.path("inst", "extdata", "bayesian-annual-example.rds"))
  dat      <- example$data
  inputs   <- example$inputs
  response <- count ~ habitat + factor(year) + (1 | site)
  scales   <- list(habitat = ~1 + (1 | year))

  priors   <- list(scale_median = 250, scale_group_sd = 0.5)

  # Keep the summaries of one slice fit that the vignette prints or cites
  keep_summaries <- function(fit, seconds) {
    list(
      seconds     = seconds,
      sampling    = fit$sampling,                       # chains, iter, warmup, seed, cores
      fit_bytes   = as.numeric(object.size(fit)),
      chain_dim   = dim(fit$chains),                    # iterations, chains, parameters
      diagnostics = fit$diagnostics,                    # stored sampler checks
      precision   = diagnose_bayes_scales(fit)[c("available", "parameters",
                      "quantiles", "chains", "iterations", "derived_passed",
                      "passed", "criteria")],
      scales      = posterior_scales(fit),              # annual scales, meters
      group_sd    = summary(fit)$group_sd
    )
  }

  # 1. Short pilot run on one worker (object name in the vignette: pilot)
  pilot_seconds <- system.time(
    pilot <- multiScale_bayes(
      response, inputs, dat, scales = scales, row_id = "site_year",
      family = "nbinom2", priors = priors,
      method = "slice", chains = 2, iter = 1000, seed = 5212
    )
  )[["elapsed"]]

  # 2. The same pilot on two workers, to record the cost of parallel startup
  pilot2_seconds <- system.time(
    pilot2 <- multiScale_bayes(
      response, inputs, dat, scales = scales, row_id = "site_year",
      family = "nbinom2", priors = priors,
      method = "slice", chains = 2, iter = 1000, seed = 5212, cores = 2
    )
  )[["elapsed"]]

  # 3. Long run on four workers (object name in the vignette: slice_long)
  long_seconds <- system.time(
    slice_long <- multiScale_bayes(
      response, inputs, dat, scales = scales, row_id = "site_year",
      family = "nbinom2", priors = priors,
      method = "slice", chains = 4, iter = 8000, seed = 5213, cores = 4
    )
  )[["elapsed"]]

  cache <- list(
    pilot = keep_summaries(pilot, pilot_seconds),
    pilot_two_cores = list(seconds = pilot2_seconds,
                           same_draws = identical(pilot$chains, pilot2$chains)),
    long  = keep_summaries(slice_long, long_seconds),
    # Precision of P(sigma <= threshold) for each year, from the long run
    probabilities = diagnose_bayes_probabilities(
      slice_long, thresholds = list(habitat = c(100, 250)), mcse_target = 0.005
    )[c("available", "probabilities", "mcse_target")],
    session = list(multiScaleR = as.character(packageVersion("multiScaleR")),
                   R = R.version.string, date = Sys.Date())
  )

  save_cache(cache, "cache-bayesian-computation.rds")
})

# ---- run the requested sections ---------------------------------------------
requested <- commandArgs(trailingOnly = TRUE)
if (!length(requested)) requested <- names(sections)
unknown <- setdiff(requested, names(sections))
if (length(unknown))
  stop("Unknown section(s): ", paste(unknown, collapse = ", "),
       ". Choose from: ", paste(names(sections), collapse = ", "), call. = FALSE)
keep <- c(ls(globalenv(), all.names = TRUE), "keep", "name", "started")
for (name in requested) {
  message("== ", name)
  started <- Sys.time()
  eval(sections[[name]], envir = globalenv())
  message("   ", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1),
          " min")
  # Clear the section's objects before the next one runs.
  rm(list = setdiff(ls(globalenv(), all.names = TRUE), keep), envir = globalenv())
}
