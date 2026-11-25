# Integration matrix: every post-fit function on one small real fit of each
# supported model class and backend.
#
# Each cell has a declared outcome. A cell is "ok" unless the model's `refuse`
# list gives a regular expression for a deliberate, documented refusal. An
# undeclared error fails the test, so a new model class, backend, or post-fit
# function must state how the two work together. The matrix also checks that
# a scale reported without a standard error is explained in the printed
# summary.
#
# Set MSR_MATRIX_DISCOVER=true to print every cell's outcome instead of
# asserting (used when adding a row or a column).
skip_on_cran()

if (requireNamespace("survival", quietly = TRUE)) {
  strata <- survival::strata; coxph <- survival::coxph
}

discover <- identical(Sys.getenv("MSR_MATRIX_DISCOVER"), "true")

matrix_draw <- function(x) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  if (inherits(x, "ggplot")) print(x) else
    for (p in x) if (inherits(p, "ggplot")) print(p)
  invisible(x)
}

# Shared small inputs: a projected 60 x 60 landscape and 60 survey points.
matrix_inputs <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    rs <- terra::subset(sim_rast(dim = 60, resolution = 10, user_seed = 555), c(1, 4))
    terra::crs(rs) <- "EPSG:32617"
    sim <- sim_dat(alpha = 0.25, beta = c(0.75, -0.5), kernel = "gaussian",
                   sigma = c(60, 100), type = "count", n_points = 60,
                   raster_stack = rs, max_D = 250, user_seed = 555)
    pts <- sim$pts
    sf::st_crs(pts) <- 32617
    prep <- function(kernel) kernel_prep(pts, rs, max_D = 250, kernel = kernel,
                                         verbose = FALSE)
    inputs <- prep("gaussian")
    dat <- data.frame(y = sim$obs, inputs$kernel_dat,
                      grp = factor(rep(seq_len(6), length.out = 60)))
    dat$yg <- as.numeric(scale(sim$obs))
    occ <- sim_dat_unmarked(alpha = 0, beta = c(1, -0.5), kernel = "gaussian",
                            sigma = c(60, 100), type = "occ", n_points = 60,
                            n_surv = 4, det = 0.5, raster_stack = rs, max_D = 250,
                            user_seed = 3)
    occ_pts <- occ$pts
    sf::st_crs(occ_pts) <- 32617
    occ_inputs <- kernel_prep(occ_pts, rs, max_D = 250, verbose = FALSE)
    cache <<- list(rs = rs, pts = pts, prep = prep, inputs = inputs, dat = dat,
                   occ = occ, occ_inputs = occ_inputs)
    cache
  }
})

# Choice sets for the step-selection rows: 40 sets of 5 steps.
matrix_steps_data <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    m <- matrix_inputs()
    rs <- m$rs
    set.seed(11)
    n_sets <- 40L; n_alt <- 5L
    sx <- rep(stats::runif(n_sets, 200, 400), each = n_alt)
    sy <- rep(stats::runif(n_sets, 200, 400), each = n_alt)
    len <- stats::rgamma(n_sets * n_alt, 2, scale = 25)
    ang <- stats::runif(n_sets * n_alt, -pi, pi)
    d <- data.frame(step_id = rep(seq_len(n_sets), each = n_alt),
                    x = sx + len * cos(ang), y = sy + len * sin(ang))
    v <- terra::extract(rs[[1L]], as.matrix(d[, c("x", "y")]))[, 1]
    score <- exp(as.numeric(scale(v)))
    d$case_ <- 0
    for (s in seq_len(n_sets)) {
      rows <- which(d$step_id == s)
      d$case_[sample(rows, 1, prob = score[rows])] <- 1
    }
    pts <- terra::vect(d, geom = c("x", "y"), crs = terra::crs(rs))
    inputs <- kernel_prep(pts, rs, max_D = 250, verbose = FALSE)
    cache <<- list(inputs = inputs, dat = data.frame(d, inputs$kernel_dat))
    cache
  }
})

# Post-fit functions (columns). `spec` supplies per-model arguments.
matrix_columns <- list(
  print = function(fit, spec) utils::capture.output(print(fit)),
  summary = function(fit, spec) utils::capture.output(print(summary(fit))),
  plot = function(fit, spec) matrix_draw(plot(fit)),
  marginal_effects = function(fit, spec) {
    out <- if (is.null(spec$type)) plot_marginal_effects(fit) else
      plot_marginal_effects(fit, type = spec$type)
    if (!length(out)) stop("no effect curves were returned")
    matrix_draw(out)
  },
  kernel_dist = function(fit, spec) kernel_dist(fit),
  predict = function(fit, spec) do.call(stats::predict, c(list(fit), spec$predict_args)),
  predict_conditional_se = function(fit, spec) do.call(stats::predict,
    c(list(fit, se.fit = TRUE, uncertainty = "conditional"), spec$predict_args)),
  predict_joint_se = function(fit, spec) do.call(stats::predict,
    c(list(fit, se.fit = TRUE, uncertainty = "joint"), spec$predict_args)),
  vcov = function(fit, spec) stats::vcov(fit),
  raster = function(fit, spec) kernel_scale.raster(spec$rs, multiScaleR = fit,
    verbose = FALSE, site_values = spec$site_values),
  raster_clamp = function(fit, spec) kernel_scale.raster(spec$rs, multiScaleR = fit,
    scale_center = TRUE, clamp = TRUE, verbose = FALSE, site_values = spec$site_values),
  raster_predict = function(fit, spec) {
    r <- kernel_scale.raster(spec$rs, multiScaleR = fit, scale_center = TRUE,
      match_buffer = TRUE, verbose = FALSE, site_values = spec$site_values)
    do.call(stats::predict, c(list(fit, newdata = r), spec$raster_args))
  },
  aic_tab = function(fit, spec) aic_tab(list(a = fit)),
  # Standard generics must describe the final model, never return NULL, and
  # never return the starting model (fitted_mod_original).
  coef = function(fit, spec) {
    value <- stats::coef(fit)
    if (is.null(value) || !length(value)) stop("coef() returned nothing")
    value
  },
  fitted = function(fit, spec) {
    value <- stats::fitted(fit)
    if (is.null(value) || is.list(value) && !is.data.frame(value))
      stop("fitted() did not return fitted values")
    value
  },
  logLik = function(fit, spec) {
    ll <- stats::logLik(fit)
    if (!is.finite(as.numeric(ll))) stop("logLik() is not finite")
    # AIC() must agree with the selection table, which counts the scales.
    table <- aic_tab(list(a = fit), AICc = FALSE)
    if (!isTRUE(all.equal(stats::AIC(fit), table$AIC[1L])))
      stop("AIC() does not match aic_tab()")
    ll
  },
  profile_sigma = function(fit, spec) matrix_draw(plot(profile_sigma(fit, n_pts = 3,
    verbose = FALSE))),
  scale_inference = function(fit, spec) scale_inference(fit, grid_points = 5L)
)

# Messages for refusals shared by every refit-backend model.
refit_refusals <- list(
  predict_joint_se = "conditional SEs only",
  vcov = "requires backend = 'profile_glm' or 'joint_nb2'",
  scale_inference = "scale_inference"
)

# Rows: one small real fit per model class and backend.
matrix_models <- list(
  list(name = "glm poisson, refit, gaussian kernel",
       fit = function(m) multiScale_optim(glm(y ~ bin1 + cont2, family = poisson,
         data = m$dat), m$inputs, verbose = FALSE),
       refuse = refit_refusals),
  list(name = "glm poisson, refit, expow kernel",
       fit = function(m) { i <- m$prep("expow")
         multiScale_optim(glm(y ~ bin1 + cont2, family = poisson,
           data = data.frame(y = m$dat$y, i$kernel_dat)), i, verbose = FALSE) },
       refuse = refit_refusals),
  list(name = "glm poisson, refit, fixed kernel",
       fit = function(m) { i <- m$prep("fixed")
         suppressWarnings(multiScale_optim(glm(y ~ bin1 + cont2, family = poisson,
           data = data.frame(y = m$dat$y, i$kernel_dat)), i, verbose = FALSE)) },
       refuse = refit_refusals),
  list(name = "glm poisson, profile_glm",
       fit = function(m) multiScale_optim(glm(y ~ bin1 + cont2, family = poisson,
         data = m$dat), m$inputs, backend = "profile_glm", verbose = FALSE),
       refuse = list(scale_inference = "NB2 with one spatial covariate")),
  list(name = "glm.nb, joint_nb2", packages = "MASS",
       fit = function(m) multiScale_optim(MASS::glm.nb(y ~ cont2, data = m$dat),
         m$inputs, backend = "joint_nb2", verbose = FALSE),
       refuse = list()),
  list(name = "lm, refit",
       fit = function(m) multiScale_optim(lm(yg ~ bin1 + cont2, data = m$dat),
         m$inputs, verbose = FALSE),
       refuse = refit_refusals),
  list(name = "gam, refit", packages = "mgcv",
       fit = function(m) multiScale_optim(mgcv::gam(y ~ bin1 + cont2, family = poisson,
         data = m$dat), m$inputs, verbose = FALSE),
       refuse = refit_refusals),
  list(name = "zeroinfl, refit", packages = "pscl",
       fit = function(m) suppressWarnings(multiScale_optim(pscl::zeroinfl(y ~ bin1 + cont2 | bin1,
         data = m$dat), m$inputs, verbose = FALSE)),
       refuse = c(refit_refusals, list(predict_conditional_se = "se.fit|standard error"))),
  list(name = "gls, refit", packages = "nlme",
       fit = function(m) multiScale_optim(nlme::gls(yg ~ bin1 + cont2, data = m$dat,
         method = "ML"), m$inputs, verbose = FALSE),
       refuse = c(refit_refusals, list(predict_conditional_se = "se.fit|standard error"))),
  list(name = "lmer, refit", packages = "lme4",
       fit = function(m) multiScale_optim(lme4::lmer(yg ~ bin1 + cont2 + (1 | grp),
         data = m$dat, REML = FALSE), m$inputs, verbose = FALSE),
       site_values = list(grp = "1"), refuse = refit_refusals),
  # nlme's `level` collides with predict()'s confidence `level`, so it goes
  # through native_args.
  list(name = "lme, refit", packages = "nlme",
       predict_args = list(native_args = list(level = 0)),
       raster_args = list(native_args = list(level = 0)),
       fit = function(m) multiScale_optim(nlme::lme(yg ~ bin1 + cont2, random = ~ 1 | grp,
         data = m$dat, method = "ML"), m$inputs, verbose = FALSE),
       refuse = c(refit_refusals, list(predict_conditional_se = "se.fit|standard error"))),
  list(name = "spaMM fitme, refit", packages = "spaMM",
       predict_args = list(), raster_args = list(),
       fit = function(m) multiScale_optim(spaMM::fitme(y ~ bin1 + cont2 + (1 | grp),
         family = poisson(), data = m$dat), m$inputs, verbose = FALSE),
       site_values = list(grp = "1"),
       refuse = c(refit_refusals, list(predict_conditional_se = "se.fit|standard error"))),
  list(name = "glmmTMB poisson, joint_mixed", packages = c("glmmTMB", "RTMB"),
       fit = function(m) multiScale_optim(glmmTMB::glmmTMB(y ~ bin1 + cont2 + (1 | grp),
         family = poisson, data = m$dat), m$inputs, backend = "joint_mixed",
         verbose = FALSE),
       site_values = list(grp = "1"),
       refuse = list(predict_joint_se = "conditional prediction SEs only",
                     vcov = "Joint covariance is unavailable for joint_mixed",
                     scale_inference = "scale_inference")),
  list(name = "glmmTMB poisson, refit", packages = "glmmTMB",
       fit = function(m) multiScale_optim(glmmTMB::glmmTMB(y ~ bin1 + cont2 + (1 | grp),
         family = poisson, data = m$dat), m$inputs, verbose = FALSE),
       site_values = list(grp = "1"), refuse = refit_refusals),
  list(name = "occu, refit", packages = "unmarked", type = "state",
       predict_args = list(type = "state"), raster_args = list(type = "state"),
       fit = function(m) multiScale_optim(unmarked::occu(~1 ~ bin1 + cont2,
         unmarked::unmarkedFrameOccu(y = m$occ$y, siteCovs = m$occ_inputs$kernel_dat)),
         m$occ_inputs, verbose = FALSE),
       refuse = refit_refusals),
  list(name = "pcount, refit", packages = "unmarked", type = "state",
       predict_args = list(type = "state"), raster_args = list(type = "state"),
       fit = function(m) {
         set.seed(7)
         y <- m$occ$y * stats::rpois(length(m$occ$y), 2)
         multiScale_optim(unmarked::pcount(~1 ~ bin1 + cont2,
           unmarked::unmarkedFramePCount(y = y, siteCovs = m$occ_inputs$kernel_dat), K = 25),
           m$occ_inputs, verbose = FALSE) },
       refuse = refit_refusals),
  list(name = "clogit, refit", packages = "survival",
       predict_args = list(type = "risk"), raster_args = list(type = "risk"),
       fit = function(m) { s <- matrix_steps_data()
         dat <- s$dat
         multiScale_optim(survival::clogit(case_ ~ bin1 + cont2 + strata(step_id),
           data = dat, model = TRUE), s$inputs, verbose = FALSE) },
       refuse = refit_refusals),
  list(name = "clogit, conditional_logit", packages = "survival",
       predict_args = list(type = "risk"), raster_args = list(type = "risk"),
       fit = function(m) { s <- matrix_steps_data()
         dat <- s$dat
         multiScale_optim(survival::clogit(case_ ~ bin1 + cont2 + strata(step_id),
           data = dat, model = TRUE), s$inputs, backend = "conditional_logit",
           verbose = FALSE) },
       refuse = list(predict_joint_se = "conditional SEs only",
                     vcov = "requires backend = 'profile_glm' or 'joint_nb2'",
                     scale_inference = "scale_inference"))
)

for (spec in matrix_models) {
  test_that(paste("integration matrix:", spec$name), {
    for (pkg in spec$packages) skip_if_not_installed(pkg)
    m <- matrix_inputs()
    spec$rs <- m$rs
    # Native fitters print notes such as lme4's singular-fit message on every
    # refit; they are not outcomes of the matrix.
    fit <- suppressMessages(spec$fit(m))
    expect_s3_class(fit, "multiScaleR")
    outcomes <- character()
    for (column in names(matrix_columns)) {
      result <- tryCatch(
        suppressWarnings(suppressMessages(matrix_columns[[column]](fit, spec))),
        error = function(e) e)
      failed <- inherits(result, "error")
      message_text <- if (failed) conditionMessage(result) else ""
      outcomes[column] <- if (failed) paste("ERROR:", substr(message_text, 1, 160)) else "ok"
      if (discover) next
      expected <- spec$refuse[[column]]
      if (is.null(expected)) {
        expect_false(failed, label = sprintf("[%s] %s failed: %s", spec$name, column,
                                             message_text))
      } else {
        expect_true(failed, label = sprintf("[%s] %s was declared a refusal (%s) but ran",
                                            spec$name, column, expected))
        if (failed) expect_match(message_text, expected,
                                 label = sprintf("[%s] %s message", spec$name, column))
      }
    }
    # A scale reported without a standard error must be explained where the
    # user reads it.
    scales <- summary(fit)$opt_scale
    if (!discover && any(!is.finite(as.data.frame(scales)$SE))) {
      printed <- paste(utils::capture.output(print(summary(fit))), collapse = " ")
      expect_match(printed, "withheld|unavailable|not estimated|hard radius|fixed|bound|profile_sigma",
                   label = sprintf("[%s] summary explains missing scale SEs", spec$name))
    }
    if (discover) {
      cat("\n== ", spec$name, "\n", sep = "")
      for (column in names(outcomes)) cat(sprintf("  %-24s %s\n", column, outcomes[[column]]))
      succeed()
    }
  })
}
