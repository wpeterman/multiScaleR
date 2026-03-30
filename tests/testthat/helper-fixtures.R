make_core_fixture <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    rs_all <- sim_rast(dim = 50, resolution = 10, user_seed = 42)
    rs <- terra::subset(rs_all, "cont1")

    sim <- sim_dat(
      alpha = 0.5,
      beta = 0.75,
      kernel = "gaussian",
      sigma = 40,
      type = "count",
      n_points = 15,
      raster_stack = rs,
      max_D = 250,
      user_seed = 42
    )

    kernel_inputs <- kernel_prep(
      pts = sim$pts,
      raster_stack = rs,
      max_D = 250,
      kernel = "gaussian",
      verbose = FALSE
    )

    df <- transform(sim$df, site = seq_len(nrow(sim$df)) / 10)

    fitted_mod <- glm(y ~ cont1 + site, family = poisson(), data = df)
    opt <- suppressWarnings(
      suppressMessages(
        multiScale_optim(
          fitted_mod = fitted_mod,
          kernel_inputs = kernel_inputs,
          par = 40 / kernel_inputs$unit_conv,
          n_cores = NULL,
          verbose = FALSE
        )
      )
    )

    plain_mod <- glm(y ~ site, family = poisson(), data = df)

    factor_df <- transform(
      df,
      group = factor(rep(c("a", "b", "c"), length.out = nrow(df)))
    )
    factor_mod <- glm(y ~ cont1 + group, family = poisson(), data = factor_df)
    opt_factor <- opt
    opt_factor$opt_mod <- factor_mod

    cache <<- list(
      rs_all = rs_all,
      rs = rs,
      sim = sim,
      kernel_inputs = kernel_inputs,
      df = df,
      fitted_mod = fitted_mod,
      opt = opt,
      plain_mod = plain_mod,
      opt_factor = opt_factor
    )

    cache
  }
})

make_simulation_fixture <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    rs <- sim_rast(dim = 50, resolution = 10, user_seed = 99)
    rs_two <- terra::subset(rs, c("bin1", "cont1"))

    sim <- sim_dat(
      alpha = 0.5,
      beta = c(0.5, -0.25),
      kernel = "gaussian",
      sigma = c(40, 60),
      type = "count",
      n_points = 15,
      raster_stack = rs_two,
      max_D = 250,
      user_seed = 99
    )

    sim_umf <- sim_dat_unmarked(
      alpha = 0.5,
      beta = c(0.5, -0.25),
      kernel = "gaussian",
      sigma = c(40, 60),
      n_points = 15,
      n_surv = 4,
      det = 0.6,
      type = "count",
      raster_stack = rs_two,
      max_D = 250,
      user_seed = 99
    )

    cache <<- list(
      rs = rs,
      rs_two = rs_two,
      sim = sim,
      sim_umf = sim_umf
    )

    cache
  }
})

make_unmarked_fixture <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    rs <- sim_rast(dim = 50, resolution = 10, user_seed = 99)
    rs_two <- terra::subset(rs, c("bin1", "cont1"))

    sim <- sim_dat_unmarked(
      alpha = 0.5,
      beta = c(0.5, -0.25),
      kernel = "gaussian",
      sigma = c(40, 60),
      n_points = 15,
      n_surv = 4,
      det = 0.6,
      type = "count",
      raster_stack = rs_two,
      max_D = 250,
      user_seed = 99
    )

    kernel_inputs <- kernel_prep(
      pts = sim$pts,
      raster_stack = rs_two,
      max_D = 250,
      kernel = "gaussian",
      verbose = FALSE
    )

    site_covs <- data.frame(
      site = seq_len(nrow(kernel_inputs$kernel_dat)),
      nuisance = seq_len(nrow(kernel_inputs$kernel_dat)) / 10,
      kernel_inputs$kernel_dat
    )

    umf <- unmarked::unmarkedFramePCount(y = sim$y, siteCovs = site_covs)
    mod <- unmarked::pcount(~1 ~ bin1 + cont1, data = umf, K = 50)

    cache <<- list(
      rs = rs,
      rs_two = rs_two,
      sim = sim,
      kernel_inputs = kernel_inputs,
      site_covs = site_covs,
      umf = umf,
      mod = mod,
      obj = structure(
        list(opt_mod = mod, scl_params = kernel_inputs$scl_params),
        class = "multiScaleR"
      )
    )

    cache
  }
})

make_zeroinfl_fixture <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    set.seed(123)
    dat <- data.frame(
      y = rpois(80, lambda = 3),
      x = rnorm(80),
      z = rnorm(80)
    )

    mod <- pscl::zeroinfl(y ~ x | z, data = dat, dist = "poisson")

    cache <<- list(
      dat = dat,
      mod = mod,
      obj = structure(
        list(opt_mod = mod, scl_params = list(mean = c(x = 0), sd = c(x = 1))),
        class = "multiScaleR"
      )
    )

    cache
  }
})

open_null_device <- function() {
  grDevices::pdf(file = NULL)
}
