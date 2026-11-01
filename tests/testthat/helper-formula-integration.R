# Existing numerical inputs only; no recovery or sampler panel.
formula_integration_fixture <- function(family = "poisson", two = FALSE) {
  f <- profile_glm_fixture()
  if (two) {
    b <- f$inputs$binned
    b$vsum$cover <- b$vsum$hab * rep(1 + sin(3 * b$bin_rep), each = nrow(f$data)) +
      outer(sin(seq_len(nrow(f$data)) * .17), b$bin_rep) * b$csum$hab
    b$csum$cover <- b$csum$hab
    f$inputs$binned <- b
    engine <- multiScaleR:::.msr_prepare_kernel_engine(b, c("hab", "cover"))
    f$data$cover <- multiScaleR:::.msr_eval_kernel_engine(engine, log(c(.25, .4)))$value[, "cover"]
  }
  response <- "y"
  fam <- switch(family, poisson = poisson(), gaussian = gaussian(),
    nbinom2 = MASS::negative.binomial(2), binomial = binomial())
  if (family == "gaussian") f$data$y <- log1p(f$data$y)
  if (family == "binomial") {
    f$data$success <- f$data$y %% 4L
    f$data$trials <- 5L
    response <- "cbind(success, trials - success)"
  }
  rhs <- if (two) "hab * cover + group + z" else "hab * group + z"
  f$model <- glm(as.formula(paste(response, "~", rhs, "+ offset(log(exposure))")),
    data = f$data, family = fam)
  f
}

# Reuse the established weighted-SQ inputs without relying on test-file order.
formula_weighted_fixture <- function(kernel = "gaussian", surface_only = FALSE) {
  f <- profile_glm_fixture(n = 120, kernel = kernel)
  b <- f$inputs$binned
  f$inputs$scale_vars <- msr_vars(hab = kernel_var("hab"),
    rough = surface_var("elev", "sq", weighted = TRUE))
  f$inputs$d_list <- rep(list(b$bin_rep), nrow(f$data))
  f$inputs$raw_cov <- lapply(seq_len(nrow(f$data)), function(i) {
    v <- b$vsum$hab[i, ] / b$csum$hab[i, ]
    data.frame(hab = v, elev = v^2 + sin(i * b$bin_rep))
  })
  f$inputs$raw_cov[[1]]$elev[1] <- NA_real_
  f$inputs$kernel_dat <- data.frame(row.names = rownames(f$data))
  if (surface_only) {
    f$inputs$scale_vars <- f$inputs$scale_vars[2, , drop = FALSE]
    f$inputs$binned <- NULL
  }
  covs <- if (surface_only) "rough" else c("hab", "rough")
  e <- multiScaleR:::.msr_prepare_surface_engine(f$inputs, covs, seq_len(nrow(f$data)))
  f$parameters <- log(if (surface_only) .4 else c(.25, .4))
  kd <- multiScaleR:::.msr_eval_kernel_engine(e, f$parameters)
  for (nm in covs) f$data[[nm]] <- kd$value[, nm]
  f$formula <- if (surface_only) y ~ rough * group + z + offset(log(exposure)) else
    y ~ hab * rough + rough * group + z + offset(log(exposure))
  f
}
