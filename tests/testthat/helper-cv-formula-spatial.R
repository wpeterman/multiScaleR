cv_formula_spatial_inputs <- function(sq = FALSE) {
  i <- 1:72; d <- seq(.025, .99, length.out = 24); ids <- paste0("formulaSpatial", i)
  values <- outer(sin(i*.73), sin(7*d)) + outer(cos(i*.41), cos(13*d)) +
    outer(sin(i*.19), d) + .1*outer(sin(i*.37), sin(19*d))
  raw <- lapply(i, function(j) matrix(values[j, ], ncol = 1, dimnames = list(NULL, "terrain")))
  distances <- lapply(i, function(j) d*(1 - (j %% 5)*.005)); names(raw) <- names(distances) <- ids
  specs <- if (sq) msr_vars(hab = surface_var("terrain", "sq", weighted = TRUE)) else msr_vars(hab = kernel_var("terrain"))
  inputs <- structure(list(kernel_dat = data.frame(hab = rep(999, length(ids)), row.names = ids), raw_cov = raw,
    d_list = distances, scale_vars = specs, kernel = "gaussian", sigma = .3, shape = NULL,
    min_D = 35, max_D = 800, unit_conv = 800, n_covs = 1L, resolution = 20, n_cols = 40L,
    binned = if (sq) NULL else .msr_build_kernel_bins(distances, raw, "hab", "terrain", 32L, ids),
    scl_params = list(mean = c(hab = 999), sd = c(hab = 100)), cell_data_stored = TRUE), class = "multiScaleR_data")
  h <- vapply(i, function(j) .msr_eval_scale_vars(distances[[j]], raw[[j]], specs, .3, NULL,
    "gaussian", 800, 20, 40L)[1L], 0)
  h <- as.numeric(scale(h))
  data <- data.frame(hab = 999, z = sin(i*.31), group = factor(rep(c("a", "b"), 36)),
    exposure = 1 + (i %% 5)/10, row.names = ids)
  eta <- .9 + .55*h + .1*data$z + .1*(data$group == "b") + log(data$exposure)
  data$count <- round(exp(eta)*rep(c(.2, .5, 1, 2.5), 18))
  data$positive <- exp(.4 + .7*h + .15*data$z + log(data$exposure))*rep(c(.4, .8, 1.5, 2), 18)
  data$continuous <- .4 + .7*h + .15*data$z + .2*cos(i*.91)
  data$binary <- as.numeric(i %% 5 < 2)
  list(data = data, inputs = inputs, folds = grouped_cv_folds(rep(letters[1:12], each = 6), 3, ids))
}


cv_formula_choice_inputs <- function() {
  # The existing profile fixture and seed are reused without any search/tuning.
  base <- profile_glm_fixture(); b <- base$inputs$binned; counts <- b$csum$hab
  ids <- paste0("formulaChoice", seq_len(nrow(base$data)))
  raw <- lapply(seq_along(ids), function(i) matrix(rep(b$vsum$hab[i, ]/counts[i, ], counts[i, ]), ncol = 1,
    dimnames = list(NULL, "hab")))
  distances <- lapply(seq_along(ids), function(i) rep(b$bin_rep, counts[i, ])); names(raw) <- names(distances) <- ids
  inputs <- structure(list(kernel_dat = data.frame(hab = base$data$hab, row.names = ids), raw_cov = raw,
    d_list = distances, scale_vars = msr_vars(hab = kernel_var("hab")), kernel = "gaussian", sigma = .25,
    shape = NULL, min_D = 50, max_D = 1000, unit_conv = 1000, cell_size = c(20, 20), resolution = 20,
    n_cols = 50L, n_covs = 1L, binned = .msr_build_kernel_bins(distances, raw, "hab", "hab", 32L, ids),
    cell_data_stored = TRUE, scl_params = list(mean = c(hab = 999), sd = c(hab = 100))), class = "multiScaleR_data")
  d <- base$data; rownames(d) <- ids; d$case_ <- rep(c(1, 0, 0), nrow(d)/3)
  d$choice <- rep(seq_len(nrow(d)/3), each = 3); folds <- grouped_cv_folds(d$choice, 3, ids)[1]
  d[folds[[1]]$assessment, "z"] <- d[folds[[1]]$assessment, "z"] + 3
  list(data = d, inputs = inputs, folds = folds)
}
