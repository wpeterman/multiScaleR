# Reuse the known conditional-logit covariate fixture. Repeated bin means are
# explicit raw cells; no new seed, response simulation or passing-seed search.
choice_spatial_cache <- new.env(parent = emptyenv())
choice_spatial_inputs <- function(kernel = "gaussian") {
  f <- profile_glm_fixture(kernel = kernel)
  b <- f$inputs$binned; ids <- paste0("choiceAlt", seq_len(nrow(f$data)))
  counts <- b$csum$hab
  raw <- lapply(seq_along(ids), function(i) matrix(rep(b$vsum$hab[i, ] / counts[i, ], counts[i, ]),
    ncol = 1L, dimnames = list(NULL, "hab")))
  distances <- lapply(seq_along(ids), function(i) rep(b$bin_rep, counts[i, ]))
  names(raw) <- names(distances) <- ids
  inputs <- structure(list(kernel_dat = data.frame(hab = f$data$hab, row.names = ids),
    raw_cov = raw, d_list = distances, scale_vars = msr_vars(hab = kernel_var("hab")),
    kernel = kernel, sigma = .25, shape = NULL, min_D = 50, max_D = 1000, unit_conv = 1000,
    cell_size = c(20, 20), resolution = 20, n_cols = 50L, n_covs = 1L,
    binned = .msr_build_kernel_bins(distances, raw, "hab", "hab", 32L, ids),
    cell_data_stored = TRUE, scl_params = list(mean = c(hab = 999), sd = c(hab = 100))), class = "multiScaleR_data")
  d <- f$data; rownames(d) <- ids
  d$case_ <- rep(c(1, 0, 0), nrow(d) / 3L); d$choice <- rep(seq_len(nrow(d)/3L), each = 3)
  d$hab <- 1e6 # Full-data values/reference must never enter training preparation.
  inputs$raster_group <- stats::setNames(rep(c("mapA", "mapB"), length.out = length(ids)), ids)
  list(data = d, inputs = inputs, folds = grouped_cv_folds(d$choice, 3, ids))
}

choice_spatial_fixture <- function() {
  if (exists("failure", choice_spatial_cache, inherits = FALSE)) stop(choice_spatial_cache$failure)
  if (!exists("value", choice_spatial_cache, inherits = FALSE)) {
    choice_spatial_cache$value <- tryCatch({
      g <- choice_spatial_inputs(); e <- choice_spatial_inputs("exp"); planned <- 0L
      form <- case_ ~ . - y - choice - exposure + offset(log(exposure)) + strata(choice)
      run <- function(fixture, folds = fixture$folds[1L], backend = "conditional_logit") {
        planned <<- planned + length(folds); stopifnot(planned <= 6L)
        arguments <- if (backend == "conditional_logit") list(profile_control = list(starts = matrix(c(.08, .25, .6), ncol = 1))) else list()
        cross_validate_clogit(form, fixture$data, folds, fixture$inputs, backend = backend,
          optim_args = arguments, keep_fits = TRUE, keep_selection = TRUE)
      }
      gaussian <- run(g, g$folds[1:2]); exponential <- run(e); refit <- run(g, backend = "refit")
      sq <- g; sq$inputs$scale_vars <- msr_vars(hab = surface_var("hab", metric = "sq", weighted = TRUE))
      sq$inputs$binned <- NULL
      weighted_sq <- run(sq)
      changed <- g; held <- g$folds[[1L]]$assessment
      for (set in unique(g$data[held, "choice"])) {
        rows <- rownames(g$data)[g$data$choice == set]; changed$data[rows, "case_"] <- c(0, 1, 0)
      }
      changed$data[held, "z"] <- changed$data[held, "z"] + 50
      changed$inputs$kernel_dat$hab <- -1e9
      changed$inputs$scl_params <- list(mean = c(hab = -1e9), sd = c(hab = .001))
      changed$inputs$binned$bin_rep[] <- 1e-5
      for (id in held) {
        changed$inputs$raw_cov[[id]] <- changed$inputs$raw_cov[[id]] * 7
        changed$inputs$d_list[[id]] <- changed$inputs$d_list[[id]] * .6
      }
      leakage <- run(changed)
      list(gaussian = gaussian, exponential = exponential, refit = refit, sq = weighted_sq,
        leakage = leakage, original = g, exp_input = e, sq_input = sq, leak_input = changed, planned = planned)
    }, error = function(e) { choice_spatial_cache$failure <- e; stop(e) })
  }
  choice_spatial_cache$value
}

test_that("spatial choice CV preserves native fits, starts, map labels and compact penalties", {
  skip_if_not_installed("survival")
  z <- choice_spatial_fixture()
  expect_identical(z$planned, 6L)
  for (case in c("gaussian", "exponential", "refit", "sq", "leakage")) {
    x <- z[[case]]
    fixture <- switch(case, exponential = z$exp_input, sq = z$sq_input, leakage = z$leak_input, z$original)
    expect_true(x$complete, info = paste(unlist(x$errors), collapse = "; "))
    if (!x$complete) next # Keep failures visible without hiding the other cases.
    compact <- x; compact$fits <- NULL
    for (criterion in c("AIC", "AICc", "BIC")) {
      a <- select_cv(list(candidate = x), criterion)
      b <- select_cv(list(candidate = compact), criterion)
      expect_true(a$complete); expect_true(b$complete)
      for (item in c("predictions", "summary", "diagnostics", "prediction_contracts", "selection", "raster_groups"))
        expect_identical(a[[item]], b[[item]])
    }
    for (name in names(x$folds)) {
      fit <- x$fits[[name]]; record <- x$selection_records[[name]]$data
      expect_equal(record$K, attr(logLik(fit$opt_mod), "df") + 1)
      expect_equal(record$nobs, length(x$folds[[name]]$analysis)/3)
      # Independent raw-cell mean/SQ reduction and frozen fitted normalization.
      # Each original distance has its own bin in this fixture, so no binned
      # approximation is being compared against a different raw-cell target.
      p <- x$predictions[x$predictions$fold == name, ]
      sigma <- fit$scale_est$Mean[1L] / fixture$inputs$unit_conv
      raw <- vapply(p$row_id, function(id) {
        d <- fixture$inputs$d_list[[id]]; v <- fixture$inputs$raw_cov[[id]][, 1L]
        # The perturbed assessment distances deliberately enter different bins.
        # Their prediction target uses frozen training representative distances,
        # never the altered raw distances or an assessment-learned reference.
        reference <- x$training[[name]]$bin_reference
        if (!is.null(reference)) d <- reference$bin_rep[.bincode(d, reference$breaks, include.lowest = TRUE)]
        w <- if (fixture$inputs$kernel == "gaussian") exp(-d^2/(2*sigma^2)) else exp(-d/sigma)
        mean <- sum(w*v)/sum(w)
        if (case == "sq") sqrt(sum(w*(v-mean)^2)/sum(w)) else mean
      }, 0)
      nd <- fixture$data[p$row_id, ]
      nd$hab <- (raw - fit$scl_params$mean[["hab"]]) / fit$scl_params$sd[["hab"]]
      nd$choice <- fixture$data[x$folds[[name]]$analysis[1L], "choice"]
      eta <- as.numeric(predict(fit$opt_mod, nd, type = "lp", reference = "zero"))
      expected <- numeric(nrow(p))
      for (g in split(seq_len(nrow(p)), p$choice)) {
        m <- max(eta[g]); expected[g] <- eta[g] - m - log(sum(exp(eta[g]-m)))
      }
      expect_equal(p$log_probability, expected, tolerance = 1e-9, info = paste(case, name))
      expect_equal(p$prediction, exp(expected), tolerance = 1e-10, info = paste(case, name))

      expect_identical(x$raster_groups[[name]], lapply(x$folds[[name]], function(rows) z$original$inputs$raster_group[rows]))
      if (identical(fit$backend, "conditional_logit")) {
        expect_equal(nrow(fit$profile_fit$start_results), 3L)
        expect_true(any(fit$profile_fit$start_results$eligible))
        expect_true(all(is.na(fit$scale_est$SE)))
        expect_false(fit$uncertainty$available)
      }
    }
  }
})

test_that("assessment changes do not alter native spatial choice training", {
  z <- choice_spatial_fixture(); a <- z$gaussian; b <- z$leakage
  expect_true(a$complete); expect_true(b$complete)
  expect_identical(coef(a$fits[[1]]$opt_mod), coef(b$fits[[1]]$opt_mod))
  expect_identical(a$fits[[1]]$scl_params, b$fits[[1]]$scl_params)
  expect_identical(a$fits[[1]]$scale_est, b$fits[[1]]$scale_est)
  expect_identical(a$training[[1]]$initial_coefficients, b$training[[1]]$initial_coefficients)
  expect_identical(a$training[[1]]$bin_reference, b$training[[1]]$bin_reference)
})

test_that("ordinary and spatial choice candidates share native event-based selection", {
  z <- choice_spatial_fixture()
  null <- cross_validate_clogit(case_ ~ strata(choice), z$original$data, z$original$folds[1:2], keep_selection = TRUE)
  models <- list(spatial = z$gaussian, null = null)
  selected <- select_cv(models, "BIC", ties = "first")
  expect_true(selected$complete)
  expect_s3_class(selected, "multiScaleR_choice_cv")
  expect_setequal(cv_tab(c(models, list(selected = selected)), "neg_log_score")$Modnames, c(names(models), "selected"))
  for (tab in selected$selection$tables) expect_equal(sort(tab$K), c(0, 4))
  expect_identical(selected$settings$selection_uncertainty, "not_integrated")
})
