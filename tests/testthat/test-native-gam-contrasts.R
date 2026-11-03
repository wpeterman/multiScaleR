# Native full-context oracles keep all fitted factor levels represented.
# No outer ecological-scale optimization or recovery study is run here.
gam_contrast_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    f <- refit_spatial_fixture(); dat <- f$data
    dat$y <- dat$y + as.integer(dat$site)^2
    dat$flag <- dat$z > 0
    old <- options(contrasts = c("contr.sum", "contr.poly"))
    on.exit(options(old), add = TRUE)
    slice <- function(value, rows) {
      if (is.list(value)) {
        result <- lapply(value, slice, rows = rows)
        attributes(result) <- attributes(value)
        return(result)
      }
      if (is.matrix(value)) result <- value[rows, , drop = FALSE] else if
          (length(dim(value)) == 1L) {
        labels <- dimnames(value)
        if (!is.null(labels)) labels <- list(labels[[1L]][rows])
        result <- array(as.numeric(value)[rows], length(rows), dimnames = labels)
      } else if (length(value) == nrow(dat) || length(value) > nrow(dat)) {
        result <- value[rows]
      } else return(value)
      for (name in setdiff(names(attributes(value)), c("dim", "dimnames", "names")))
        attr(result, name) <- attr(value, name)
      result
    }
    models <- list()
    for (kind in c("single", "single_factor", "gaulss", "gaulss_factor",
        "single_logical", "gaulss_logical")) {
      transformed <- grepl("factor", kind, fixed = TRUE)
      group <- if (identical(kind, "single_logical")) "flag" else
        if (grepl("logical", kind, fixed = TRUE)) "I(z > 0)" else
        if (transformed) "factor(group)" else "group"
      first <- as.formula(paste("y ~ rough + z + z:", group, "+ s(mean,k=4)"),
        env = asNamespace("mgcv"))
      second <- as.formula(paste("~ fixed + fixed:", group), env = asNamespace("mgcv"))
      multi <- startsWith(kind, "gaulss")
      model <- mgcv::gam(if (multi) list(first, second) else first,
        data = dat, method = "REML", family = if (multi) mgcv::gaulss() else gaussian())
      expect_identical(model$outer.info$conv, "full convergence")
      expect_true(all(eigen(model$outer.info$hess, symmetric = TRUE,
        only.values = TRUE)$values > 0))
      source <- dat[c("rough", "z", "group", "mean", "fixed", "flag")]
      outputs <- list()
      for (type in c("link", "response", "terms", "lpmatrix")) {
        outputs[[type]] <- predict(model, source, type = type, se.fit = type != "lpmatrix")
      }
      # Augment a full context to obtain independent single-level slice oracles.
      map_data <- source; map_data$group <- factor("b", levels = levels(dat$group))
      map_data$z <- .2
      map_data$flag <- TRUE
      map_full <- predict(model, rbind(source, map_data), type = "response", se.fit = TRUE)
      map_expected <- slice(map_full, nrow(dat) + seq_len(nrow(dat)))
      defaults <- lapply(source, function(v) if (is.numeric(v)) mean(v) else v[1L])
      curve_data <- as.data.frame(lapply(defaults, rep, times = 5), check.names = FALSE)
      curve_data$z <- seq(min(dat$z), max(dat$z), length.out = 5)
      curve_full <- predict(model, rbind(source, curve_data), type = "response", se.fit = TRUE)
      curve_expected <- slice(curve_full, nrow(dat) + seq_len(5))
      selected <- attr(if (multi) model$pterms[[1L]] else model$pterms, "term.labels")[1L]
      control <- list(terms = selected, block.size = 1)
      if (!is.null(model$Vc)) control$unconditional <- TRUE
      expected_control <- control; expected_control$block.size <- nrow(source)
      selected_expected <- do.call(predict, c(list(object = model, newdata = source,
        type = "link", se.fit = TRUE), expected_control))
      models[[kind]] <- list(model = model, multi = multi, outputs = outputs,
        map = map_expected, curve = curve_expected, controls = control,
        selected = selected_expected,
        site_values = if (identical(kind, "single_logical")) list(z = .2, flag = TRUE) else
          if (grepl("logical", kind, fixed = TRUE)) list(z = .2) else
          list(group = "b", z = .2), native_df = as.numeric(attr(logLik(model), "df")))
    }
    cache <<- list(f = f, data = dat, models = models, slice = slice)
    cache
  }
})

test_that("native GAM coding stays fixed across package prediction routes", {
  skip_if_not_installed("mgcv")
  f <- gam_contrast_fixture(); dat <- f$data
  old_options <- options(contrasts = c("contr.treatment", "contr.helmert"))
  on.exit(options(old_options), add = TRUE)
  provider <- utils::getS3method("predict", "gam", envir = asNamespace("mgcv"))
  original_body <- body(provider); original_environment <- environment(provider)
  compare <- function(actual, expected, exact_rows = FALSE) {
    if (is.list(expected)) {
      expect_identical(names(actual), names(expected))
      for (name in names(expected)) compare(actual[[name]], expected[[name]], exact_rows)
    } else {
      expect_identical(dim(actual), dim(expected))
      if (exact_rows) expect_identical(dimnames(actual), dimnames(expected))
      else if (length(dim(expected)) > 1L) expect_identical(colnames(actual), colnames(expected))
      expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-8)
      expect_true(all(is.finite(actual)))
    }
  }
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  for (kind in names(f$models)) {
    case <- f$models[[kind]]; model <- case$model; original <- model
    object <- f$f$object; object$opt_mod <- model
    expect_equal(.msr_selection_parameter_count(object), case$native_df + 2L, tolerance = 1e-8)
    for (rows in list(seq_len(nrow(dat)), c(1L,17L,24L), 17L)) {
      nd <- dat[rows, , drop = FALSE]
      for (type in names(case$outputs)) {
        expected <- f$slice(case$outputs[[type]], rows)
        direct <- .msr_predict_gam(model, newdata = nd, type = type,
          se.fit = type != "lpmatrix", block.size = 1)
        compare(direct, expected, TRUE)
        actual <- predict(object, nd, f$f$inputs, include_fixed = TRUE,
          type = type, se.fit = type != "lpmatrix", uncertainty = "conditional",
          native_args = list(block.size = 1))
        compare(actual, expected)
      }
    }
    selected <- predict(object, dat, f$f$inputs, include_fixed = TRUE, type = "link",
      se.fit = TRUE, uncertainty = "conditional", native_args = case$controls)
    compare(selected, case$selected)
    rows <- c(1L,17L,24L); nd <- dat[rows, , drop = FALSE]
    snap <- prepare_prediction(object, nd, f$f$inputs, include_fixed = TRUE,
      component = "response", native_args = list(block.size = 1))
    snap <- unserialize(serialize(snap, NULL))
    expect_null(snap$opt_mod)
    compare(predict(snap, se.fit = TRUE, uncertainty = "conditional"),
      f$slice(case$outputs$response, rows))
    grid <- kernel_scale.raster(f$f$rasters, multiScaleR = object,
      scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE,
      site_values = case$site_values, verbose = FALSE)
    map <- predict(object, grid, type = "response", se.fit = TRUE, uncertainty = "conditional")
    expect_equal(unname(as.matrix(terra::extract(map, f$f$points, ID = FALSE))),
      unname(cbind(case$map$fit, case$map$se.fit)), tolerance = 1e-8)
    reference <- dat; row.names(reference) <- row.names(model$model)
    plots <- plot_marginal_effects(object, reference_data = reference, link = TRUE,
      component = if (case$multi) 1L else NULL, length.out = 5)
    expected_fit <- if (case$multi) case$curve$fit[,1L] else as.numeric(case$curve$fit)
    expected_se <- if (case$multi) case$curve$se.fit[,1L] else as.numeric(case$curve$se.fit)
    expect_equal(plots$z$data$fit, unname(expected_fit), tolerance = 1e-8)
    expect_equal(plots$z$data$lwr, unname(expected_fit + qnorm(.025)*expected_se), tolerance = 1e-8)
    expect_equal(plots$z$data$upr, unname(expected_fit + qnorm(.975)*expected_se), tolerance = 1e-8)
    expect_identical(model, original)
  }
  expect_identical(utils::getS3method("predict","gam",envir=asNamespace("mgcv")), provider)
  expect_identical(body(provider), original_body)
  expect_identical(environment(provider), original_environment)
})

test_that("GAM coding provenance and provider boundaries fail explicitly", {
  skip_if_not_installed("mgcv")
  f <- gam_contrast_fixture(); dat <- f$data
  old_options <- options(contrasts = c("contr.treatment", "contr.helmert"))
  on.exit(options(old_options), add = TRUE)
  case <- f$models$gaulss; object <- f$f$object; object$opt_mod <- case$model
  reference <- dat; row.names(reference) <- row.names(case$model$model)
  contrasts(reference$group) <- contr.treatment(3)
  expect_error(plot_marginal_effects(object, reference_data = reference,
    component = 1L), "contrasts differ")
  model <- f$models$gaulss_factor$model
  bad <- dat[1:3, , drop = FALSE]; bad$group <- factor("new")
  expect_error(.msr_predict_gam(model, newdata = bad), "new levels")
  fake <- function(object, newdata, ...) rep(0, nrow(newdata))
  expect_error(.msr_predict_gam_saved(case$model, list(newdata = dat), fake),
    "Cannot verify")
  corrupt <- case$model
  corrupt$contrasts <- list(group = "contr.sum", group = "contr.treatment")
  expect_error(.msr_predict_gam(corrupt, newdata = dat), "conflicting")
})
