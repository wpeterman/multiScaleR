# Six fixed software fixture fits, with failures cached rather than retried.
# No recovery, sampler, profile, timing or scientific interval admission panel.
unmarked_trained_cache <- new.env(parent = emptyenv())
unmarked_trained_fixture <- function() {
  if (exists("failure", unmarked_trained_cache, inherits = FALSE)) stop(unmarked_trained_cache$failure)
  if (!exists("value", unmarked_trained_cache, inherits = FALSE)) {
    unmarked_trained_cache$value <- tryCatch({
      had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
      if (had) seed <- get(".Random.seed", .GlobalEnv)
      on.exit(if (had) assign(".Random.seed", seed, .GlobalEnv) else
        if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
      set.seed(712)
      site <- factor(rep(seq_len(30L), 3L)); Year <- factor(rep(1:3, each = 30L))
      habitat <- rep(rnorm(30L), 3L)
      temperature <- runif(length(site), 10, 20); effort <- runif(length(site), .8, 1.2)
      N <- rpois(length(site), exp(1.1 + .4 * habitat + .2 * (as.integer(Year) - 1)))
      y <- replicate(3L, rbinom(length(N), N, .7))
      sc <- data.frame(site, effort, temperature, Year, unused = seq_along(site), habitat)
      set.seed(907); present <- rbinom(nrow(sc), 1, .65)
      yo <- replicate(3L, present * rbinom(nrow(sc), 1, .6))
      f <- ~I(scale(temperature)^2) + Year ~ I(base::scale(habitat)^2) + Year +
        offset(.02 * base::scale(effort))
      mixed <- ~I(scale(temperature)^2) + Year ~ I(base::scale(habitat)^2) + Year +
        offset(.02 * base::scale(effort)) + (1 | site)
      fits <- warnings <- snapshots <- list(); calls <- 0L
      for (name in c("P", "NB", "ZIP", "occu", "P_TMB", "NB_TMB")) {
        data <- sc; response <- if (name == "occu") yo else y
        if (name == "P") contrasts(data$Year) <- contr.sum(3L)
        if (name %in% c("NB", "occu")) response[4L, ] <- NA
        frame <- if (name == "occu") unmarked::unmarkedFrameOccu(response, siteCovs = data) else
          unmarked::unmarkedFramePCount(response, siteCovs = data)
        tmb <- grepl("_TMB$", name); messages <- character()
        calls <- calls + 1L; stopifnot(calls <= 6L)
        fit <- withCallingHandlers(if (name == "occu") unmarked::occu(f, frame, engine = "C") else
          unmarked::pcount(if (tmb) mixed else f, frame, K = 40,
            mixture = sub("_TMB$", "", name), engine = if (tmb) "TMB" else "C", threads = 1L),
          warning = function(w) {messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")})
        fits[[name]] <- fit; warnings[[name]] <- messages
        snapshots[[name]] <- list(estimates = fit@estimates, opt = fit@opt,
          negLogLike = fit@negLogLike, AIC = fit@AIC, formula = fit@formula,
          formlist = fit@formlist, data = fit@data, TMB_par = fit@TMB$par)
      }
      list(fits = fits, warnings = warnings, snapshots = snapshots, calls = calls, rows = c(2, 8, 17, 31, 48, 69, 86))
    }, error = function(e) {unmarked_trained_cache$failure <- e; stop(e)})
  }
  unmarked_trained_cache$value
}

# An independent numeric normalization reference. Detection training covariates
# have three observation rows per site, so their sample SD differs from state.
unmarked_trained_reference <- function(model, type, rows, backTransform = TRUE, level = .95, re.form = NULL) {
  raw <- unmarked:::get_covariates(model, type)
  scaled <- function(x, original) (x - mean(original))/sd(original)
  rows$sH <- scaled(rows$habitat, raw$habitat)
  rows$sT <- scaled(rows$temperature, raw$temperature)
  rows$sE <- scaled(rows$effort, raw$effort)
  form <- if (type == "state") ~I(sH^2) + Year else ~I(sT^2) + Year
  if (!is.null(attr(raw$Year, "contrasts"))) contrasts(rows$Year) <- attr(raw$Year, "contrasts")
  X <- model.matrix(form, rows, na.action = na.pass)
  offset <- if (type == "state") .02 * rows$sE else NULL
  if (is.null(re.form) && unmarked:::has_random(model@formlist[[type]]))
    X <- cbind(X, unmarked:::get_Z(model@formlist[[type]], raw, rows))
  unmarked:::predict_by_chunk(model, type, level, X, offset, 70L, backTransform, re.form)
}

test_that("unmarked trained point tables retain native component targets and SEs", {
  f <- unmarked_trained_fixture()
  expect_identical(f$calls, 6L)
  provider <- body(unmarked:::make_mod_matrix)
  for (name in names(f$fits)) {
    model <- f$fits[[name]]; object <- structure(list(opt_mod = model), class = "multiScaleR")
    expect_identical(model@opt$convergence, 0L)
    k <- .msr_parameter_count(model)
    expect_equal(k, c(P = 8, NB = 9, ZIP = 9, occu = 8, P_TMB = 9, NB_TMB = 10)[[name]])
    if (name %in% c("NB", "occu")) expect_match(f$warnings[[name]], "Site.*4.*dropped") else
      expect_identical(f$warnings[[name]], character())
    rows <- model@data@siteCovs[f$rows, , drop = FALSE]
    predictor <- .msr_unmarked_prediction_function(model)
    for (type in c("state", "det")) for (nd in list(rows, rows[1L, , drop = FALSE], rows[c(2, 2, 1), , drop = FALSE])) {
      for (back in c(TRUE, FALSE)) {
        expected <- unmarked_trained_reference(model, type, nd, backTransform = back, level = .9)
        actual <- predictor(model, type, nd, backTransform = back, level = .9)
        expect_equal(actual, expected, tolerance = 1e-10)
        point <- .msr_native_prediction_rows(object, nd, type, TRUE, list(backTransform = back, level = .9))
        expect_equal(as.data.frame(point), expected, tolerance = 1e-10, ignore_attr = TRUE)
        expect_identical(attr(point, "prediction_contract")$uncertainty$kernel_parameters, "fixed_at_estimate")
      }
    }
    expect_identical(.msr_parameter_count(model), k)
    before <- f$snapshots[[name]]
    expect_identical(model@estimates, before$estimates)
    expect_identical(model@opt, before$opt)
    expect_identical(model@negLogLike, before$negLogLike)
    expect_identical(model@AIC, before$AIC)
    expect_identical(model@formlist, before$formlist)
    expect_identical(model@data, before$data)
    expect_identical(model@TMB$par, before$TMB_par)
  }
  expect_identical(body(unmarked:::make_mod_matrix), provider)
})

test_that("unmarked frames preserve training bases and native observation layout", {
  f <- unmarked_trained_fixture()
  for (model in f$fits) {
    predictor <- .msr_unmarked_prediction_function(model)
    frame <- model@data[f$rows, ]
    for (type in c("state", "det")) {
      rows <- unmarked:::get_covariates(frame, type)
      expected <- unmarked_trained_reference(model, type, rows)
      actual <- predictor(model, type, frame)
      expect_equal(actual, expected, tolerance = 1e-10)
      expect_equal(nrow(actual), length(f$rows) * if (type == "det") 3L else 1L)
      if (length(model@sitesRemoved)) expect_warning(original <- predictor(model, type), "Site.*dropped") else
        original <- predictor(model, type)
      raw <- unmarked:::get_covariates(model, type)
      keep <- if (type == "det") rep(!seq_len(nrow(model@data@y)) %in% model@sitesRemoved, each = 3L) else
        !seq_len(nrow(model@data@y)) %in% model@sitesRemoved
      expect_equal(original, unmarked_trained_reference(model, type, raw[keep, , drop = FALSE]), tolerance = 1e-10)
    }
  }
})

test_that("unmarked snapshots and raster chunks reuse the frozen training basis", {
  f <- unmarked_trained_fixture()
  for (name in names(f$fits)) {
    model <- f$fits[[name]]; object <- structure(list(opt_mod = model), class = "multiScaleR")
    rows <- model@data@siteCovs[f$rows, , drop = FALSE]
    for (type in c("state", "det")) {
      expected <- unmarked_trained_reference(model, type, rows)
      snap <- prepare_prediction(object, rows, component = type)
      restored <- unserialize(serialize(snap, NULL))
      expect_equal(predict(restored), expected, tolerance = 1e-10)
      expect_lt(as.numeric(object.size(restored)), 128 * 1024)
      expect_false(any(vapply(restored, is.environment, logical(1))))
      thin <- list(object = list(opt_mod = model), type = type, se.fit = TRUE,
        native_args = list(), logical_names = character(), unmarked_predictor = .msr_unmarked_prediction_function(model))
      value <- .msr_native_raster_eval(thin, rows)
      expect_equal(as.data.frame(value), expected, tolerance = 1e-10, ignore_attr = TRUE)
      chunks <- rbind(.msr_native_raster_eval(thin, rows[1:3, , drop = FALSE]),
        .msr_native_raster_eval(thin, rows[4:7, , drop = FALSE]))
      expect_equal(as.matrix(chunks), as.matrix(value), tolerance = 1e-10, ignore_attr = TRUE)
    }
  }
})

test_that("unmarked reconstruction refuses changed or missing original designs", {
  f <- unmarked_trained_fixture(); model <- f$fits$P
  bad <- model; bad@data@siteCovs$habitat[1L] <- bad@data@siteCovs$habitat[1L] + .1
  expect_error(.msr_unmarked_prediction_function(bad), "do not match")
  bad <- model; bad@data@siteCovs$habitat <- NULL
  expect_error(.msr_unmarked_prediction_function(bad), "Original unmarked covariates")
  bad <- model; bad@formlist$state <- ~I(base::scale(habitat)^2) + Year + (scale(effort) | site)
  expect_error(.msr_unmarked_prediction_function(bad), "Learned random-effect")
  bad <- model; bad@nllFun <- function(x) 0
  expect_error(.msr_unmarked_prediction_function(bad), "Original fitted")
  predictor <- .msr_unmarked_prediction_function(model)
  rows <- model@data@siteCovs[f$rows, , drop = FALSE]
  rows$Year <- factor("new")
  expect_error(predictor(model, "state", rows), "new level")
  rows$Year <- model@data@siteCovs$Year[f$rows]
  rows$habitat <- NULL
  expect_error(.msr_native_prediction_rows(structure(list(opt_mod = model), class = "multiScaleR"), rows, "state", TRUE), "explicitly contain")
})

test_that("all nonnull native unmarked re.form controls mean zero effects", {
  f <- unmarked_trained_fixture()
  for (name in c("P_TMB", "NB_TMB")) {
    model <- f$fits[[name]]; rows <- model@data@siteCovs[f$rows, , drop = FALSE]
    object <- structure(list(opt_mod = model), class = "multiScaleR")
    for (control in list(NA, ~0, ~(1 | site))) {
      expected <- unmarked_trained_reference(model, "state", rows, re.form = control)
      rows$site <- NULL
      actual <- .msr_native_prediction_rows(object, rows, "state", TRUE, list(re.form = control))
      expect_equal(as.data.frame(actual), expected, tolerance = 1e-10, ignore_attr = TRUE)
      contract <- attr(actual, "prediction_contract")
      expect_identical(contract$effects$response, "zero_deviation_reference")
      expect_identical(contract$integration$random_effects, "set_to_zero_not_integrated")
    }
  }
})

test_that("all four unmarked fixed learners retain trained interactions without fitting", {
  f <- unmarked_trained_fixture(); data <- f$fits$P@data@siteCovs
  expressions <- c("base::scale(habitat)", "stats::poly(habitat, 2)",
    "splines::ns(habitat, df = 3)", "splines::bs(habitat, df = 3)")
  for (text in expressions) {
    learner <- eval(str2lang(text), data)
    formula <- as.formula(paste("~I(scale(temperature)) ~ Year * I(", text,
      ") + offset(.02 * base::scale(effort))"))
    forms <- unmarked:::split_formula(formula); names(forms) <- c("det", "state")
    frame <- unmarked::unmarkedFramePCount(matrix(0, nrow(data), 3), siteCovs = data)
    dm <- unmarked:::getDesign(frame, forms)
    likelihood <- local({original <- dm; dm <- original; function(p) 0})
    # This metadata-only prototype has no fit, estimates or covariance.
    prototype <- methods::new("unmarkedFitPCount", data = frame, formula = formula,
      formlist = forms, nllFun = likelihood, sitesRemoved = integer(), K = 40, mixture = "P")
    basis <- .msr_unmarked_prediction_basis(prototype, "state")
    expect_identical(basis$columns, colnames(dm$X_state))
    expect_equal(.msr_unmarked_fixed_design(basis, data)$X, dm$X_state, tolerance = 1e-12)
    for (indices in list(f$rows, f$rows[1L], f$rows[c(1L, 1L, 3L)])) {
      rows <- data[indices, , drop = FALSE]
      reference <- if (grepl("scale", text, fixed = TRUE))
        (rows$habitat - mean(data$habitat))/sd(data$habitat) else predict(learner, rows$habitat)
      rows$basis <- reference
      oracle <- model.matrix(~Year * I(basis), rows)
      value <- .msr_unmarked_fixed_design(basis, data[indices, , drop = FALSE])
      expect_equal(unname(value$X), unname(oracle), tolerance = 1e-12, ignore_attr = TRUE)
      expect_equal(as.numeric(value$offset),
        .02 * (rows$effort - mean(data$effort))/sd(data$effort), tolerance = 1e-12)
    }
  }
})

test_that("unmarked appendData and missing-site controls retain native semantics", {
  f <- unmarked_trained_fixture()
  for (model in f$fits) {
    predictor <- .msr_unmarked_prediction_function(model)
    rows <- model@data@siteCovs[f$rows, , drop = FALSE]
    for (type in c("state", "det")) {
      reference <- unmarked_trained_reference(model, type, rows)
      appended <- predictor(model, type, rows, appendData = TRUE)
      attr(reference, "row.names") <- attr(rows, "row.names")
      expect_equal(appended[, names(reference)], reference, tolerance = 1e-10)
      expect_identical(rownames(appended), rownames(rows))
      expect_equal(appended[, names(rows)], rows, ignore_attr = TRUE)
      raw <- unmarked:::get_covariates(model, type)
      expected <- unmarked_trained_reference(model, type, raw)
      expect_equal(predictor(model, type, na.rm = FALSE), expected, tolerance = 1e-10)
    }
  }
})

test_that("unmarked zero-effect frames omit random-only inputs and retain row layout", {
  f <- unmarked_trained_fixture()
  for (name in c("P_TMB", "NB_TMB")) {
    model <- f$fits[[name]]; predictor <- .msr_unmarked_prediction_function(model)
    frame <- model@data[f$rows, ]; frame@siteCovs$site <- NULL
    expect_error(predictor(model, "state", frame), "explicitly contain.*site")
    for (type in c("state", "det")) for (control in list(NA, ~0, ~(1 | site))) {
      raw <- unmarked:::get_covariates(frame, type)
      expected <- unmarked_trained_reference(model, type, raw, re.form = control)
      expect_equal(predictor(model, type, frame, re.form = control), expected, tolerance = 1e-10)
      object <- structure(list(opt_mod = model), class = "multiScaleR")
      snap <- prepare_prediction(object, frame, component = type, native_args = list(re.form = control))
      expect_equal(predict(unserialize(serialize(snap, NULL))), expected, tolerance = 1e-10)
      expect_identical(snap$prediction_contract[[type]]$effects$response, "zero_deviation_reference")
    }
  }
})


test_that("unmarked fitted contrast coding survives changed session defaults", {
  f <- unmarked_trained_fixture()
  old <- options("contrasts"); on.exit(options(old), add = TRUE)
  for (model in f$fits) for (type in c("state", "det")) {
    options(old)
    rows <- model@data@siteCovs[f$rows, , drop = FALSE]
    expected <- unmarked_trained_reference(model, type, rows)
    options(contrasts = c("contr.helmert", "contr.poly"))
    actual <- .msr_unmarked_prediction_function(model)(model, type, rows)
    expect_equal(actual, expected, tolerance = 1e-10)
    expect_identical(options("contrasts")$contrasts, c("contr.helmert", "contr.poly"))
  }
})

test_that("ordinary factor and random-only designs enter the checked adapter", {
  data <- data.frame(x = sin(seq_len(18)), group = factor(rep(1:6, 3)),
    Year = factor(rep(1:3, each = 6)))
  frame <- unmarked::unmarkedFramePCount(matrix(0, 18, 3), siteCovs = data)
  for (formula in list(~1 ~ Year + x, ~1 ~x + (1 | group))) {
    forms <- unmarked:::split_formula(formula); names(forms) <- c("det", "state")
    dm <- unmarked:::getDesign(frame, forms)
    likelihood <- local({original <- dm; dm <- original; function(p) 0})
    prototype <- methods::new("unmarkedFitPCount", data = frame, formula = formula,
      formlist = forms, nllFun = likelihood, sitesRemoved = integer(), K = 20, mixture = "P")
    expect_false(identical(.msr_unmarked_prediction_function(prototype), unmarked::predict))
  }
})


test_that("lost main-factor contrast metadata is recovered from original fitted values", {
  old <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(old), add = TRUE)
  raw <- data.frame(Year = factor(rep(1:3, each = 6)), x = sin(seq_len(18)))
  tt <- terms(~Year * x); frame <- model.frame(tt, raw)
  full <- model.matrix(tt, frame)
  keep <- seq_len(nrow(raw)) != 4L
  fitted <- full[keep, , drop = FALSE]
  expect_null(attr(fitted, "contrasts"))
  options(contrasts = c("contr.helmert", "contr.poly"))
  contrast <- .msr_unmarked_training_contrasts(tt, frame, fitted, keep)
  expect_equal(unname(contrast$Year), unname(contr.treatment(3)), tolerance = 0)
  basis <- list(terms = tt, xlevels = list(Year = levels(raw$Year)),
    contrasts = contrast, columns = colnames(fitted))
  expect_equal(unname(.msr_unmarked_fixed_design(basis, raw)$X), unname(full),
    tolerance = 0, ignore_attr = TRUE)
  changed <- fitted; changed[1L, 2L] <- changed[1L, 2L] + .1
  expect_error(.msr_unmarked_training_contrasts(tt, frame, changed, keep), "unambiguously")
  absent <- raw$Year != "3"
  expect_error(.msr_unmarked_training_contrasts(tt, frame, full[absent, , drop = FALSE], absent), "unambiguously")
  interaction <- terms(~Year:x)
  mf <- model.frame(interaction, raw); X <- model.matrix(interaction, mf)
  attr(X, "contrasts") <- NULL
  expect_error(.msr_unmarked_training_contrasts(interaction, mf, X, rep(TRUE, nrow(raw))), "unambiguously")
})


# Two additional fixed SOFTWARE fits. These do not assess recovery or interval
# calibration. Native warnings/failures are retained and never retried.
unmarked_zero_slope_cache <- new.env(parent = emptyenv())
unmarked_zero_slope_fixture <- function() {
  if (exists("failure", unmarked_zero_slope_cache, inherits = FALSE)) stop(unmarked_zero_slope_cache$failure)
  if (!exists("value", unmarked_zero_slope_cache, inherits = FALSE)) {
    unmarked_zero_slope_cache$value <- tryCatch({
      inherited <- unmarked_trained_fixture()
      data <- inherited$fits$NB@data@siteCovs
      contrasts(data$Year) <- contr.treatment(3L)
      had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
      if (had) seed <- get(".Random.seed", .GlobalEnv)
      on.exit(if (had) assign(".Random.seed", seed, .GlobalEnv) else
        if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
      set.seed(9713)
      data$wind <- sin(seq_len(nrow(data)) * .31) + cos(seq_len(nrow(data)) * .17)
      standard <- function(x) as.numeric(base::scale(x))
      field <- rnorm(30L, sd = .55)
      mu <- exp(.85 + .12 * standard(data$habitat)^2 + .15 * (as.integer(data$Year) - 1L) +
        .02 * standard(data$effort) + standard(data$wind) * field[data$site])
      N <- rnbinom(nrow(data), mu = mu, size = 3)
      y <- replicate(3L, rbinom(length(N), N, .65))
      K <- max(120L, max(y) + 40L); stopifnot(K <= 400L)
      frame <- unmarked::unmarkedFramePCount(y, siteCovs = data)
      formula <- ~I(scale(temperature)^2) + Year ~ I(base::scale(habitat)^2) + Year +
        offset(.02 * base::scale(effort)) + (0 + base::scale(wind) || site)
      fits <- warnings <- before <- list(); calls <- 0L
      for (mixture in c("P", "NB")) {
        messages <- character(); calls <- calls + 1L; stopifnot(calls <= 2L)
        model <- withCallingHandlers(unmarked::pcount(formula, frame, K = K,
          mixture = mixture, engine = "TMB", threads = 1L), warning = function(w) {
            messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
          })
        fits[[mixture]] <- model; warnings[[mixture]] <- messages
        before[[mixture]] <- list(opt = model@opt, estimates = model@estimates,
          data = model@data, AIC = model@AIC, negLogLike = model@negLogLike, TMB_par = model@TMB$par)
      }
      list(fits = fits, warnings = warnings, before = before, calls = calls, K = K)
    }, error = function(e) {unmarked_zero_slope_cache$failure <- e; stop(e)})
  }
  unmarked_zero_slope_cache$value
}

test_that("zero-effect admission does not require unused learned random slopes", {
  f <- unmarked_zero_slope_fixture(); expect_identical(f$calls, 2L)
  rows <- unmarked_trained_fixture()$rows
  for (name in names(f$fits)) {
    model <- f$fits[[name]]; expect_identical(model@opt$convergence, 0L)
    expect_identical(f$warnings[[name]], character())
    expect_equal(.msr_parameter_count(model), if (name == "P") 9 else 10)
    expect_error(.msr_unmarked_prediction_function(model), "Learned random-effect")
    data <- model@data@siteCovs[rows, , drop = FALSE]
    data$wind <- data$site <- NULL
    expected <- unmarked_trained_reference(model, "state", data, re.form = NA)
    expect_true(all(is.finite(as.matrix(expected))))
    zero <- .msr_unmarked_prediction_function(model, re.form = NA)
    expect_equal(zero(model, "state", data), expected, tolerance = 1e-10)
    expect_error(zero(model, "state", data, re.form = NULL), "cannot supply conditional")
    object <- structure(list(opt_mod = model), class = "multiScaleR")
    point <- .msr_native_prediction_rows(object, data, "state", TRUE, list(re.f = NA))
    expect_equal(as.data.frame(point), expected, tolerance = 1e-10, ignore_attr = TRUE)
    expect_identical(attr(point, "prediction_contract")$effects$response, "zero_deviation_reference")
    snap <- prepare_prediction(object, data, component = "state", native_args = list(re.f = NA))
    expect_equal(predict(unserialize(serialize(snap, NULL))), expected, tolerance = 1e-10)
    expect_identical(snap$provenance$native_args$re.form, NA)
    expect_identical(model@opt, f$before[[name]]$opt)
    expect_identical(model@estimates, f$before[[name]]$estimates)
    expect_identical(model@data, f$before[[name]]$data)
    expect_identical(model@AIC, f$before[[name]]$AIC)
    expect_identical(model@negLogLike, f$before[[name]]$negLogLike)
    expect_identical(model@TMB$par, f$before[[name]]$TMB_par)
  }
})

test_that("ordinary unmarked curves share original native point targets", {
  inherited <- unmarked_trained_fixture(); extra <- unmarked_zero_slope_fixture()
  models <- c(inherited$fits, stats::setNames(extra$fits, c("slope_P", "slope_NB")))
  old <- options("contrasts"); on.exit(options(old), add = TRUE)
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  for (model in models) {
    options(old)
    data <- model@data@siteCovs
    defaults <- lapply(data, function(x) if (is.numeric(x)) mean(x) else x[1L])
    reference <- list()
    for (variable in c("habitat", "Year")) {
      x <- if (variable == "Year") data$Year[match(levels(data$Year), as.character(data$Year))] else
        seq(min(data$habitat), max(data$habitat), length.out = 7L)
      rows <- as.data.frame(lapply(defaults, rep, times = length(x)))
      rows[[variable]] <- x
      reference[[variable]] <- unmarked_trained_reference(model, "state", rows, re.form = NA)
    }
    object <- structure(list(opt_mod = model,
      scl_params = list(mean = c(habitat = 0), sd = c(habitat = 1))), class = "multiScaleR")
    before <- list(opt = model@opt, estimates = model@estimates, AIC = model@AIC, K = .msr_parameter_count(model))
    for (coding in list(old$contrasts, c("contr.helmert", "contr.poly"))) {
      options(contrasts = coding)
      plots <- plot_marginal_effects(object, type = "state", length.out = 7L)
      expect_setequal(names(plots), c("habitat", "Year"))
      contract <- attr(plots, "prediction_contract")
      expect_identical(contract$effects$response, "zero_deviation_reference")
      expect_identical(contract$integration$random_effects, "set_to_zero_not_integrated")
      expect_identical(contract$uncertainty$kernel_parameters, "fixed_at_estimate")
      expect_false(contract$metadata_establishes_eligibility)
      for (variable in names(reference)) {
        expect_equal(plots[[variable]]$data$fit, reference[[variable]]$Predicted, tolerance = 1e-8)
        expect_equal(plots[[variable]]$data$lwr, reference[[variable]]$lower, tolerance = 1e-8)
        expect_equal(plots[[variable]]$data$upr, reference[[variable]]$upper, tolerance = 1e-8)
      }
    }
    expect_identical(model@opt, before$opt)
    expect_identical(model@estimates, before$estimates)
    expect_identical(model@AIC, before$AIC)
    expect_equal(.msr_parameter_count(model), before$K)
  }
})

test_that("ordinary S4 control aliases match native semantics before target admission", {
  f <- unmarked_trained_fixture()
  for (model in f$fits) {
    rows <- model@data@siteCovs[f$rows, , drop = FALSE]
    expected <- unmarked_trained_reference(model, "state", rows, backTransform = FALSE, level = .9, re.form = NA)
    object <- structure(list(opt_mod = model), class = "multiScaleR")
    controls <- list(re.f = NA, backT = FALSE, lev = .9)
    canonical <- .msr_native_control_names(model, controls)
    expect_identical(names(canonical), c("re.form", "backTransform", "level"))
    expect_identical(unname(canonical), unname(controls))
    actual <- .msr_native_prediction_rows(object, rows, "state", TRUE, controls)
    expect_equal(as.data.frame(actual), expected, tolerance = 1e-10, ignore_attr = TRUE)
    snap <- prepare_prediction(object, rows, component = "state", native_args = controls)
    expect_equal(predict(snap), expected, tolerance = 1e-10)
    expect_identical(snap$prediction_contract$state$output_scale, "native_link")
    expect_error(.msr_native_control_names(model, list(n = FALSE), newdata = FALSE), "matches multiple")
  }
})

test_that("checked ordinary contexts survive two-worker prediction transfer", {
  f <- unmarked_trained_fixture(); model <- f$fits$P
  rows <- model@data@siteCovs[f$rows, , drop = FALSE]
  thin <- list(object = list(opt_mod = model), type = "state", se.fit = TRUE,
    native_args = list(re.form = NA), logical_names = character(),
    unmarked_predictor = .msr_unmarked_prediction_function(model, re.form = NA))
  expected <- .msr_native_raster_eval(thin, rows)
  workers <- parallel::makePSOCKcluster(2L, outfile = "")
  on.exit(parallel::stopCluster(workers), add = TRUE)
  cluster_prep(model, workers)
  result <- parallel::parLapply(workers, list(1:3, 4:7), function(i, thin, rows)
    multiScaleR:::.msr_native_raster_eval(thin, rows[i, , drop = FALSE]), thin, rows)
  expect_equal(as.matrix(do.call(rbind, result)), as.matrix(expected), tolerance = 1e-10, ignore_attr = TRUE)
  paths <- parallel::clusterEvalQ(workers, normalizePath(getNamespaceInfo("multiScaleR", "path"), winslash = "/"))
  expect_true(all(vapply(paths, identical, logical(1), normalizePath(getNamespaceInfo("multiScaleR", "path"), winslash = "/"))))
})

# Reuse the SAME eight failure-cached software fixture sites above. No new fit,
# sampler, recovery, profile or timing constructor is introduced here.
test_that("ordinary detection curves retain observation-trained native targets", {
  inherited <- unmarked_trained_fixture(); extra <- unmarked_zero_slope_fixture()
  models <- c(inherited$fits, stats::setNames(extra$fits, c("slope_P", "slope_NB")))
  old <- options("contrasts"); on.exit(options(old), add = TRUE)
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  for (model in models) {
    options(old)
    raw <- unmarked:::get_covariates(model, "det")
    defaults <- lapply(raw, function(x) if (is.numeric(x)) mean(x) else x[1L])
    reference <- list()
    for (variable in c("temperature", "Year")) {
      values <- if (variable == "Year") raw$Year[match(levels(raw$Year), as.character(raw$Year))] else
        seq(min(raw$temperature), max(raw$temperature), length.out = 7L)
      rows <- as.data.frame(lapply(defaults, rep, times = length(values)))
      rows[[variable]] <- values
      reference[[variable]] <- unmarked_trained_reference(model, "det", rows, re.form = NA)
    }
    object <- structure(list(opt_mod = model,
      scl_params = list(mean = c(habitat = 0), sd = c(habitat = 1))), class = "multiScaleR")
    before <- list(opt = model@opt, estimates = model@estimates, data = model@data,
      AIC = model@AIC, nll = model@negLogLike, K = .msr_parameter_count(model), par = model@TMB$par)
    for (coding in list(old$contrasts, c("contr.helmert", "contr.poly"))) {
      options(contrasts = coding)
      plots <- plot_marginal_effects(object, type = "det", length.out = 7L)
      expect_setequal(names(plots), c("temperature", "Year"))
      contract <- attr(plots, "prediction_contract")
      # One related contract comparison keeps every distinction explicit without
      # multiplying assertion instances for identical metadata fields.
      expect_identical(list(component = contract$component,
        response = contract$effects$response, integration = contract$integration$random_effects,
        kernel = contract$uncertainty$kernel_parameters, intervals = contract$uncertainty$intervals,
        observation_noise = contract$uncertainty$observation_noise,
        eligible = contract$metadata_establishes_eligibility),
        list(component = "det", response = "zero_deviation_reference",
          integration = "set_to_zero_not_integrated", kernel = "fixed_at_estimate",
          intervals = "native_conditional_confidence", observation_noise = FALSE, eligible = FALSE))
      for (variable in names(reference)) {
        expect_equal(unname(as.matrix(plots[[variable]]$data[c("fit", "lwr", "upr")])),
          unname(as.matrix(reference[[variable]][c("Predicted", "lower", "upper")])), tolerance = 1e-8)
      }
    }
    expect_identical(list(opt = model@opt, estimates = model@estimates, data = model@data,
      AIC = model@AIC, nll = model@negLogLike, K = .msr_parameter_count(model), par = model@TMB$par), before)
  }
})

test_that("ordinary lambda plotting aliases and zero-only curve metadata agree", {
  inherited <- unmarked_trained_fixture(); extra <- unmarked_zero_slope_fixture()
  models <- c(inherited$fits, stats::setNames(extra$fits, c("slope_P", "slope_NB")))
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  for (model in models) {
    object <- structure(list(opt_mod = model,
      scl_params = list(mean = c(habitat = 0), sd = c(habitat = 1))), class = "multiScaleR")
    state <- plot_marginal_effects(object, type = "state", length.out = 7L)
    alias <- plot_marginal_effects(object, type = "lambda", length.out = 7L)
    expect_identical(names(alias), names(state))
    expect_identical(attr(alias, "prediction_contract"), attr(state, "prediction_contract"))
    for (name in names(state)) expect_equal(alias[[name]]$data, state[[name]]$data, tolerance = 1e-10)
    if (unmarked:::has_random(model@formlist$state)) {
      reduced <- model
      reduced@data@siteCovs$site <- reduced@data@siteCovs$wind <- NULL
      clone <- object; clone$opt_mod <- reduced
      # Required conditional metadata stays guarded; zero curves need no unused Z.
      expect_error(.msr_unmarked_prediction_function(reduced), "Learned random|Original unmarked covariates")
      for (type in c("state", "det")) {
        expected <- if (type == "state") state else plot_marginal_effects(object, type = type, length.out = 7L)
        actual <- plot_marginal_effects(clone, type = type, length.out = 7L)
        expect_identical(names(actual), names(expected))
        for (name in names(expected)) expect_equal(actual[[name]]$data, expected[[name]]$data, tolerance = 1e-10)
      }
      expect_identical(reduced@opt, model@opt)
      expect_identical(reduced@estimates, model@estimates)
      expect_identical(reduced@TMB$par, model@TMB$par)
      expect_identical(reduced@negLogLike, model@negLogLike)
      expect_identical(.msr_parameter_count(reduced), .msr_parameter_count(model))
    }
  }
  bad <- inherited$fits$P; bad@data@siteCovs$temperature <- NULL
  object <- structure(list(opt_mod = bad,
    scl_params = list(mean = c(habitat = 0), sd = c(habitat = 1))), class = "multiScaleR")
  expect_error(plot_marginal_effects(object, type = "det"), "observation covariates.*temperature")
})

test_that("unqualified native subclasses keep their original metadata", {
  raw <- data.frame(y = 1:12, x = sin(1:12), unused = cos(1:12))
  frame <- model.frame(y ~ x + unused - unused, raw)
  # These terms-only prototypes deliberately have an extra provider component.
  # They have no fitted likelihood, coefficients, covariance or inference target.
  for (kind in list(c("custom_glm", "glm", "lm"), c("custom_cox", "coxph"),
      c("gam", "glm", "lm"), c("coxph.penal", "coxph"))) {
    model <- structure(list(terms = terms(frame), model = frame,
      provider_components = list(dispersion = ~unused)), class = kind)
    expect_identical(.msr_native_active_metadata(model), model)
  }
})


# Reuse the same eight deterministic software fits; no additional fit sites.
test_that("requested unmarked point contexts ignore unrelated component metadata", {
  f <- unmarked_trained_fixture()
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  for (name in names(f$fits)) for (type in c("state", "det")) {
    model <- f$fits[[name]]
    rows <- model@data@siteCovs[f$rows, , drop = FALSE]
    expected <- unmarked_trained_reference(model, type, rows)
    unrelated <- if (type == "state") "temperature" else c("habitat", "effort")
    reduced <- model
    for (variable in unrelated) reduced@data@siteCovs[[variable]] <- NULL
    keep <- if (type == "state") c("habitat", "effort", "Year",
      if (unmarked:::has_random(model@formlist$state)) "site") else c("temperature", "Year")
    rows <- rows[keep]
    before <- list(opt = reduced@opt, estimates = reduced@estimates,
      nll = reduced@negLogLike, AIC = reduced@AIC, K = .msr_parameter_count(reduced),
      data = reduced@data, TMB_par = reduced@TMB$par)
    scoped <- .msr_unmarked_prediction_function(reduced, component = type)
    expect_equal(scoped(reduced, type, rows), expected, tolerance = 1e-10)
    object <- structure(list(opt_mod = reduced), class = "multiScaleR")
    point <- .msr_native_prediction_rows(object, rows, type, TRUE)
    expect_equal(as.data.frame(point), expected, tolerance = 1e-10, ignore_attr = TRUE)
    contract <- attr(point, "prediction_contract")
    expect_identical(list(component = contract$component,
      intervals = contract$uncertainty$intervals, noise = contract$uncertainty$observation_noise,
      kernel = contract$uncertainty$kernel_parameters, eligible = contract$metadata_establishes_eligibility),
      list(component = type, intervals = "native_conditional_confidence", noise = FALSE,
        kernel = "fixed_at_estimate", eligible = FALSE))
    snap <- prepare_prediction(object, rows, component = type)
    expect_equal(predict(unserialize(serialize(snap, NULL))), expected, tolerance = 1e-10)
    thin <- list(object = list(opt_mod = reduced), type = type, se.fit = TRUE,
      native_args = list(), logical_names = character(), unmarked_predictor = scoped)
    chunks <- rbind(.msr_native_raster_eval(thin, rows[1:3, , drop = FALSE]),
      .msr_native_raster_eval(thin, rows[4:7, , drop = FALSE]))
    expect_equal(as.matrix(chunks), as.matrix(expected), tolerance = 1e-10, ignore_attr = TRUE)
    curve_object <- structure(list(opt_mod = model,
      scl_params = list(mean = c(habitat = 0), sd = c(habitat = 1))), class = "multiScaleR")
    expected_curves <- plot_marginal_effects(curve_object, type = type, length.out = 5L)
    curve_object$opt_mod <- reduced
    curves <- plot_marginal_effects(curve_object, type = type, length.out = 5L)
    expect_identical(names(curves), names(expected_curves))
    expect_identical(attr(curves, "prediction_contract"), attr(expected_curves, "prediction_contract"))
    for (variable in names(curves))
      expect_equal(curves[[variable]]$data, expected_curves[[variable]]$data, tolerance = 1e-10)
    # Native frame removal still requires both components, even for one output.
    expect_error(.msr_unmarked_prediction_function(reduced), "Original unmarked covariates")
    expect_error(prepare_prediction(object, model@data, component = type), "Original unmarked covariates")
    expect_error(scoped(reduced, type, model@data), "explicit data-frame")
    expect_error(scoped(reduced, if (type == "state") "det" else "state", rows), "prepared type")
    expect_identical(list(opt = reduced@opt, estimates = reduced@estimates,
      nll = reduced@negLogLike, AIC = reduced@AIC, K = .msr_parameter_count(reduced),
      data = reduced@data, TMB_par = reduced@TMB$par), before)
  }
})

test_that("unmarked requested components keep their own guards and native ZIP target", {
  f <- unmarked_trained_fixture()
  for (name in names(f$fits)) {
    model <- f$fits[[name]]
    rows <- model@data@siteCovs[f$rows, , drop = FALSE]
    for (type in c("state", "det")) {
      required <- if (type == "state") "habitat" else "temperature"
      changed <- model
      changed@data@siteCovs[[required]][1L] <- changed@data@siteCovs[[required]][1L] + .1
      expect_error(.msr_unmarked_prediction_function(changed, component = type), "do not match")
      changed@data@siteCovs[[required]] <- NULL
      expect_error(.msr_unmarked_prediction_function(changed, component = type), "Original unmarked covariates")
      scoped <- .msr_unmarked_prediction_function(model, component = type)
      absent <- rows; absent[[required]] <- NULL
      expect_error(scoped(model, type, absent), "explicitly contain")
    }
    expect_error(.msr_unmarked_prediction_function(model, component = c("state", "det")), "one native")
    if (name == "ZIP") {
      # Retain both the native overall mean and its psi variance contribution.
      scoped <- .msr_unmarked_prediction_function(model, component = "state")
      for (response in c(FALSE, TRUE)) for (level in c(.8, .95))
        expect_equal(scoped(model, "state", rows, backTransform = response, level = level),
          unmarked_trained_reference(model, "state", rows, backTransform = response, level = level),
          tolerance = 1e-10)
      nuisance <- .msr_unmarked_prediction_function(model, component = "psi")
      expect_equal(nuisance(model, "state", rows),
        unmarked_trained_reference(model, "state", rows), tolerance = 1e-10)
      # Nuisance requests retain the original full-context admission boundary.
      reduced <- model; reduced@data@siteCovs$temperature <- NULL
      expect_error(.msr_unmarked_prediction_function(reduced, component = "psi"), "Original unmarked covariates")
      reduced <- model
      reduced@formlist$state <- ~I(base::scale(habitat)^2) + Year + (scale(effort) | site)
      expect_error(.msr_unmarked_prediction_function(reduced, component = "psi"), "Learned random-effect")
    }
  }
})

test_that("unrelated learned random slopes do not block conditional detection", {
  f <- unmarked_zero_slope_fixture()
  rows <- unmarked_trained_fixture()$rows
  for (model in f$fits) {
    data <- model@data@siteCovs[rows, , drop = FALSE]
    expected <- unmarked_trained_reference(model, "det", data)
    reduced <- model
    for (variable in c("wind", "site", "habitat", "effort")) reduced@data@siteCovs[[variable]] <- NULL
    data <- data[c("temperature", "Year")]
    scoped <- .msr_unmarked_prediction_function(reduced, component = "det")
    expect_equal(scoped(reduced, "det", data), expected, tolerance = 1e-10)
    object <- structure(list(opt_mod = reduced), class = "multiScaleR")
    point <- .msr_native_prediction_rows(object, data, "det", TRUE)
    expect_equal(as.data.frame(point), expected, tolerance = 1e-10, ignore_attr = TRUE)
    expect_identical(attr(point, "prediction_contract")$integration$random_effects, "native_method")
    expect_error(.msr_unmarked_prediction_function(reduced, component = "state"), "Learned random-effect")
    expect_error(.msr_unmarked_prediction_function(reduced), "Learned random-effect")
    expect_error(prepare_prediction(object, model@data, component = "det"), "Learned random-effect")
    zero <- .msr_unmarked_prediction_function(reduced, re.form = NA, component = "det")
    expect_error(zero(reduced, "det", data, re.form = NULL), "cannot supply conditional")
    expect_identical(list(opt = reduced@opt, estimates = reduced@estimates,
      nll = reduced@negLogLike, AIC = reduced@AIC, K = .msr_parameter_count(reduced), par = reduced@TMB$par),
      list(opt = model@opt, estimates = model@estimates, nll = model@negLogLike,
        AIC = model@AIC, K = .msr_parameter_count(model), par = model@TMB$par))
  }
})

test_that("one requested unmarked basis is reused across prediction chunks", {
  f <- unmarked_trained_fixture(); model <- f$fits$P
  rows <- model@data@siteCovs[f$rows, , drop = FALSE]
  original <- .msr_unmarked_prediction_basis; calls <- character()
  local_mocked_bindings(.msr_unmarked_prediction_basis = function(object, type, verify_random = TRUE) {
    calls <<- c(calls, type)
    original(object, type, verify_random)
  }, .package = "multiScaleR")
  scoped <- .msr_unmarked_prediction_function(model, component = "det")
  expected <- unmarked_trained_reference(model, "det", rows)
  actual <- rbind(scoped(model, "det", rows[1:3, , drop = FALSE]),
    scoped(model, "det", rows[4:7, , drop = FALSE]))
  rownames(actual) <- NULL
  expect_equal(actual, expected, tolerance = 1e-10)
  expect_identical(calls, "det")
})
