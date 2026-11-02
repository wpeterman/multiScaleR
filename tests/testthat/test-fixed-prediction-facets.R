# Deterministic adapter conformance on established software fixtures. These
# tests assess numerical contracts, not scale recovery or interval calibration.
test_that("prepared standard splines retain native knots with factors and logic", {
  skip_if_not_installed("splines")
  dat <- data.frame(y = seq_len(12), z = c(-2, -1.7, -1.2, -.8, -.4, 0, .3, .7, 1, 1.5, 2.2, 3),
    group = factor(rep(c("a", "b"), 6)), site = factor(rep(seq_len(4), each = 3)),
    a = rep(c(TRUE, FALSE, FALSE), 4), b = rep(c(FALSE, TRUE), 6))
  nd <- dat[c(10, 2, 6, 10), , drop = FALSE]
  for (kind in c("ns", "bs")) {
    make <- function(qualified) {
      rhs <- paste0(if (qualified) "splines::" else "", kind,
        "(z,df=3) * group + I(a | b) + (1 | site)")
      form <- as.formula(paste("y ~", rhs), env = asNamespace("splines"))
      .msr_bayes_formula(form, dat, response = TRUE, active_frame = TRUE)
    }
    b <- make(FALSE); q <- make(TRUE)
    expect_equal(model.matrix(q$terms, q$frame), model.matrix(b$terms, b$frame),
      tolerance = 1e-12, ignore_attr = TRUE)
    qt <- delete.response(q$terms); bt <- delete.response(b$terms)
    X <- model.matrix(qt, model.frame(qt, nd))
    expected <- model.matrix(bt, model.frame(bt, nd))
    expect_equal(X, expected, tolerance = 1e-12, ignore_attr = TRUE)
    # The independent native basis predictor uses its saved training attributes.
    oracle <- predict(q$frame[[2L]], nd$z)
    expect_equal(as.numeric(X[, 1L + seq_len(ncol(oracle)), drop = FALSE]),
      as.numeric(oracle), tolerance = 1e-12)
    expect_identical(q$random[[1L]]$group, "site")
    expect_identical(.msr_fixed_prediction_terms(q$terms, q$frame), q$terms)
  }
})

test_that("qualified fixed normalization agrees across public SQ prediction facets", {
  for (kernel in c("gaussian", "exp")) {
    f <- refit_spatial_fixture(kernel = kernel)
    dat <- f$data
    set.seed(442) # Preserve the established raster-formula response fixture.
    dat$y <- rpois(nrow(dat), exp(3 + 1.5 * dat$rough + .1 * dat$z))
    model <- glm(y ~ rough * group + base::scale(z) + offset(log(exposure)),
      dat, family = poisson())
    retained <- character()
    object <- withCallingHandlers(multiScale_optim(model, f$inputs,
      backend = "profile_glm", n_cores = 1, verbose = FALSE), warning = function(w) {
        text <- conditionMessage(w)
        if (grepl("Joint Wald uncertainty unavailable", text, fixed = TRUE)) {
          retained <<- c(retained, text); invokeRestart("muffleWarning")
        }
      })
    if (length(retained)) message(paste(retained, collapse = "\n"))
    nd <- dat; nd$group <- factor("b", levels = levels(dat$group))
    nd$z <- .2; nd$exposure <- 2
    design <- .msr_profile_prediction_design(object, nd, f$inputs)
    manual <- design$X
    column <- match("base::scale(z)", colnames(manual))
    expect_false(is.na(column))
    manual[, column] <- (nd$z - mean(dat$z)) / sd(dat$z)
    expect_equal(design$X, manual, tolerance = 1e-12)
    beta <- coef(object$opt_mod)[colnames(manual)]
    # Prediction must use its requested numerical observed-information covariance.
    # Native IRLS and closed Poisson information are independent accuracy checks,
    # using the existing 1e-6 SE contract, not bitwise covariance identity.
    V <- vcov(object, type = "conditional")
    train <- model.matrix(object$opt_mod)[, colnames(manual), drop = FALSE]
    offset <- model.offset(object$opt_mod$model)
    index <- object$uncertainty$beta_index
    target_beta <- object$uncertainty$par[index]
    weight <- exp(drop(train %*% target_beta) + offset)
    information <- crossprod(train, train * weight)
    expect_equal(object$uncertainty$hessian[index, index, drop = FALSE],
      information, tolerance = 1e-6, ignore_attr = TRUE)
    native_V <- vcov(object$opt_mod)[colnames(manual), colnames(manual), drop = FALSE]
    eta <- drop(manual %*% beta) + log(nd$exposure)
    se <- sqrt(pmax(0, rowSums((manual %*% V) * manual)))
    link <- predict(object, nd, f$inputs, type = "link", se.fit = TRUE,
      uncertainty = "conditional")
    point <- predict(object, nd, f$inputs, se.fit = TRUE, uncertainty = "conditional")
    expect_equal(as.numeric(link$fit), as.numeric(eta), tolerance = 1e-10)
    expect_equal(as.numeric(link$se.fit), as.numeric(se), tolerance = 1e-10)
    native_se <- sqrt(pmax(0, rowSums((manual %*% native_V) * manual)))
    expect_equal(as.numeric(link$se.fit), as.numeric(native_se), tolerance = 1e-6)
    message("FIXED_BASIS_INFO ", kernel, " native_SE_relative_gap=",
      format(max(abs(link$se.fit / native_se - 1)), digits = 17))
    expect_equal(as.numeric(point$fit), as.numeric(exp(eta)), tolerance = 1e-10)
    expect_equal(as.numeric(point$se.fit), as.numeric(exp(eta) * se), tolerance = 1e-10)
    snapshot <- prepare_prediction(object, nd, f$inputs)
    path <- tempfile(fileext = ".rds"); saveRDS(snapshot, path)
    restored <- readRDS(path); unlink(path)
    expect_equal(predict(restored, se.fit = TRUE, uncertainty = "conditional"), point)
    expect_equal(predict(restored, type = "link", se.fit = TRUE,
      uncertainty = "conditional"), link)
    inputs <- kernel_scale.raster(f$rasters, multiScaleR = object,
      scale_center = TRUE, match_buffer = TRUE, site_values = list(group = "b", z = .2),
      offset_values = c(exposure = 2), verbose = FALSE)
    map <- predict(object, inputs, se.fit = TRUE, uncertainty = "conditional")
    actual <- terra::extract(map, f$points, ID = FALSE)
    expect_equal(actual$fit, as.numeric(point$fit), tolerance = 1e-8)
    expect_equal(actual$se.fit, as.numeric(point$se.fit), tolerance = 1e-8)
    expect_error(predict(object, inputs, se.fit = TRUE, uncertainty = "joint"), "conditional SEs")
    expect_equal(.msr_selection_parameter_count(object),
      as.numeric(attr(logLik(object$opt_mod), "df")) + 1L)
  }
})

test_that("conditional raw and retained routes share frozen public predictions", {
  skip_if_not_installed("survival")
  coxph <- survival::coxph
  f <- profile_glm_fixture(); dat <- f$data
  dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
  dat$choice <- rep(seq_len(nrow(dat) / 3L), each = 3L)
  env <- new.env(parent = asNamespace("survival")); env$dat <- dat
  env$form <- as.formula("case_ ~ hab + group + base::scale(z) + offset(log(exposure)) + strata(choice)", env = env)
  model <- evalq(survival::clogit(form, data = dat, model = TRUE, x = TRUE, y = TRUE), env)
  raw <- .msr_clogit_contract(model, f$inputs)
  expect_identical(raw$native_data_mode, "raw_variables")
  saved <- model; saved$call$data <- quote(unavailable_original_data)
  # Expire the raw call environment; the evaluated frame and prediction calls remain.
  environment(saved$formula) <- asNamespace("survival")
  attr(saved$terms, ".Environment") <- asNamespace("survival")
  retained <- .msr_clogit_contract(saved, f$inputs)
  expect_identical(retained$native_data_mode, "retained_model_frame")
  make <- function(contract) {
    native <- .msr_clogit_native(contract, log(.3))
    contract$objective <- list(engine = contract$engine, X_fixed = contract$X_fixed,
      offset = contract$offset, response_design = contract$response_design)
    structure(list(backend = "conditional_logit", opt_mod = native,
      kernel_inputs = f$inputs, profile_fit = list(contract = contract, par = log(.3)),
      diagnostics = .empty_diagnostics(), uncertainty = list(available = FALSE)),
      class = "multiScaleR")
  }
  a <- make(raw); b <- make(retained)
  nd <- dat[rev(seq_len(12)), , drop = FALSE]; nd$z <- .2
  # A separate bare-call native fit supplies the provider's own centering and
  # SE convention. Compare like native outputs, not native risk SEs to our delta SEs.
  native_dat <- dat
  native_dat$hab <- .msr_eval_kernel_engine(raw$engine, log(.3))$value[, 1L]
  bare_form <- as.formula("case_ ~ hab + group + scale(z) + offset(log(exposure)) + strata(choice)",
    env = asNamespace("survival"))
  bare_native <- survival::clogit(bare_form, data = native_dat,
    model = TRUE, x = TRUE, y = TRUE)
  native_nd <- nd; native_nd$hab <- native_dat$hab[rev(seq_len(12))]
  oracle <- predict(bare_native, native_nd, type = "lp", se.fit = TRUE, reference = "zero")
  for (fit in list(a$opt_mod, b$opt_mod)) {
    direct <- predict(fit, native_nd, type = "lp", se.fit = TRUE, reference = "zero")
    expect_equal(as.numeric(direct$fit), as.numeric(oracle$fit), tolerance = 1e-10)
    expect_equal(as.numeric(direct$se.fit), as.numeric(oracle$se.fit), tolerance = 1e-10)
  }
  design <- .msr_profile_prediction_design(b, nd, f$inputs)
  X <- design$X; X[, "base::scale(z)"] <- (nd$z - mean(dat$z)) / sd(dat$z)
  beta <- coef(b$opt_mod)[colnames(X)]; V <- vcov(b$opt_mod)[colnames(X), colnames(X), drop = FALSE]
  eta <- drop(X %*% beta) + log(nd$exposure)
  se <- sqrt(pmax(0, rowSums((X %*% V) * X)))
  for (type in c("link", "risk", "response")) {
    expected <- predict(a, nd, f$inputs, type = type, se.fit = TRUE, uncertainty = "conditional")
    actual <- predict(b, nd, f$inputs, type = type, se.fit = TRUE, uncertainty = "conditional")
    expect_equal(actual, expected, tolerance = 1e-10)
    prepared <- prepare_prediction(b, nd, f$inputs)
    expect_equal(predict(prepared, type = type, se.fit = TRUE,
      uncertainty = "conditional"), actual, tolerance = 1e-10)
  }
  link <- predict(b, nd, f$inputs, type = "link", se.fit = TRUE, uncertainty = "conditional")
  risk <- predict(b, nd, f$inputs, type = "risk", se.fit = TRUE, uncertainty = "conditional")
  expect_equal(as.numeric(link$fit), as.numeric(eta), tolerance = 1e-10)
  expect_equal(as.numeric(link$se.fit), as.numeric(se), tolerance = 1e-10)
  expect_equal(as.numeric(risk$fit), as.numeric(exp(eta)), tolerance = 1e-10)
  expect_equal(as.numeric(risk$se.fit), as.numeric(exp(eta) * se), tolerance = 1e-10)
  no_sets <- nd; no_sets$choice <- NULL
  prepared <- prepare_prediction(b, no_sets, f$inputs)
  expect_equal(as.numeric(predict(prepared, type = "link")), as.numeric(eta), tolerance = 1e-10)
  expect_error(predict(prepared, type = "response"), "choice")
  expect_error(predict(prepared, type = "link", se.fit = TRUE, uncertainty = "joint"), "conditional SEs")
})
