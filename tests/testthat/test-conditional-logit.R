# Reuse the existing binned numerical fixture; no new simulation panel.
clogit_fixture <- function(kernel = "gaussian") {
  coxph <- survival::coxph
  f <- profile_glm_fixture(kernel = kernel)
  dat <- f$data
  dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
  dat$choice <- rep(seq_len(nrow(dat) / 3L), each = 3L)
  form <- stats::as.formula("case_ ~ hab + z + group + offset(log(exposure)) + strata(choice)",
                           env = asNamespace("survival"))
  f$model <- survival::clogit(form, data = dat, model = TRUE, x = TRUE, y = TRUE)
  f$data <- dat
  f
}

test_that("conditional likelihood and derivatives match the native one-case likelihood", {
  skip_if_not_installed("survival")
  for (kernel in c("gaussian", "exp")) {
    f <- clogit_fixture(kernel)
    c <- .msr_clogit_contract(f$model, f$inputs)
    o <- .msr_clogit_objective(c)
    p <- c(log(.3), rep(.15, o$nb))
    numerical <- vapply(seq_along(p), function(i) {
      a <- b <- p; a[i] <- a[i] + 1e-5; b[i] <- b[i] - 1e-5
      (o$fn(a) - o$fn(b)) / 2e-5
    }, 0)
    expect_equal(o$gr(p), numerical, tolerance = 1e-6, ignore_attr = TRUE)
    n <- .msr_clogit_native(c, log(.3))
    beta <- coef(n)[c(c$fixed_columns, c$covariates)]
    expect_equal(o$fn(c(log(.3), beta)), -as.numeric(logLik(n)), tolerance = 1e-10)
    expect_lt(max(abs(tail(o$gr(c(log(.3), beta)), length(beta)))), 1e-6)
  }
})

test_that("conditional contracts reject broken choice sets and unsupported data", {
  skip_if_not_installed("survival")
  expect_error(.msr_clogit_groups(c(1, 1, 2), c(1, 0, 1)), "two alternatives")
  expect_error(.msr_clogit_groups(c(1, 1), c(1, 1)), "exactly one")
  expect_error(.msr_clogit_groups(c(1, NA), c(1, 0)), "missing")
  f <- clogit_fixture()
  m <- f$model; m$call$weights <- quote(w)
  expect_error(.msr_clogit_contract(m, f$inputs), "weights")
  m <- f$model; m$na.action <- 1L
  expect_error(.msr_clogit_contract(m, f$inputs), "Missing-row")
  m <- f$model; m$model[[1L]][, 1L] <- 2
  expect_error(.msr_clogit_contract(m, f$inputs), "constant-time")
  b <- f$inputs; b$binned$point_ids[1] <- "unmatched"
  expect_error(.msr_clogit_contract(f$model, b), "row names")
  m <- f$model; m$coefficients[1] <- NA
  expect_error(.msr_clogit_contract(m, f$inputs), "nonaliased")
  m <- f$model; class(m) <- c("coxph.penal", class(m))
  expect_error(.msr_clogit_contract(m, f$inputs), "special Cox")
  z <- .msr_clogit_softmax(c(1000, 1001, -1000, -1001), list(1:2, 3:4))
  expect_true(all(is.finite(z$log_probability)))
  expect_equal(sum(z$probability), 2)
})

test_that("joint conditional fitting keeps starts and native prediction contracts", {
  skip_if_not_installed("survival")
  f <- clogit_fixture()
  starts <- matrix(c(.08, .25, .6), ncol = 1)
  fit <- multiScale_optim(f$model, f$inputs, backend = "conditional_logit",
    n_cores = 1, verbose = FALSE, profile_control = list(starts = starts))
  expect_equal(exp(fit$profile_fit$starts), starts, ignore_attr = TRUE)
  expect_equal(nrow(fit$profile_fit$start_results), 3L)
  expect_true(any(fit$profile_fit$start_results$eligible))
  records <- fit$profile_fit$start_results
  expect_true(all(records$projected_gradient <= records$joint_projected_gradient))
  refined <- records$terminal_coefficients_refined
  expect_equal(records$projected_gradient[refined], records$cold_projected_gradient[refined])
  expect_true(all(records$native_gap[refined] <= 1e-5 &
                  records$coefficient_gap[refined] <= 1e-4 &
                  records$covariance_gap[refined] <= 1e-4))
  expect_true(all(is.na(fit$scale_est$SE)))
  expect_false(fit$uncertainty$available)
  expect_equal(fit$optim_results$value, -as.numeric(logLik(fit$opt_mod)), tolerance = 1e-5)
  p <- predict(fit)
  expect_equal(ncol(fit$profile_fit$contract$frame), 0L)
  expect_false(any(c("model", "data", "engine", "X_fixed", "y") %in%
                     names(fit$profile_fit$contract)))
  expect_identical(fit$fitted_mod_original, f$model)
  legacy <- fit
  legacy$profile_fit$contract <- .msr_clogit_contract(f$model, f$inputs)
  legacy$profile_fit$contract$objective <- fit$profile_fit$contract$objective
  expect_identical(predict(legacy), p)
  expect_identical(predict(legacy, se.fit = TRUE, uncertainty = "conditional"),
                   predict(fit, se.fit = TRUE, uncertainty = "conditional"))
  expect_lt(length(serialize(fit, NULL)), length(serialize(legacy, NULL)))
  expect_equal(as.numeric(rowsum(p, f$data$choice)), rep(1, 80), tolerance = 1e-12)
  risk <- predict(fit$opt_mod, type = "risk", reference = "zero")
  reference <- risk / as.numeric(rowsum(risk, f$data$choice))[f$data$choice]
  expect_equal(unname(p), unname(reference), tolerance = 1e-10)
  expect_error(predict(fit, se.fit = TRUE), "conditional SEs")
  expect_true(all(is.finite(predict(fit, se.fit = TRUE, uncertainty = "conditional")$se.fit)))
  expect_equal(predict(fit, newdata = f$data, kernel_inputs = f$inputs), p, tolerance = 1e-12)
  nd <- f$data; nd$choice <- NULL
  expect_error(predict(fit, newdata = nd, kernel_inputs = f$inputs), "choice-set IDs")
  expect_equal(predict(unserialize(serialize(fit, NULL))), p)
  expect_error(update(fit$opt_mod), "explicit data")
  nd <- f$data
  nd$hab <- .msr_eval_kernel_engine(fit$profile_fit$contract$objective$engine, fit$profile_fit$par)$value[,1]
  expect_equal(as.numeric(logLik(update(fit$opt_mod, data=nd))), as.numeric(logLik(fit$opt_mod)), tolerance=1e-10)
  nd <- f$data; nd$group <- "new_level"
  expect_error(predict(fit, newdata = nd, kernel_inputs = f$inputs), "new level")
  # New data retain fitted transforms and normalization under row reordering.
  order <- rev(seq_len(nrow(f$data)))
  expect_equal(unname(predict(fit, newdata = f$data[order, ], kernel_inputs = f$inputs)),
               unname(p[order]), tolerance = 1e-12)
  # A wrapper carries the same statistical model and preserves movement metadata.
  wrapped <- structure(list(model = f$model, sl_ = "retained", ta_ = "retained"), class = "fit_clogit")
  w <- multiScale_optim(wrapped, f$inputs, backend = "conditional_logit", n_cores = 1,
    verbose = FALSE, profile_control = list(starts = starts))
  expect_identical(w$opt_mod$sl_, wrapped$sl_)
  expect_identical(w$opt_mod$ta_, wrapped$ta_)
  expect_equal(predict(w), p)
  s <- summary(w)
  expect_s3_class(s, "summary_multiScaleR")
  expect_true(all(is.finite(s$opt_dist$Mean)))
  expect_true(all(is.na(s$opt_dist[, c("2.5%", "97.5%")])))
  expect_equal(kernel_dist(fit)$Mean, s$opt_dist$Mean)
})

test_that("simultaneous scales have correct conditional derivatives", {
  skip_if_not_installed("survival")
  f <- clogit_fixture()
  b <- f$inputs$binned; j <- rev(seq_along(b$bin_rep))
  b$vsum$second <- b$vsum$hab[, j]; b$csum$second <- b$csum$hab[, j]
  f$inputs$binned <- b
  dat <- f$data
  dat$second <- .msr_eval_kernel_engine(.msr_prepare_kernel_engine(b), log(c(.2,.4)))$value[,2]
  coxph <- survival::coxph; strata <- survival::strata
  m <- survival::clogit(case_ ~ hab + second + z + strata(choice), data=dat,model=TRUE)
  con <- .msr_clogit_contract(m,f$inputs); obj <- .msr_clogit_objective(con)
  p <- c(log(c(.2,.4)),rep(.1,obj$nb))
  numerical <- vapply(seq_along(p),function(i) {
    a <- b <- p; a[i] <- a[i]+1e-5; b[i] <- b[i]-1e-5
    (obj$fn(a)-obj$fn(b))/2e-5
  },0)
  expect_equal(obj$gr(p),numerical,tolerance=1e-6,ignore_attr=TRUE)
  fit <- multiScale_optim(m,f$inputs,backend='conditional_logit',n_cores=1,verbose=FALSE,
    profile_control=list(starts=rbind(c(.1,.6),c(.6,.1))))
  expect_equal(nrow(fit$scale_est),2L)
  expect_equal(nrow(fit$profile_fit$start_results),2L)
})

test_that("failed native audits remain explicit for every independent start", {
  skip_if_not_installed("survival")
  f <- clogit_fixture()
  local_mocked_bindings(.msr_clogit_native = function(...) stop("forced audit failure"))
  e <- tryCatch(multiScale_optim(f$model,f$inputs,backend='conditional_logit',n_cores=1,
    verbose=FALSE,profile_control=list(starts=matrix(c(.1,.6),ncol=1))),error=identity)
  expect_s3_class(e,'msr_clogit_fit_error')
  expect_equal(nrow(e$start_results),2L)
  expect_false(any(e$start_results$eligible))
  expect_true(all(grepl('forced audit failure',e$start_results$message)))
})

test_that("terminal refinement cannot conceal coefficient disagreement", {
  skip_if_not_installed("survival")
  f <- clogit_fixture()
  native <- .msr_clogit_native
  local_mocked_bindings(.msr_clogit_native = function(...) {
    z <- native(...)
    z$coefficients <- z$coefficients + .01
    z
  })
  e <- tryCatch(multiScale_optim(f$model, f$inputs, backend = "conditional_logit",
    n_cores = 1, verbose = FALSE,
    profile_control = list(starts = matrix(c(.1, .6), ncol = 1))), error = identity)
  expect_s3_class(e, "msr_clogit_fit_error")
  expect_equal(nrow(e$start_results), 2L)
  expect_false(any(e$start_results$eligible))
  expect_false(any(e$start_results$terminal_coefficients_refined))
  expect_true(all(e$start_results$coefficient_gap > 1e-4))
})
