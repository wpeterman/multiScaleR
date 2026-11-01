test_that("native formula products retain contrasts, nesting and literal names", {
  n <- 18L
  dat <- data.frame(y = sin(seq_len(n)), h = seq(-1, 1, length.out = n),
    k = cos(seq_len(n) * .7), z = sin(seq_len(n) * .3),
    f = factor(rep(c("a", "wood land", "x:y"), 6)),
    g = ordered(rep(c("low", "high"), 9)), effort = seq(1, 2, length.out = n), check.names = FALSE)
  contrasts(dat$f) <- contr.sum(3)
  dat[["habitat distance"]] <- dat$h; dat[["land class"]] <- dat$f
  forms <- c("y ~ h+k+z+f", "y ~ h*f", "y ~ h:f", "y ~ 0+h:f", "y ~ h*k*f",
    "y ~ (h+k+z+f)^2", "y ~ h/f", "y ~ h %in% f", "y ~ 0+h*k", "y ~ h*f-f",
    "y ~ h*k-h-k-1", "y ~ h+poly(z,2)+f+offset(log(effort))",
    "y ~ `habitat distance` * `land class`", "y ~ h+k+I(z^2)",
    "y ~ h*f+k %in% f+z:g", "y ~ h^2", "y ~ h+k-k", "y ~ -1+h*f",
    "y ~ . - `habitat distance` - `land class` - effort - g - f - z")
  for (text in forms) {
    mf <- model.frame(as.formula(text), dat); tt <- terms(mf)
    covs <- intersect(c("h", "k", "habitat distance"), .msr_active_formula_variables(tt))
    contrasts <- attr(model.matrix(tt, mf), "contrasts")
    plan <- .msr_compile_response_design(tt, mf, covs, contrasts = contrasts)
    K <- as.matrix(dat[covs]); K[1, ] <- 0
    nd <- dat; for (j in seq_along(covs)) nd[[covs[j]]] <- K[, j]
    expected <- model.matrix(tt, model.frame(tt, nd), contrasts.arg = contrasts)
    X <- .msr_response_design_matrix(plan, K)
    expect_identical(colnames(X), colnames(expected), info = text)
    expect_equal(X, expected, tolerance = 1e-12, ignore_attr = TRUE, info = text)
    beta <- seq_len(ncol(X)) / ncol(X)
    expect_equal(.msr_response_design_eta(plan, K, beta), drop(X %*% beta),
      tolerance = 1e-12, ignore_attr = TRUE, info = text)
    numerical <- vapply(seq_along(covs), function(j) {
      a <- b <- K; a[, j] <- a[, j] + 1e-6; b[, j] <- b[, j] - 1e-6
      drop((.msr_response_design_matrix(plan, a) - .msr_response_design_matrix(plan, b)) %*% beta) / 2e-6
    }, numeric(n))
    expect_equal(.msr_response_kernel_sensitivity(plan, K, beta), numerical,
      tolerance = 1e-7, ignore_attr = TRUE, info = text)
  }
  expect_no_error(.msr_validate_kernel_expressions(terms(y ~ I(h^2)), "h"))
  expect_error(.msr_validate_kernel_expressions(terms(y ~ I(y^2)), "y"), "response expression")
})

test_that("grouped response products retain complete AD gradients", {
  skip_if_not_installed("RTMB")
  dat <- data.frame(y = seq_len(12), h = sin(1:12), k = cos(1:12),
    f = factor(rep(letters[1:3], 4)))
  mf <- model.frame(y ~ h*k*f, dat)
  plan <- .msr_compile_response_design(terms(mf), mf, c("h", "k"))
  K <- as.matrix(dat[c("h", "k")]); K[1, ] <- 0
  beta <- seq_len(plan$p) / plan$p
  X <- .msr_response_design_matrix(plan, K); eta <- drop(X %*% beta)
  obj <- RTMB::MakeADFun(function(par) {
    v <- .msr_response_design_eta(plan, par$kernel, par$beta)
    sum(v * v) / 2
  }, list(kernel = K, beta = beta), silent = TRUE)
  expect_identical(unname(obj$par), c(as.vector(K), beta))
  expect_equal(obj$fn(obj$par), sum(eta * eta) / 2, tolerance = 1e-12)
  expected <- c(as.vector(.msr_response_kernel_sensitivity(plan, K, beta) * eta),
    drop(crossprod(X, eta)))
  expect_equal(as.numeric(obj$gr(obj$par)), as.numeric(expected), tolerance = 1e-10)
})
