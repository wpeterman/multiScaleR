curve_numeric_reference <- function(a, b) is.numeric(a) && is.numeric(b) &&
  identical(dim(a), dim(b)) && all(is.finite(a)) && all(is.finite(b)) &&
  isTRUE(all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE))

test_that("reference numeric equality retains complete tolerance decisions", {
  compare <- multiScaleR:::.msr_curve_numeric_equivalent
  for (n in c(0L, 1L, 4095L, 4096L, 5000L)) {
    a <- seq_len(n) / 100
    for (x in list(a, as.integer(a), matrix(a, nrow = n, ncol = 1L))) {
      expect_identical(compare(x, x), curve_numeric_reference(x, x))
      if (length(x)) for (delta in c(1e-12, 1e-10, 1e-6)) {
        b <- x; b[length(b)] <- b[length(b)] + delta
        expect_identical(compare(x, b), curve_numeric_reference(x, b))
      }
    }
  }
  a <- rep(1, 4096)
  for (position in c(1L, 4095L, 4096L)) for (ratio in c(1 - 1e-6, 1 + 1e-6)) {
    b <- a; b[position] <- b[position] + 1e-10 * ratio
    expect_identical(compare(a, b), curve_numeric_reference(a, b))
  }
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    b <- a; b[4096L] <- bad
    expect_false(compare(a, b))
    expect_false(compare(b, b))
  }
  expect_true(compare(rep(0, 4096), rep(-0, 4096)))
  expect_identical(compare(1:5, 1:6), curve_numeric_reference(1:5, 1:6))
  expect_false(compare(a, as.character(a)))
  a <- matrix(a, 1024L, 4L)
  expect_false(compare(a, t(a)))
  b <- a; attr(b, "extra") <- "not part of numeric equality"
  expect_identical(compare(a, b), curve_numeric_reference(a, b))
  b <- a; colnames(b) <- letters[1:4]
  expect_identical(compare(a, b), curve_numeric_reference(a, b))
  b <- a + 0i
  expect_identical(compare(b, b), curve_numeric_reference(b, b))
})

# Restore exactly the pre-test method table, including the absence of a method.
with_curve_comparison_method <- function(kind, method, expr) {
  table <- get(".__S3MethodsTable__.", envir = asNamespace("base"))
  key <- paste0("all.equal.", kind)
  existed <- exists(key, envir = table, inherits = FALSE)
  previous <- if (existed) get(key, envir = table, inherits = FALSE) else NULL
  on.exit(if (existed) assign(key, previous, envir = table) else
    rm(list = key, envir = table), add = TRUE)
  registerS3method("all.equal", kind, method, envir = asNamespace("base"))
  force(expr)
}

test_that("reference shortcuts preserve custom numeric comparison dispatch", {
  compare <- multiScaleR:::.msr_curve_numeric_equivalent
  a <- matrix(seq_len(8192), 2048L, 4L)
  with_curve_comparison_method("matrix", function(...) "custom mismatch", {
    expect_false(curve_numeric_reference(a, a))
    expect_false(compare(a, a))
  })
  with_curve_comparison_method("matrix", function(...) stop("custom comparison error"), {
    expect_error(curve_numeric_reference(a, a), "custom comparison error")
    expect_error(compare(a, a), "custom comparison error")
  })
  a <- structure(a, class = "msr_numeric_probe")
  with_curve_comparison_method("msr_numeric_probe", function(...) "custom class mismatch", {
    expect_false(curve_numeric_reference(a, a))
    expect_false(compare(a, a))
  })
})

test_that("large raw-reference designs retain native rows and complete verification", {
  n <- 2048L
  dat <- data.frame(z = seq(-2, 2, length.out = n),
    group = factor(rep(c("a", "b"), length.out = n)),
    exposure = 1 + seq_len(n) / n)
  dat$y <- .2 * dat$z + .3 * (dat$group == "b") + .1 * dat$z^2 +
    log(dat$exposure) + sin(seq_len(n) * sqrt(2))
  model <- lm(y ~ z * group + I(z^2) + offset(log(exposure)), data = dat)
  object <- structure(list(opt_mod = model,
    scl_params = list(mean = numeric(), sd = numeric())), class = "multiScaleR")
  reference <- dat; row.names(reference) <- row.names(model$model)
  original <- multiScaleR:::.msr_curve_numeric_equivalent
  large_exact <- 0L
  testthat::local_mocked_bindings(.msr_curve_numeric_equivalent = function(a, b) {
    if (is.matrix(a) && length(a) >= 4096L && identical(a, b))
      large_exact <<- large_exact + 1L
    original(a, b)
  }, .package = "multiScaleR")
  # A formula environment points back to this evolving test scope. Preserve
  # its identity separately and compare every other fitted-model byte.
  environment_identity <- function(x) list(model = attr(x$terms, ".Environment"),
    frame = attr(attr(x$model, "terms"), ".Environment"))
  stable_model <- function(x) {
    attr(x$terms, ".Environment") <- baseenv()
    frame_terms <- attr(x$model, "terms")
    attr(frame_terms, ".Environment") <- baseenv()
    attr(x$model, "terms") <- frame_terms
    serialize(x, NULL)
  }
  before <- stable_model(model)
  environments_before <- environment_identity(model)
  value <- multiScaleR:::.msr_curve_reference_context(object, model, reference)
  expect_gt(large_exact, 0L)
  expect_true(value$contract$native_design_and_offsets_verified)
  expect_identical(value$contract$equivalence_tolerance, 1e-10)
  expect_identical(row.names(value$data), row.names(model$model))
  expect_identical(stable_model(model), before)
  expect_identical(environment_identity(model), environments_before)
  reference$z[n] <- reference$z[n] + .01
  expect_error(multiScaleR:::.msr_curve_reference_context(object, model, reference),
    "fitted feature")
})
