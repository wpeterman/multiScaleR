# Closures stored in returned objects carried source references whenever the
# package was loaded with keep.source = TRUE (covr, interactive load_all()).
# Each serialized copy then held a whole source file: 4.6 MB under covr, which
# failed the transport-size tests in the 0.7.48 coverage run on CI.

test_that(".msr_drop_srcref() strips source references and keeps the closure", {
  fn <- eval(parse(text = "function(x) {\n  x + 1\n}", keep.source = TRUE)[[1L]])
  expect_false(is.null(attr(fn, "srcref")))
  env <- new.env()
  environment(fn) <- env
  attr(fn, "tag") <- "kept"

  out <- .msr_drop_srcref(fn)
  expect_null(attr(out, "srcref"))
  expect_null(attr(body(out), "srcref"))
  expect_null(attr(body(out), "wholeSrcref"))
  expect_null(attr(body(out), "srcfile"))
  expect_identical(environment(out), env)
  expect_identical(attr(out, "tag"), "kept")
  expect_equal(out(1), 2)

  # Closures without source references come back untouched.
  plain <- eval(parse(text = "function(x) x", keep.source = FALSE)[[1L]])
  expect_identical(.msr_drop_srcref(plain), plain)
})

test_that("stored spline programs carry no source references", {
  skip_if_not_installed("lpSolve")
  skip_if_not_installed("gmp")
  programs <- raw_binomial_fixture()$context$design$expression_programs
  closures <- unlist(lapply(programs, function(p) c(list(p$value), p$derivative)))
  expect_true(length(closures) > 0L)
  expect_true(all(vapply(closures, function(g) is.null(attr(g, "srcref")), logical(1))))
})
