test_that("prepared SQ preserves radius boundaries and finite-cell masks", {
  specs <- msr_vars(a = surface_var("z", "sq"),
                    b = surface_var("z", "sq"))
  d <- c(2, 0, 1, 1, 3, NA, Inf, 4)
  for (z in list(c(2, 1, NA, 4, 8, 3, 4, 5), rep(4, 8),
                 1e10 + c(1, 2, 3, 4, 5, 6, 7, 8))) {
    values <- data.frame(z = z)
    cache <- .msr_prepare_surface_sq(list(d), list(values), specs)[[1]]
    expect_identical(cache$a, cache$b)
    for (radius in c(.01, .99, 1, 1.01, 2, 3, 4, 20)) {
      reference <- .surface_metric_by_buffer(d, values, radius, "sq", FALSE)[[1]]
      expect_equal(.msr_eval_surface_neighbors(cache$a, radius), reference,
                   tolerance = 1e-10)
    }
  }
  expect_null(.msr_prepare_surface_sq(list(d), list(data.frame(z = 1:8)),
                                     specs, budget = 1)[[1]])
  fixed <- msr_vars(a = surface_var("z", "sq", radius = 100))
  expect_null(.msr_prepare_surface_sq(list(d), list(data.frame(z = 1:8)), fixed))
  weighted <- msr_vars(a = surface_var("z", "sq", weighted = TRUE))
  expect_null(.msr_prepare_surface_sq(list(d), list(data.frame(z = 1:8)), weighted))
})

test_that("SQ context maps independent radii without changing scalar results", {
  old <- options(multiScaleR.experimental_sq_cache = TRUE)
  on.exit(options(old), add = TRUE)
  specs <- msr_vars(a = surface_var("z", "sq"), b = surface_var("other", "sq"))
  d <- c(0, .2, .4, .4, .8, 1)
  v <- data.frame(z = c(2, 3, NA, 5, 8, 9), other = c(NA, 2, 3, 4, 5, 6))
  cache <- .msr_prepare_surface_neighbors(list(d), list(v), specs, 1, 3)[[1]]
  run <- function(prepared) .msr_eval_scale_vars(d, v, specs, c(.4, .8), NULL,
    "gaussian", 1, 1, 3, validate = FALSE, prepared_neighbors = prepared)
  expect_equal(run(cache), run(NULL), tolerance = 1e-12)
})

test_that("experimental SQ cache leaves ordinary context preparation unchanged", {
  old <- options(multiScaleR.experimental_sq_cache = FALSE)
  on.exit(options(old), add = TRUE)
  specs <- msr_vars(a = surface_var("z", "sq"))
  d <- list(c(0, 1, 2)); v <- list(data.frame(z = c(1, 3, 2)))
  expect_identical(.msr_prepare_surface_neighbors(d, v, specs, 1, 3),
                   .msr_prepare_surface_neighbors_base(d, v, specs, 1, 3))
})
