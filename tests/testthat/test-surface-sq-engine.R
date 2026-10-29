test_that("weighted SQ scale derivatives match scalar metrics and finite differences", {
  inputs <- list(scale_vars = msr_vars(a = surface_var("z", "sq", weighted = TRUE),
                                      b = surface_var("other", "sq", weighted = TRUE)),
    d_list = rep(list(c(0, .1, .2, .4, .8, 1.2)), 6),
    raw_cov = lapply(1:6, function(i) data.frame(z = c(1, i, 3, i^2, 4, 5),
                                               other = c(i, 3, NA, 5, i^2, 2))))
  for (kernel in c("gaussian", "exp")) for (standardize in c(FALSE, TRUE)) {
    inputs$kernel <- kernel
    e <- .msr_prepare_surface_engine(inputs, c("b", "a"), c(5, 2, 4, 1), standardize)
    p <- log(c(.4, .7))
    got <- .msr_eval_kernel_engine(e, p)
    streamed <- e; streamed$implementation <- "stream"
    expect_equal(.msr_eval_kernel_engine(streamed, p), got, tolerance = 1e-12)
    expected <- vapply(c(other = .4, z = .7), function(sigma) numeric(4), numeric(4))
    for (j in 1:2) for (i in 1:4) {
      row <- e$rows[i]; source <- c("other", "z")[j]
      expected[i, j] <- .surface_weighted_by_buffer(inputs$d_list[[row]],
        inputs$raw_cov[[row]][, source, drop = FALSE], kernel, exp(p[j]), NULL, "sq", FALSE)[[1]]
    }
    expect_equal(unname(got$raw), unname(expected), tolerance = 1e-12)
    for (j in 1:2) {
      plus <- minus <- p; plus[j] <- plus[j] + 1e-5; minus[j] <- minus[j] - 1e-5
      numerical <- (.msr_eval_kernel_engine(e, plus)$value -
                    .msr_eval_kernel_engine(e, minus)$value) / 2e-5
      expect_equal(got$d_scale[, j], numerical[, j], tolerance = 1e-6)
      expect_equal(numerical[, 3-j], rep(0, e$n))
    }
  }
  inputs$kernel <- "gaussian"
  inputs$scale_vars$weighted[1] <- FALSE
  expect_error(.msr_prepare_surface_engine(inputs, c("a", "b"), 1:6), "weighted SQ")
})

test_that("compiled surface derivative handles constant and tiny support", {
  buffers <- list(list(d = c(0, 1, 1e200), v = c(3, 3, 3)))
  for (kernel in c("gaussian", "exp")) {
    expect_equal(unname(msr_surface_sq_derivative_cpp(buffers, 1e-100, kernel)),
                 matrix(c(0, 0), 1, 2))
    expect_error(msr_surface_sq_derivative_cpp(
      list(list(d = c(1, 2), v = c(1, 2))), 1e-100, kernel), "effective surface support")
  }
  expect_error(msr_surface_sq_derivative_cpp(list(list(d = 1, v = NA_real_)),
                                           1, "gaussian"), "finite cells")
})

test_that("compiled SQ preserves centered variance with large offsets and unequal weights", {
  inputs <- list(kernel = "gaussian",
    scale_vars = msr_vars(a = surface_var("z", "sq", weighted = TRUE),
                         b = surface_var("other", "sq", weighted = TRUE)),
    d_list = rep(list(c(0, 1e-6, .1, 1, 100, Inf)), 4),
    raw_cov = lapply(1:4, function(i) data.frame(
      z = 1e8 + c(0, .01, .03 * i, .04, .05, .06),
      other = c(NA, i, 4, 8, 16, 32))))
  for (kernel in c("gaussian", "exp")) {
    inputs$kernel <- kernel
    for (standardize in c(FALSE, TRUE)) {
      e <- .msr_prepare_surface_engine(inputs, c("b", "a"), c(4, 1, 3, 2),
                                       standardize, "reference")
      stream <- e; stream$implementation <- "stream"
      p <- log(c(.2, .5))
      reference <- .msr_eval_surface_engine(e, p)
      compiled <- .msr_eval_surface_engine(stream, p)
      expect_equal(compiled, reference, tolerance = 1e-12)
      expect_true(all(is.finite(compiled$d_scale)))
      expect_true(all(compiled$raw >= 0))
    }
  }
})

test_that("surface buffer reuse preserves masks and copy-on-write isolation", {
  inputs <- list(kernel = "gaussian",
    scale_vars = msr_vars(a = surface_var("z", "sq", weighted = TRUE)),
    d_list = list(c(0, 1, 2), c(0, 1, 2)),
    raw_cov = list(data.frame(z = c(1, 2, 4)), data.frame(z = c(2, NA, 6))))
  e <- .msr_prepare_surface_engine(inputs, "a", 1:2, FALSE)
  before <- .msr_eval_surface_engine(e, 0)
  expect_identical(e$buffers$a[[1]]$d, inputs$d_list[[1]])
  expect_identical(e$buffers$a[[2]]$v, c(2, 6))
  inputs$d_list[[1]][1] <- 20
  inputs$raw_cov[[1]]$z[1] <- 100
  expect_identical(.msr_eval_surface_engine(e, 0), before)
})
