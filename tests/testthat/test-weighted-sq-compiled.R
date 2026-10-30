test_that('compiled weighted SQ preserves the centered reference calculation', {
  values <- cbind(a=c(1,2,4,NA,Inf,7), b=1e10+c(1,2,4,NA,6,7), constant=rep(5,6))
  d <- c(0,.1,.5,1,2,3)
  for (kernel in c('gaussian','exp','fixed','expow')) for (s in c(.01,.5,100)) {
    w <- .msr_kernel_shape_weight(d,s,1.3,kernel)
    ref <- vapply(seq_len(ncol(values)),function(j) .surface_weighted_metric(values[,j],w,'sq'),0)
    got <- .surface_weighted_by_buffer(d,values,kernel,s,1.3,'sq',validate=FALSE)
    expect_equal(unname(got),ref,tolerance=1e-12)
    expect_identical(names(got),colnames(values))
  }
  expect_true(is.na(.surface_weighted_by_buffer(rep(100,6),matrix(1:6), 'gaussian',1e-5,NULL,'sq',FALSE)))
  expect_equal(.surface_weighted_by_buffer(0,matrix(7),'gaussian',1,NULL,'sq',FALSE),7*0)
  expect_error(.surface_weighted_by_buffer(0,matrix(7:8),'gaussian',1,NULL,'sq'), 'one row')
  expect_error(.surface_weighted_by_buffer(d,values,'gaussian',1,NULL,'sq'), 'finite|non-finite')
})

test_that("compiled surface values and derivatives are invariant to additive origins", {
  buffer <- list(list(d = c(0, .1, .5, 1), v = c(1, 2, 4, 7)))
  shifted <- buffer
  shifted[[1L]]$v <- shifted[[1L]]$v + 1e10
  for (kernel in c("gaussian", "exp")) for (sigma in c(.01, .5, 100)) {
    expected <- msr_surface_sq_derivative_cpp(buffer, sigma, kernel)
    actual <- msr_surface_sq_derivative_cpp(shifted, sigma, kernel)
    expect_equal(actual, expected, tolerance = 1e-12)
  }
})
