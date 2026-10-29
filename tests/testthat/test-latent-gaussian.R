test_that("proper Gaussian blocks have normalized densities and exact coordinate maps", {
  A <- matrix(c(1, 0, .5, 0, 1, -.3), 3, 2)
  C <- matrix(c(1, .4, .4, 2), 2)
  block <- multiScaleR:::.msr_gaussian_block(A, C)
  z <- c(.3, -.7); log_sd <- log(.6)
  b <- .6 * drop(t(chol(C)) %*% z)
  expect_equal(multiScaleR:::.msr_gaussian_effect(block, z, log_sd), drop(A %*% b))
  dense <- length(b) / 2 * log(2*pi) + .5 * as.numeric(determinant(.6^2*C, logarithm=TRUE)$modulus) +
    .5 * drop(t(b) %*% solve(.6^2*C, b))
  nll <- multiScaleR:::.msr_gaussian_centered_nll(block, b, log_sd)
  expect_equal(nll, dense, tolerance = 1e-12)
  jacobian <- length(z) * log_sd + sum(log(diag(t(chol(C)))))
  expect_equal(nll - jacobian, -sum(dnorm(z, log = TRUE)), tolerance = 1e-12)
  expect_equal(multiScaleR:::.msr_gaussian_effect(block, z, -Inf), rep(0, 3))
  expect_error(multiScaleR:::.msr_gaussian_block(A, matrix(1,2,2)), "positive definite")
  expect_error(multiScaleR:::.msr_gaussian_block(A, matrix(c(1,2,0,1),2)), "symmetric")
  expect_error(multiScaleR:::.msr_gaussian_block(A, diag(c(1,1e-14))), "well conditioned")
  expect_error(multiScaleR:::.msr_gaussian_effect(block, 1, log_sd), "dimension")
  expect_error(multiScaleR:::.msr_gaussian_block(replace(A,1,NA)), "finite")
})

test_that("Gaussian block derivatives match independent dense precision calculations", {
  skip_if_not_installed("RTMB")
  A <- diag(3); C <- matrix(c(2,.2,.1,.2,1,.3,.1,.3,1.5),3)
  block <- multiScaleR:::.msr_gaussian_block(A,C)
  Q <- solve(C); logdet <- .5 * as.numeric(determinant(C, logarithm=TRUE)$modulus)
  independent <- function(p) {
    b <- p$b; ell <- p$ell
    .5 * exp(-2*ell) * sum(b * drop(Q %*% b)) + 3*ell + 1.5*log(2*pi) + logdet
  }
  p <- list(b = c(.2,-.4,.3), ell = log(.7))
  first <- RTMB::MakeADFun(function(p) multiScaleR:::.msr_gaussian_centered_nll(block,p$b,p$ell),p,silent=TRUE)
  second <- RTMB::MakeADFun(independent,p,silent=TRUE)
  for (quantity in c("fn","gr","he"))
    expect_equal(first[[quantity]](first$par),second[[quantity]](second$par),tolerance=1e-10)
  expect_equal(as.numeric(first$gr(first$par))[1:3], exp(-2*p$ell)*drop(Q%*%p$b),tolerance=1e-10)
  expect_equal(first$he(first$par)[1:3,1:3], exp(-2*p$ell)*Q,tolerance=1e-10)
})

test_that("IID projection adapters reproduce scale and response grouping effects", {
  f <- bayesian_scale_fixture(n_site=4,n_year=3)
  original <- prepare_bayesian_fixture(f)$specification
  p <- original$parameters; p$scale_z[] <- c(-.2,.1,.3)
  p$response_z[] <- seq_along(p$response_z)/10
  expected <- multiScaleR:::.msr_bayes_components(original,p)
  projected <- original
  for (component in c("scale","response")) {
    terms <- if(component=="scale") projected$scales$habitat$random else projected$response$random
    for(k in seq_along(terms)) {
      term <- terms[[k]]
      A <- diag(length(term$levels))[term$index,,drop=FALSE]
      term$gaussian_block <- multiScaleR:::.msr_gaussian_block(A)
      terms[[k]] <- term
    }
    if(component=="scale") projected$scales$habitat$random <- terms else projected$response$random <- terms
  }
  expect_true(multiScaleR:::.msr_has_gaussian_blocks(projected))
  expect_false(multiScaleR:::.msr_has_gaussian_blocks(original))
  expect_equal(multiScaleR:::.msr_bayes_components(projected,p),expected,tolerance=1e-12)
  for(ad in c(FALSE,TRUE)) {
    if(ad) {
      skip_if_not_installed("RTMB")
      a <- RTMB::MakeADFun(multiScaleR:::.msr_bayes_objective(original,TRUE),p,silent=TRUE)
      b <- RTMB::MakeADFun(multiScaleR:::.msr_bayes_objective(projected,TRUE),p,silent=TRUE)
      for(what in c("fn","gr","he")) expect_equal(a[[what]](a$par),b[[what]](b$par),tolerance=1e-10)
    } else expect_equal(multiScaleR:::.msr_bayes_objective(original)(p),multiScaleR:::.msr_bayes_objective(projected)(p))
  }
})
