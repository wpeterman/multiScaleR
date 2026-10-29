test_that("conditional mode searches certify exterior tails and expand beyond twelve SDs", {
  for(target in c(-15,15)) {
    sigma <- .01
    log_density <- function(z)dnorm(target,z,sigma,log=TRUE)+dnorm(z,log=TRUE)
    bound <- -.5*log(2*pi)-log(sigma)
    search <- multiScaleR:::.msr_factorized_modes(log_density,bound)
    expect_gt(search$radius,12)
    expect_equal(search$best$minimum,target/(1+sigma^2),tolerance=1e-6)
    expect_lt(search$exterior_log_bound,-search$best$objective)
  }
  ordinary <- multiScaleR:::.msr_factorized_modes(function(z)dnorm(z,log=TRUE),0)
  expect_equal(ordinary$radius,12)
  expect_equal(ordinary$best$minimum,0,tolerance=1e-7)
  expect_error(multiScaleR:::.msr_factorized_modes(function(z)NA_real_,0),"invalid")
  expect_error(multiScaleR:::.msr_factorized_modes(function(z)0,Inf),"finite")
  # An unresolved exterior remains a numerical failure rather than a clipped mode.
  expect_error(multiScaleR:::.msr_factorized_modes(function(z)-1e6,0),"could not certify")
})
