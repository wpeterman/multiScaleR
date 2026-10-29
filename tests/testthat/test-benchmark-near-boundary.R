boundary_helpers <- function() {
  root <- file.path('..','..','tools/benchmarks/joint_scale_assessment/phase18')
  skip_if_not(file.exists(file.path(root,'boundary-policy/policy.R')))
  e <- new.env(parent=globalenv())
  sys.source(file.path(root,'policy.R'),e)
  sys.source(file.path(root,'boundary-policy/policy.R'),e); e
}
boundary_example <- function() {
  info <- list(available=FALSE,boundary=TRUE,weak=TRUE,multiple_modes=TRUE,
    poisson_limit=FALSE,poisson_compatible=FALSE,reasons=c('parameter at an optimization bound',
      'log-kernel standard error exceeds one'))
  other <- info; other$boundary <- FALSE; other$reasons <- other$reasons[-1]
  fit <- list(best=list(value=0,par=0),records=data.frame(eligible=c(TRUE,TRUE),value=c(0,0)),
    fits=list(list(par=0),list(par=.2)))
  warm <- fit; warm$best$par <- .2
  list(reference=list(full=fit),warm=list(full=warm),audit=list(reference=list(information=info),warm=list(information=other)))
}
test_that('near-optimal boundary comparison preserves raw selected flags', {
  e <- boundary_helpers(); z <- boundary_example()
  x <- e$compare_boundary_policy(z,0,1)
  expect_true(x$passed); expect_true(x$explained); expect_false(x$original_match)
  expect_identical(x$raw_information$reference,z$audit$reference$information)
  expect_identical(x$raw_information$warm,z$audit$warm$information)
  expect_true(x$raw_information$reference$boundary)
  expect_false(x$raw_information$warm$boundary)
  expect_identical(x$effective_available,list(reference=FALSE,warm=FALSE))
})
test_that('availability, Poisson and distinct boundary evidence cannot be excused', {
  e <- boundary_helpers()
  z <- boundary_example(); z$audit$warm$information$available <- TRUE
  expect_false(e$compare_boundary_policy(z,0,1)$passed)
  z <- boundary_example(); z$audit$warm$information$poisson_limit <- TRUE
  expect_false(e$compare_boundary_policy(z,0,1)$passed)
  z <- boundary_example(); z$warm$full$fits[[1]]$par <- 1
  expect_false(e$compare_boundary_policy(z,0,1)$passed)
  z <- boundary_example(); z$warm$full$records$eligible[1] <- FALSE
  expect_false(e$compare_boundary_policy(z,0,1)$passed)
  z <- boundary_example(); z$warm$full$records$value[1] <- .1
  expect_false(e$compare_boundary_policy(z,0,1)$passed)
  z <- boundary_example(); z$audit$warm$information$reasons <- c(z$audit$warm$information$reasons,'unknown failure')
  expect_false(e$compare_boundary_policy(z,0,1)$passed)
  z <- boundary_example(); z$warm$full$records$eligible[] <- FALSE
  expect_error(e$compare_boundary_policy(z,0,1),'No eligible')
})
