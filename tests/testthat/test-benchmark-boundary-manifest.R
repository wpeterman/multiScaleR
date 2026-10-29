test_that('all-empty archive columns stay character for new reference fitting', {
  path <- file.path('..','..','tools/benchmarks/joint_scale_assessment/phase18/boundary-policy/revision2/manifest.R')
  skip_if_not(file.exists(path)); e <- new.env(); sys.source(path,e)
  f <- tempfile(fileext='.csv'); on.exit(unlink(f),add=TRUE)
  x <- data.frame(archive=c('',''),old_path=c(NA_character_,NA_character_),
    path=c('a','b'),fixture=c('x','y'),parent_path=c('p','q'),index=1:2)
  write.csv(x,f,row.names=FALSE); j <- e$read_boundary_jobs(f)
  expect_identical(j$archive,c('','')); expect_false(any(nzchar(j$archive)))
  expect_true(is.numeric(j$index))
  x$archive[2] <- 'saved.rds'; write.csv(x,f,row.names=FALSE)
  expect_identical(e$read_boundary_jobs(f)$archive,c('','saved.rds'))
  x$archive[1] <- NA_character_; write.csv(x,f,row.names=FALSE)
  expect_error(e$read_boundary_jobs(f),'Archive paths must be character')
})
