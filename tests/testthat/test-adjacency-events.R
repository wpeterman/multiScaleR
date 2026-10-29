test_that("entry-distance adjacency preserves exact cutoffs and class order", {
  # Unsorted distances change which class appears first as the radius grows.
  d <- c(.8,.2,.4,.2,0,.9,.4,.6,.7,.3,.5,.1)
  cells <- c(1:8,10,11,12,16)
  values <- c(3,2,1,3,1,2,3,1,NA,2,1,3)
  p <- msr_adjacency_prepare_cpp(d,values,cells,4L)
  events <- msr_adjacency_events_prepare_cpp(p)
  stored <- serialize(events,NULL)
  metrics <- c(.msr_adjacency_metrics(),"clumpy")
  for (radius in c(.01,.1,.2,.200000001,.399999999,.4,.8,1)) for (base in c(2,exp(1))) {
    args <- list(radii=rep(radius,length(metrics)),metrics=metrics,
      bases=rep(base,length(metrics)),focals=rep(2,length(metrics)))
    expect_identical(do.call(msr_adjacency_events_eval_cpp,c(list(events),args)),
                     do.call(msr_adjacency_eval_cpp,c(list(p),args)))
  }
  expect_identical(serialize(events,NULL),stored)
  expect_identical(.msr_adjacency_eval(unserialize(stored),c(.2,.8),c("ai","ent"),c(2,2),c(NA,NA)),
                   .msr_adjacency_eval(p,c(.2,.8),c("ai","ent"),c(2,2),c(NA,NA)))
  empty <- msr_adjacency_events_prepare_cpp(msr_adjacency_prepare_cpp(NA_real_,NA_real_,1,4L))
  expect_true(is.na(.msr_adjacency_eval(empty,1,"ent",2,NA_real_)))
  bad <- p;bad$right[1] <- length(p$d)
  expect_error(msr_adjacency_events_prepare_cpp(bad),"index")
  bad <- events;bad$pair_offset[1] <- -1L
  expect_error(.msr_adjacency_eval(bad,1,"ai",2,NA_real_),"dimensions")
  expect_error(.msr_adjacency_eval(events,0,"ai",2,NA_real_),"radius")
  expect_error(.msr_adjacency_eval(events,1,"ent",1,NA_real_),"base")
  expect_error(.msr_adjacency_eval(events,1,"clumpy",2,NA_real_),"class")
})
