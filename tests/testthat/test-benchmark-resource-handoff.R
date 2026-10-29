test_that("benchmark resource reads retry without relaxing dispatch gates", {
  root<-normalizePath(file.path("..",".."),winslash="/",mustWork=FALSE)
  helper<-file.path(root,"tools/benchmarks/joint_scale_assessment/phase12/resource.R")
  skip_if_not(file.exists(helper),"Development benchmark helper is not packaged")
  e<-new.env(parent=baseenv())
  sys.source(file.path(root,"tools/benchmarks/joint_scale_assessment/phase11/dispatch.R"),e)
  sys.source(helper,e)
  now<-Sys.time();deadline<-now+60;calls<-0L;sleeps<-0L
  reader<-function(path) {calls<<-calls+1L;if(calls<3L)stop("locked");c(as.numeric(now),12)}
  pause<-function(x)sleeps<<-sleeps+1L
  z<-e$resource_state(tempdir(),read_record=reader,sleep=pause)
  expect_identical(calls,3L);expect_identical(sleeps,2L)
  expect_identical(e$dispatch_gate(deadline,z),"ok")
  for(record in list(numeric(),1,c(NA,12),c(Inf,12),c(1,2,3))) {
    calls<-0L
    z<-e$resource_state(tempdir(),read_record=function(path){calls<<-calls+1L;record},sleep=function(x)NULL)
    expect_identical(calls,5L)
    expect_identical(e$dispatch_gate(deadline,z),"stale_resource_monitor")
  }
  record<-function(age,ram)c(as.numeric(Sys.time())-age,ram)
  read<-function(x)e$resource_state(tempdir(),read_record=function(path)x,sleep=function(x)NULL)
  expect_identical(e$dispatch_gate(deadline,read(record(21,12))),"stale_resource_monitor")
  expect_identical(e$dispatch_gate(deadline,read(record(-10,12))),"stale_resource_monitor")
  expect_identical(e$dispatch_gate(deadline,read(record(0,7.9))),"low_ram")
  expect_identical(e$dispatch_gate(deadline,read(record(0,12)),TRUE),"stop_requested")
  expect_identical(e$dispatch_gate(now-1,read(record(0,12))),"deadline")
  # A valid low-RAM or stale record is returned immediately, never retried away.
  calls<-0L
  z<-e$resource_state(tempdir(),read_record=function(path){calls<<-calls+1L;record(0,7)},sleep=pause)
  expect_identical(calls,1L);expect_identical(e$dispatch_gate(deadline,z),"low_ram")
  expect_identical(e$dispatch_gate(deadline,e$resource_state(tempfile(),delay=0)),"stale_resource_monitor")
})
