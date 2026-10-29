test_that('weighted reuse requires exact source, scale and shape agreement', {
 s<-msr_vars(a=surface_var('z','sa',weighted=TRUE),b=surface_var('z','sq',weighted=TRUE),c=surface_var('z','ssk',weighted=TRUE),e=surface_var('z','sku',weighted=TRUE))
 d<-c(0,.1,.5,1,2,3);v<-data.frame(z=c(1,2,4,NA,6,9))
 run<-function(sig,shape=NULL,kernel='gaussian',spec=s,values=v) .msr_shared_weighted_metrics(d,values,spec,spec$covariate,sig,shape,kernel,spec$covariate,FALSE)
 for(k in c('gaussian','exp','fixed','expow')) {
  got<-run(rep(.7,4),rep(1.3,4),k)
  ref<-vapply(s$metric,function(m).surface_weighted_by_buffer(d,v,k,.7,1.3,m,FALSE)[[1]],0)
  expect_equal(unname(got),unname(ref),tolerance=1e-12)
 }
 expect_length(run(c(.7,.7+1e-14,.8,.9)),0)
 expect_length(run(rep(.7,4),1:4,'expow'),0)
 different<-s;different$source<-c('z','z2','z3','z4')
 v$z2<-v$z;v$z3<-replace(v$z,1,NA);v$z4<-replace(v$z,2,NA)
 expect_length(run(rep(.7,4),spec=different,values=v),0)
 # Same source name on a later buffer must recompute, including a changed mask.
 changed<-v['z'];changed$z[1]<-NA
 got<-run(rep(.7,4),values=changed)
 ref<-vapply(s$metric,function(m).surface_weighted_by_buffer(d,changed,'gaussian',.7,NULL,m,FALSE)[[1]],0)
 expect_equal(unname(got),unname(ref),tolerance=1e-12)
 for(z in list(rep(NA_real_,6),rep(4,6),1e10+c(1,2,3,4,5,6))) {
  got<-run(rep(.7,4),values=data.frame(z=z))
  ref<-vapply(s$metric,function(m).surface_weighted_by_buffer(d,matrix(z),'gaussian',.7,NULL,m,FALSE)[[1]],0)
  expect_equal(unname(got),unname(ref),tolerance=1e-12)
 }
 duplicate<-s[1:2,];duplicate$metric<-'sq'
 expect_equal(unname(run(rep(.7,2),spec=duplicate)),rep(.surface_weighted_by_buffer(d,v['z'],'gaussian',.7,NULL,'sq',FALSE)[[1]],2))
 expect_true(all(is.na(run(rep(1e-10,4),values=data.frame(z=c(NA,1:5))))))
 expect_error(.msr_shared_weighted_metrics(d,data.frame(z=c(Inf,1:5)),s,s$covariate,rep(.7,4),NULL,'gaussian',s$covariate,TRUE),'finite')
})

