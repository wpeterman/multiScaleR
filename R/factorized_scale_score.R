# Independent cold-audit derivatives from conditional score moments.
.msr_factorized_score_adfun <- function(spec, parameters=spec$parameters,
                                        rel.tol=1e-9, step=1) {
  .msr_factorized_validate(spec)
  parameters$scale_z <- parameters$response_z <- numeric()
  if(!length(spec$scales[[1L]]$random))
    return(.msr_factorized_adfun(spec,parameters,rel.tol,step))
  sizes<-lengths(parameters);blocks<-rep(names(parameters),sizes)
  theta<-unlist(parameters,use.names=FALSE);names(theta)<-blocks
  unpack<-function(x) {
    p<-parameters;at<-0L
    for(i in seq_along(p)) {
      n<-sizes[i];p[[i]]<-if(n)x[at+seq_len(n)] else numeric();at<-at+n
    }
    p
  }
  groups<-lapply(seq_along(spec$scales[[1L]]$random[[1L]]$levels),function(i) {
    g<-.msr_factorized_group(spec,i);p<-parameters;p$scale_z<-0
    obj<-RTMB::MakeADFun(function(q)
      sum(.msr_factorized_response_nll(g,q,TRUE))-RTMB::dnorm(q$scale_z,log=TRUE),p,silent=TRUE)
    list(spec=g,obj=obj,outer=which(names(obj$par)!="scale_z"),z=which(names(obj$par)=="scale_z"))
  })
  state<-new.env(parent=emptyenv());state$theta<-NULL;state$order<--1L
  state$last.par<-theta;state$parList<-function(x,par=NULL)unpack(x)
  evaluate<-function(x,order) {
    if(!identical(as.numeric(x),state$theta)) {
      state$theta<-as.numeric(x);state$order<--1L;state$records<-NULL
      state$last.par<-x
    }
    if(state$order>=order)return(invisible(NULL))
    p<-unpack(x);d<-length(x)
    if(is.null(state$records))state$records<-lapply(groups,function(g) {
      cache<-new.env(hash=TRUE,parent=emptyenv())
      node<-function(z,order=0L) {
        key<-sprintf("%a",z)
        a<-if(exists(key,cache,inherits=FALSE))get(key,cache,inherits=FALSE) else list()
        v<-g$obj$par;v[g$outer]<-x;v[g$z]<-z
        if(is.null(a$value))a$value<-as.numeric(g$obj$fn(v))
        if(order>=1L&&is.null(a$gradient))a$gradient<-as.numeric(g$obj$gr(v))[g$outer]
        if(order>=2L&&is.null(a$hessian))a$hessian<-g$obj$he(v)[g$outer,g$outer,drop=FALSE]
        assign(key,a,cache);a
      }
      bound<-if(spec$family=="gaussian")
        -length(g$spec$y)*(p$log_dispersion+.5*log(2*pi)) else 0
      search<-.msr_factorized_modes(function(z)-node(z)$value,bound)
      shift<--search$best$objective
      cuts<-sort(unique(c(-Inf,seq(-8,8,by=step),
        vapply(search$candidates,`[[`,0,"minimum"),Inf)))
      integrate<-function(fun,absolute)sum(vapply(seq_len(length(cuts)-1L),function(j)
        stats::integrate(function(z)vapply(z,fun,0),cuts[j],cuts[j+1L],rel.tol=rel.tol,
          abs.tol=absolute,subdivisions=500L,stop.on.error=TRUE)$value,0))
      integral<-integrate(function(z)exp(-node(z)$value-shift),0)
      if(!is.finite(integral)||integral<=0)stop("Invalid conditional normalization.",call.=FALSE)
      lognorm<-log(integral)+shift
      weight<-function(z)exp(-node(z)$value-lognorm)
      list(node=node,integrate=integrate,weight=weight,lognorm=lognorm,
        center=node(search$best$minimum,1L)$gradient)
    })
    state$value<--sum(vapply(state$records,`[[`,0,"lognorm"))
    if(order>=1L&&state$order<1L) {
      state$gradient<-numeric(d)
      for(i in seq_along(state$records)) {
        r<-state$records[[i]]
        mean_difference<-vapply(seq_len(d),function(j)r$integrate(function(z) {
          w<-r$weight(z);if(w==0)return(0)
          w*(r$node(z,1L)$gradient[j]-r$center[j])
        },1e-10),0)
        state$records[[i]]$mean_difference<-mean_difference
        state$gradient<-state$gradient+r$center+mean_difference
      }
    }
    if(order>=2L&&state$order<2L) {
      state$hessian<-matrix(0,d,d)
      for(r in state$records) {
        h<-matrix(0,d,d)
        for(j in seq_len(d))for(k in seq_len(j)) {
          h[j,k]<-r$integrate(function(z) {
            w<-r$weight(z);if(w==0)return(0)
            a<-r$node(z,2L);delta<-a$gradient-r$center
            w*(a$hessian[j,k]-delta[j]*delta[k])
          },1e-8)+r$mean_difference[j]*r$mean_difference[k]
          h[k,j]<-h[j,k]
        }
        state$hessian<-state$hessian+h
      }
    }
    state$order<-max(state$order,order)
    invisible(NULL)
  }
  list(par=theta,env=state,fn=function(x){evaluate(x,0L);state$value},
    gr=function(x){evaluate(x,1L);state$gradient},
    he=function(x){evaluate(x,2L);state$hessian},
    hessian_method="Conditional joint-Hessian expectation minus score covariance")
}

