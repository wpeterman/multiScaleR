# Compiler-free joint likelihood. Uses the optional RTMB dependency;
# no changes to the Bayesian objective, priors, or parameterization.
.msr_mixed_rtmb <- function(engine, X, Xzi, y, offset, group, parameters,
                           family = c('nbinom2','poisson'), fixed_scale=FALSE,
                           response_design=NULL) {
 family<-match.arg(family)
 stopifnot(engine$kernel %in% c('gaussian','exp'),isTRUE(engine$standardize),
  nrow(X)==length(y),nrow(Xzi)==length(y),length(offset)==length(y),
  engine$n==length(y),length(group)==length(y),all(group>=1L),
  all(group<=length(parameters$b)),all(is.finite(y)),all(y>=0),all(y==floor(y)))
 power<-if(engine$kernel=='gaussian')2 else 1
 means<-if(inherits(engine,'msr_surface_engine'))engine$means else engine
 distance<-if(!is.null(means))matrix(rep(means$distance^power,each=engine$n),engine$n) else NULL
 n<-engine$n;p<-engine$p
 if(!is.null(response_design)) {
  stopifnot(identical(dim(response_design$multipliers),c(n,response_design$p)),
    response_design$n==n,response_design$q==p,
    length(parameters$beta)==response_design$p)
  # Additional response-design planning allowance, not a peak-tape bound.
  total<-.msr_total_ram_bytes()
  estimate<-1024^3+4*as.numeric(utils::object.size(response_design))+
    512*as.double(n)*response_design$p
  if(!is.finite(total)||total<=0)stop('Cannot detect RAM for expanded mixed differentiation.',call.=FALSE)
  if(estimate>.75*total)stop('Expanded mixed response exceeds the conservative differentiation RAM budget.',call.=FALSE)
 }
 # Every builder has its own cold inner optimizer and scale tape.
 fn<-function(par) {
  expanded<-!is.null(response_design)
  eta<-if(expanded)offset+par$b[group] else
    drop(X %*% par$beta[seq_len(ncol(X))])+offset+par$b[group]
  kernel<-if(expanded)stats::setNames(vector('list',p),engine$covariates) else NULL
  for(j in seq_len(p)) {
   name<-engine$covariates[j]
   if(inherits(engine,'msr_surface_engine') && name %in% names(engine$buffers)) {
    # Center before squaring; RTMB differentiates both the weighted mean and SQ.
    x<-do.call(c,lapply(engine$buffers[[name]],function(buffer) {
     if(all(buffer$v==buffer$v[1L]))return(0*par$kernel_par[j])
     w<-exp(-(buffer$d^power)*exp(-power*par$kernel_par[j])/if(power==2)2 else 1)
     mu<-sum(w*buffer$v)/sum(w)
     sqrt(sum(w*(buffer$v-mu)^2)/sum(w))
    }))
   } else {
    k<-match(name,means$covariates)
    w<-exp(-distance*exp(-power*par$kernel_par[j])/if(power==2)2 else 1)
    x<-RTMB::rowSums(w*means$V[[k]])/RTMB::rowSums(w*means$C[[k]])
   }
   centered<-x-sum(x)/n
   x<-centered/sqrt(sum(centered^2)/(n-1))
   if(expanded)kernel[[j]]<-x else eta<-eta+x*par$beta[ncol(X)+j]
  }
  if(expanded)eta<-eta+.msr_response_design_eta(response_design,kernel,par$beta)
  ll<-if(family=='nbinom2')RTMB::dnbinom_robust(y,eta,2*eta-par$log_theta,log=TRUE)
      else RTMB::dpois(y,exp(eta),log=TRUE)
  if(ncol(Xzi)>0L) {
   zeta<-drop(Xzi %*% par$gamma)
   lp<-RTMB::dbinom_robust(rep(1,n),rep(1,n),zeta,log=TRUE)
   lq<-RTMB::dbinom_robust(rep(0,n),rep(1,n),zeta,log=TRUE)
   zero<-which(y==0)
   ll<-ll+lq
   ll[zero]<-RTMB::logspace_add(lp[zero],ll[zero])
  }
  -sum(ll)-sum(RTMB::dnorm(par$b,0,exp(par$log_sd),log=TRUE))
 }
 map<-list()
 if(fixed_scale)map$kernel_par<-factor(rep(NA,p))
 if(family=='poisson')map$log_theta<-factor(NA)
 RTMB::MakeADFun(fn,parameters,random='b',map=map,silent=TRUE)
}


# Bridge from a validated native fit to the joint RTMB objective.
# Deliberately narrow: log Poisson/NB2, one conditional random intercept,
# constant dispersion, and native numeric kernel products. Native
# fitted objects remain responsible for final prediction and covariance.
.msr_native_mixed_design <- function(model, engine, kernel_names,
                                        kernel_par, engine_row_ids, check_kernel = TRUE) {
 fail <- function(message) stop(message, call. = FALSE)
 if (!inherits(model, 'glmmTMB')) fail('Expected a glmmTMB fit.')
 if (model$fit$convergence != 0L || !isTRUE(model$sdr$pdHess) ||
     !is.finite(as.numeric(stats::logLik(model)))) fail('Native fit is not eligible.')
 info <- model$modelInfo
 fam <- stats::family(model)
 if (!fam$family %in% c('poisson', 'nbinom2') || fam$link != 'log')
   fail('Only log Poisson and NB2 are supported.')
 if (isTRUE(info$REML) || length(info$map) || length(info$priors))
   fail('REML, parameter maps, and priors are unsupported.')
 re <- info$reTrms$cond
 if (length(re$cnms) != 1L ||
     !identical(unname(re$cnms[[1L]]), '(Intercept)') ||
     length(info$reTrms$zi$cnms)) fail('One conditional random intercept is required.')
 if (length(info$reStruc$condReStruc) != 1L ||
     info$reStruc$condReStruc[[1L]]$blockNumTheta != 1L)
   fail('Unsupported random-effect structure.')
 mf <- stats::model.frame(model)
 if (!identical(as.character(engine_row_ids), rownames(mf)) || engine$n != nrow(mf))
   fail('Kernel rows must exactly match the retained native model rows.')
 weights <- stats::model.weights(mf)
 if (!is.null(weights) && any(weights != 1)) fail('Non-unit weights are unsupported.')
 y <- stats::model.response(mf)
 if (!is.numeric(y) || !is.null(dim(y))) fail('A numeric count response is required.')
 if (length(kernel_names) != engine$p || anyDuplicated(kernel_names))
   fail('Kernel names do not match the engine.')
 # Native fixed terms retain expanded products and quoted predictor names.
 tt <- info$reTrms$cond$terms$fixed
 .msr_validate_kernel_expressions(tt, kernel_names)
 offsets <- attr(tt, 'offset')
 if (length(offsets) && any(kernel_names %in%
     all.vars(as.call(c(quote(list), as.list(attr(tt, 'variables'))[-1L][offsets])))))
   fail('Kernel-dependent offsets are unsupported.')
 if (any(kernel_names %in% all.vars(info$allForm$ziformula)))
   fail('Kernel-dependent zero inflation is unsupported.')
 dt <- info$reTrms$disp$terms$fixed
 if (fam$family == 'nbinom2' &&
     (length(attr(dt, 'term.labels')) || !identical(attr(dt, 'intercept'), 1L)))
   fail('Only constant dispersion is supported.')
 if (length(attr(info$reTrms$zi$terms$fixed, 'offset')))
   fail('Zero-inflation offsets are unsupported.')
 Xfull <- stats::model.matrix(model, component = 'cond')
 Xzi <- stats::model.matrix(model, component = 'zi')
 expanded <- .msr_response_design_required(tt, kernel_names, colnames(Xfull))
 response_design <- if(expanded) {
   frame <- .msr_response_frame_subset(attr(mf,'terms'),mf,tt)
   .msr_compile_response_design(tt,frame,kernel_names,columns=colnames(Xfull),
     contrasts=info$contrasts)
 } else NULL
 kv <- .msr_eval_kernel_engine(engine, unname(kernel_par))$value
 expected <- if(expanded).msr_response_design_matrix(response_design,kv) else kv
 observed <- if(expanded)Xfull else Xfull[,kernel_names,drop=FALSE]
 if (check_kernel && !(identical(dim(observed),dim(expected)) &&
     isTRUE(all.equal(as.numeric(observed),as.numeric(expected),tolerance=1e-10))))
   fail('Kernel values or their row order disagree with the native design.')
 static <- if(expanded)lengths(response_design$dependence)==0L else
   !colnames(Xfull)%in%kernel_names
 X <- Xfull[,static,drop=FALSE]
 fx <- glmmTMB::fixef(model)
 order <- if(expanded)colnames(Xfull) else c(colnames(X), kernel_names)
 if (!setequal(names(fx$cond), order) || !identical(names(fx$zi), colnames(Xzi)))
   fail('Native coefficients do not match the retained design.')
 group_name <- names(re$cnms)[1L]
 random <- glmmTMB::ranef(model)$cond[[group_name]]
 group <- match(as.character(re$flist[[group_name]]), rownames(random))
 if (anyNA(group) || length(group) != length(y)) fail('Invalid random-effect row mapping.')
 offset <- stats::model.offset(mf)
 if (is.null(offset)) offset <- rep(0, length(y))
 parameters <- list(kernel_par = unname(kernel_par), beta = unname(fx$cond[order]),
   gamma = unname(fx$zi), b = unname(random[, 1L]),
   log_sd = log(unname(attr(glmmTMB::VarCorr(model)$cond[[group_name]], 'stddev'))),
   log_theta = if (fam$family == 'nbinom2') log(stats::sigma(model)) else 0)
 out <- list(engine = engine, X = X, Xzi = Xzi, y = y, offset = offset,
      group = group, parameters = parameters, family = fam$family)
 if(!is.null(response_design)) {
   attr(response_design$native_terms,'.Environment')<-baseenv()
   out$response_design<-response_design
 }
 out
}


# Internal multi-start fitter. All scales are in engine log-distance units.
# native_refit must build a fresh glmmTMB fit, including its full Hessian,
# at supplied scales. It must not carry fitted starts between calls.
# Its full covariance is retained for audits, not calibrated scale inference.
.msr_native_joint <- function(engine, kernel_names, engine_row_ids,
                             starts, lower, upper, native_refit) {
 p <- engine$p
 if (!is.matrix(starts) || !is.numeric(starts) || ncol(starts) != p ||
     nrow(starts) < 1L || any(!is.finite(starts)) ||
     length(lower) != p || length(upper) != p ||
     any(!is.finite(c(lower, upper))) || any(lower >= upper) ||
     any(sweep(starts, 2L, lower, '<')) ||
     any(sweep(starts, 2L, upper, '>')) || !is.function(native_refit))
   stop('Supply finite independent log-scale starts, bounds, and a cold native refit function.', call. = FALSE)
 eligible_native <- function(m) inherits(m, 'glmmTMB') &&
   identical(m$fit$convergence, 0L) && isTRUE(m$sdr$pdHess) &&
   is.finite(as.numeric(stats::logLik(m)))
 outcomes <- vector('list', nrow(starts))
 for (i in seq_len(nrow(starts))) {
   warnings <- character()
   outcomes[[i]] <- tryCatch(withCallingHandlers({
     initial <- native_refit(unname(starts[i, ]))
     args <- .msr_native_mixed_design(initial, engine, kernel_names,
                                          starts[i, ], engine_row_ids)
     obj <- do.call(.msr_mixed_rtmb, args)
     ki <- which(names(obj$par) == 'kernel_par')
     if (length(ki) != p) stop('Joint scale parameter mapping failed.')
     lo <- rep(-Inf, length(obj$par)); hi <- rep(Inf, length(obj$par))
     lo[ki] <- lower; hi[ki] <- upper
     opt <- stats::optim(obj$par, obj$fn, obj$gr, method = 'L-BFGS-B',
       lower = lo, upper = hi,
       control = list(maxit = 300, factr = 1e5, pgtol = 1e-7))
     kp <- unname(opt$par[ki])
     value <- obj$fn(opt$par)
     sd <- TMB::sdreport(obj, par.fixed = opt$par)
     terminal <- native_refit(kp)
     # Revalidate the terminal design too; callbacks must preserve contracts.
     terminal_args <- .msr_native_mixed_design(terminal, engine, kernel_names,
                                                   kp, engine_row_ids)
     same_design <- all(vapply(c('X', 'Xzi', 'y', 'offset', 'group', 'family'),
       function(nm) identical(args[[nm]], terminal_args[[nm]]), TRUE))
     # Static multipliers, native columns and kernel dependencies must agree
     # across independently constructed native fits, without comparing envs.
     fields<-c('multipliers','dependence','groups','covariates','columns',
       'kernel_frame_columns','n','p','q','assignment','contrasts',
       'expression_programs','expression_frame_columns','expression_nodes')
     same_design<-same_design && identical(is.null(args$response_design),
       is.null(terminal_args$response_design)) &&
       identical(args$response_design[fields],terminal_args$response_design[fields])
     delta <- if (eligible_native(terminal)) abs(value + as.numeric(stats::logLik(terminal))) else Inf
     ok <- opt$convergence == 0L && isTRUE(sd$pdHess) &&
       eligible_native(terminal) && same_design && is.finite(delta) && delta < 1e-5
     list(eligible = ok, log_scale = kp, objective = value, native_delta = delta,
          design_unchanged = same_design, optim = opt, joint_sdreport = sd,
          native = terminal)
   }, warning = function(w) {
     warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
   }), error = function(e) list(eligible = FALSE, error = conditionMessage(e)))
   outcomes[[i]]$warnings <- warnings
   outcomes[[i]]$start <- starts[i, ]
 }
 good <- which(vapply(outcomes, function(x) isTRUE(x$eligible), TRUE))
 selected <- if (length(good)) good[which.min(vapply(outcomes[good], function(x) x$objective, 0))] else NA_integer_
 # Keep failures and every independent start in the return value, even if all fail.
 list(outcomes = outcomes, selected = selected,
      model = if (!is.na(selected)) outcomes[[selected]]$native else NULL,
      uncertainty = 'Native coefficient/prediction covariance conditional on selected scales; scale uncertainty uncalibrated.')
}

