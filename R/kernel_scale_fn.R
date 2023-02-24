#' @title Kernel scaling function
#' @description Function for internal use with optim
#' @param par list of parameters
#' @param kernel_inputs object created from `kernel_prep`
#' @param fitted_mod fitted model object
#' @param mod_return Default: NULL
#' @return Either estimated parameters or the fitted model using provided parameters
#' @details For internal use
#' @rdname kernel_scale_fn
#' @export
#' @keywords internal
#' @importFrom insight get_data
#' @importFrom stats formula logLik

kernel_scale_fn <- function(par,
                            kernel_inputs,
                            fitted_mod,
                            mod_return = NULL){


  # D <- kernel_inputs$D
  d_list <- kernel_inputs$d_list

  cov_df <- kernel_inputs$raw_cov
  kernel <- kernel_inputs$kernel
  n_ind <- length(d_list)
  mod <- fitted_mod

  if(any(class(mod) == 'gls')){
    mod_class <- 'gls'
    dat <- nlme::getData(mod)
    covs <- all.vars(formula(mod)[-2])
    covs <- covs[which(covs %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(covs)
  } else if(any(grepl("^unmarked", class(mod)))) {
    mod_class <- 'unmarked'
    dat <- mod@data@siteCovs
    covs <- all.vars(mod@formula)
    covs <- covs[which(covs %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(covs)
  } else if(any(class(mod) == 'glm')) {
    mod_class <- 'glm'
    dat <- fitted_mod$data
    covs <- all.vars(formula(mod)[-2])
    covs <- covs[which(covs %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(covs)
  } else {
    dat <- insight::get_data(mod)
    covs <- all.vars(formula(mod)[-2])
    covs <- covs[which(covs %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(covs)
  }

  if(is.null(dat)){
    stop("Data from original model not saved to data frame. Try using `glm`.\n UDPATE FUNCTION TO GENERALIZE!!!")
  }

  scale_type_opt <- function(d,
                             kernel = 'expow',
                             sigma,
                             shape = NULL,
                             r_stack.df = NULL,
                             output = NULL) {
    if(kernel == 'exp'){
      # w0 <- ((1/(2*pi*sigma^2)) * exp(-(d / sigma))) * as.numeric(d!=0)  ## Neg. Exponential
      w0 <- ((1/(2*pi*sigma^2)) * exp(-outer(d, sigma, "/")))

    } else if(kernel == 'fixed'){
      # w0 <- as.numeric(d != 0 & outer(d, sigma, "<"))  ## Fixed
      w0 <- as.numeric(outer(d, sigma, "<"))  ## Fixed

    } else if(kernel == 'gaussian'){
      # w0_ <- ((1/(pi*sigma[1]^2)) * exp(-(d^2 / (2*sigma[1]^2)))) #* as.numeric(d!=0)  ## Original
      w0 <- ((1/(pi*sigma^2)) * exp(-outer(d^2, (2*sigma^2), FUN = "/"))) #* as.numeric(d!=0)  ## Original

    } else {
      # w0 <- ((shape/((2*pi*sigma^2)*gamma(2/shape))) * exp(-(d^shape / (sigma^shape)))) * as.numeric(d!=0)  ## Exponential power
      w0 <- ((shape/((2*pi*sigma^2)*gamma(2/shape)))  * exp(-(outer(d,shape,"^") / (sigma^shape)))) ## Exponential power
    }

    w <- apply(w0, 2, function(x) x/sum(x))

    if(!is.null(output)){
      return(w)
    } else {
      # if(is.null(cov_subset)){
      #   sum(r_stack.df[,layer] * w)
      #
      # } else {
      #   sum(r_stack.df[cov_subset, layer+2] * w)
      #
      # }

      # colSums(r_stack.df[r_stack.df$ID == ID, 2:(nlayer+1)] * w) ## terra solution
      colSums(r_stack.df * w)  ## exact_extract solution
    }
  }

  cov.w <- vector('list', n_ind)

  # n_cores <- 4
  #create the cluster
  # my.cluster <- parallel::makeCluster(
  #   n_cores,
  #   type = "PSOCK"
  # )

  #register it to be used by %dopar%
  # doParallel::registerDoParallel(cl = my.cluster)

  #check if it is registered (optional)
  # foreach::getDoParRegistered()

  # cov.w_ <- foreach(i = 1:n_ind) %dopar% {
  for(i in 1:n_ind){
    sigma <- par[1:n_covs]
    if(kernel == 'expow'){
      shape <- par[(n_covs + 1):(n_covs * 2)]
    } else {
      shape <- NULL
    }

    cov.w[[i]] <-
      scale_type_opt(d_list[[i]],
                     kernel = kernel,
                     sigma = sigma,
                     shape = shape,
                     r_stack.df = cov_df[[i]][,covs])

  } ## End for loop
  # }
  # parallel::stopCluster(cl = my.cluster)


  # cov.w[[i]] <- apply(D, 1, function(x) {
  #   w0 <- exp(-x^2 / (2*sigma^2)) ## x > sigma --> fixed buffer solution?
  #   # w0[w0==1] <- 0
  #   w0[w0==1] <- 0
  #   w <- w0/sum(w0)
  #   sum(cov_df[,i] * w)
  # }

  df <- data.frame(do.call(rbind, cov.w))
  colnames(df) <- covs
  if(mod_class == 'unmarked'){
    umf <- mod@data
    scl_df <- scale(df)
    # scl_df <- (df)
    umf@siteCovs[,covs] <- scl_df
    mod <- update(mod, data = umf)
  } else {
    scl_df <- scale(df)
    # scl_df <- (df)
    dat[,covs] <- scl_df
    mod <- update(mod, data = dat)
  }


  if(is.null(mod_return)){
    obj <- logLik(mod)[1] * -1
  } else {
    obj <- list(mod = mod,
                scl_params = list(mean = attr(scl_df, "scaled:center"),
                                  sd = attr(scl_df, "scaled:scale")))
    # obj <- list(mod = mod,
    #             scl_params = NULL)
  }
  return(obj)

}
