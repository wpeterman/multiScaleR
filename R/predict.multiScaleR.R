#' @title multiScaleR prediction
#' @description Function to make spatial predictions of scale-optimized model
#' @param object Fitted model `multiScaleR` object
#' @param link PARAM_DESCRIPTION, Default: NULL
#' @param raster_stack PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details Function not yet completed or operational
#' @examples
#' ## Not Run
#'
#' @rdname predict.multiScaleR
#' @keywords internal
#' @export
predict.multiScaleR <- function(object,
                                link = NULL,
                                raster_stack,
                                ...){
  if(is.null(link)){
    stop("Specify the link function used in fitting your model.")
  }

  ## Scale rasters



  if(any(class(fitted_mod) == 'gls')){
    mod_class <- 'gls'
    link <- 'identity'
    mod_vars <- all.vars(formula(fitted_mod)[-2])
    r_vars <- mod_vars[which(mod_vars %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(r_vars)
  } else if(any(grepl("^unmarked", class(fitted_mod)))) {
    mod_class <- 'unmarked'
    link <- 'identity'
    mod_vars <- all.vars(formula(fitted_mod@formula))
    r_vars <- mod_vars[which(mod_vars %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(r_vars)
    fitType <- fitted_mod@fitType
  } else {
    mod_class <- 'glm'
    link <- 'identity'
    mod_vars <- all.vars(formula(fitted_mod)[-2])
    r_vars <- mod_vars[which(mod_vars %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(r_vars)
  }

  link <- x$opt_mod$family
}
