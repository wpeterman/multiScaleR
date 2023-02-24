#' @title Scale Distance
#' @description Function to estimate the effective distance encompassing a specified probability density of the kernel density function
#' @param x \code{\link[multiScaleR]{multiScaleR}} object
#' @param prob Density probability cutoff for calculating distance, Default: 0.9
#' @return Distance
#' @details None
#' @examples
#' kernel_dist(x)
#'
#' @rdname kernel_dist
#' @export
kernel_dist <- function(model,
                        prob = 0.9){
  if(class(model) != 'multiScaleR'){
    stop("Provide a fitted `multiScaleR` model object")
  }

  d <- seq(1, 1e6, length.out = 1e6)
  ci_ <- summary(x)$opt_scale
  dist_list <- vector('list', nrow(ci_))

  for(i in 1:nrow(ci_)){

    wt_mn <- scale_type(d = d,
                        kernel = model$kernel_inputs$kernel,
                        sigma = ci_[[1]],
                        shape = model$shape_est[[1]],
                        output = 'wts')


    wt_l <- scale_type(d = d,
                       kernel = model$kernel_inputs$kernel,
                       sigma = ci_[[3]],
                       shape = model$shape_est[[1]],
                       output = 'wts')

    wt_u <- scale_type(d = d,
                       kernel = model$kernel_inputs$kernel,
                       sigma = ci_[[4]],
                       shape = model$shape_est[[1]],
                       output = 'wts')

    scale_mn <- DescTools::Quantile(d,
                                    weights = wt_mn,
                                    probs = prob)
    scale_l <- DescTools::Quantile(d,
                                   weights = wt_l,
                                   probs = prob)
    scale_u <- DescTools::Quantile(d,
                                   weights = wt_u,
                                   probs = prob)
    dist_list[[i]] <- data.frame(mn = scale_mn,
                                 l = scale_l,
                                 u = scale_u)
  }
  dist_out <- do.call(rbind, dist_list)
  rownames(dist_out) <- rownames(ci_)
  colnames(dist_out) <- colnames(ci_)[c(1,3,4)]
}
