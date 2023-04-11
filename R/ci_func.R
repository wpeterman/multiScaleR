#' @title CI Function
#' @description Internal function for calculating confidence intervals
#' @param x Estimates
#' @param df Degrees of Freedom
#' @param min_D minimum distance, Default: NULL
#' @param names Names of fitted variables, Default: NULL
#' @param as_dist Calculate distance confidence interval; Default = FALSE
#' @param ... Not used
#' @return data frame with calculated confidence intervals

#' @rdname ci_func
#' @keywords internal
#'
ci_func <- function(x,
                    df,
                    min_D = NULL,
                    names = NULL,
                    as_dist = FALSE,
                    ...){

  out <- vector('list', nrow(x))
  for(i in 1:nrow(x)){
    if(is.nan(x[i,2]) | (x[i,2] == 'Inf')){
      ci <- quantile(rnorm(10000, 10, 0.1), c(0.025, 0.975))
      ci[ci > 0] <- NaN
      out[[i]] <- ci
    } else {
      ci <- c("2.5%" = x[i,1] - qt(0.975, df = df) * x[i,2],
              "97.5%" = x[i,1] + qt(0.975, df = df) * x[i,2])
      # ci <- quantile(rnorm(10000, x[i,1], x[i,2]), c(0.025, 0.975)) ## Original
      if(!is.null(min_D)){
        ci[ci < min_D] <- min_D
        out[[i]] <- ci
      } else {
        # ci[ci < 0] <- 0
        out[[i]] <- ci
      }
    }
  }

  out_df <- cbind(x, do.call(rbind, out))
  rownames(out_df) <- names
  return(out_df)
}
