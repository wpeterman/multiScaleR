#' @title Print multiScaleR data
#' @description Print summary method for multiScaleR data objects
#' @param x Object of class `multiScaleR_data`
#' @param ... Not used
#' @return Summary of data
#' @usage
#' print(x)
#' @rdname print.multiScaleR_data
#' @export

print.multiScaleR_data <- function(x, ...){
  cat("\nThere are ")
  cat(paste0(nrow(x$kernel_dat)," observations at ", ncol(x$kernel_dat), ' spatial covariate(s): \n'))
  cat(colnames(x$kernel_dat))
  cat("\n\nThe specified kernel is:\n")
  cat(x$kernel)
  # cat("\n\nSparse Matrix contains: ")
  # cat(paste0(length(x$D@x), ' elements\n'))
  cat("\n\nNumber of elements: \n")
  cat(paste0(length(x$d_list[[1]])))
  cat("\nMinimum Distance:\n")
  cat(x$min_D)
  cat("\nMaximum Distance:\n")
  cat(x$max_D)
  cat("\nUnit Conversion:\n")
  cat(x$unit_conv)
}
