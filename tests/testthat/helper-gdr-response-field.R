# A changed, explicit three-state numeric layout for fixed-point conformance.
# No fitted, simulated or sampled outcome is selected here.
gdr_response_field_fixture <- function(K = 2L) {
  f <- gdr_gaussian_fixture(TRUE, TRUE, K, unused = TRUE)
  index <- as.integer(c(1, 2, 3, 1, 2, 3, 1, 2))
  f$data$block <- index
  f$captured$data$Z_lambda <- Matrix::sparseMatrix(i = seq_along(index), j = index,
    x = rep(1, length(index)), dims = c(8L, 3L))
  locations <- rbind(c(430000, 4400000), c(430040, 4400010), c(430020, 4400070))
  f$locations <- locations
  f$coords <- locations[index, , drop = FALSE]
  rownames(f$coords) <- f$data$row_ids$survey
  f$crs <- 32617
  f
}
