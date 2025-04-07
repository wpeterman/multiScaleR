#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector scale_type_sparse(NumericVector d,
                                std::string kernel = "gaussian",
                                Nullable<NumericVector> sigma_ = R_NilValue,
                                Nullable<NumericVector> shape_ = R_NilValue,
                                Nullable<RObject> r_stack_df = R_NilValue,
                                Nullable<std::string> output = R_NilValue) {

  int n = d.size();  // Number of distances
  sp_mat r_stack;    // Sparse matrix to store raster stack
  int ncol = 0;      // Number of columns in raster stack

  // Check constraints on r_stack_df and sigma
  if (sigma_.isNotNull() && r_stack_df.isNull() && output.isNull()) {
    stop("r_stack_df cannot be NULL when using sigma unless output is provided.");
  }

  if (r_stack_df.isNotNull()) {
    r_stack = as<sp_mat>(r_stack_df); // Convert R object to sparse matrix
    ncol = r_stack.n_cols;  // Number of columns (layers)
  }

  // Ensure sigma is provided and matches ncol
  if (sigma_.isNull()) {
    stop("sigma must be provided and must have the same length as the number of columns in r_stack_df.");
  }

  NumericVector sigma = NumericVector(sigma_);
  if (sigma.size() == 1 && ncol > 1) {
    sigma = rep(sigma[0], ncol); // Expand sigma to match ncol if it's a scalar
  } else if (sigma.size() != ncol && ncol != 0) {
    stop("sigma must have the same length as the number of columns in r_stack_df.");
  }

  // Handle shape: required only if kernel == "expow"
  NumericVector shape(ncol, 1.0);
  if (kernel == "expow") {
    if (shape_.isNull()) {
      stop("shape must be provided when kernel is 'expow'.");
    }
    NumericVector shape_input = NumericVector(shape_);
    if (shape_input.size() == 1 && ncol > 1) {
      shape_input = rep(shape_input[0], ncol); // Expand shape if needed
    } else if (shape_input.size() != ncol) {
      stop("shape must have the same length as the number of columns in r_stack_df.");
    }
    shape = shape_input;
  }

  mat w0(n, ncol, fill::zeros);  // Matrix to hold weights

  // Apply the kernel functions column-wise using sigma[j] and shape[j] where needed
  for (int j = 0; j < ncol; j++) {
    for (int i = 0; i < n; i++) {
      if (kernel == "exp") {
        w0(i, j) = exp(-d[i] / sigma[j]) * (1 / (2 * M_PI * sigma[j] * sigma[j]));
      } else if (kernel == "fixed") {
        w0(i, j) = (d[i] < sigma[j]) ? 1.0 : 0.0;
      } else if (kernel == "gaussian") {
        w0(i, j) = exp(- (d[i] * d[i]) / (2 * sigma[j] * sigma[j])) * (1 / (M_PI * sigma[j] * sigma[j]));
      } else if (kernel == "expow") {
        w0(i, j) = exp(-pow(d[i], shape[j]) / pow(sigma[j], shape[j])) *
          (shape[j] / ((2 * M_PI * sigma[j] * sigma[j]) * std::tgamma(2.0 / shape[j])));
      } else {
        stop("Invalid kernel type.");
      }
    }
  }

  // Normalize weights for each column using sparse matrix-specific operations
  for (int j = 0; j < ncol; j++) {
    double sum_w = sum(w0.col(j));
    if (sum_w == 0) {
      warning("Sum of weights for column %d is zero. Normalization skipped.", j);
    } else {
      w0.col(j) /= sum_w;
    }
  }

  // Handle output selection
  if (output.isNotNull()) {
    return wrap(w0);  // Return weights matrix
  }


  // Compute weighted raster values using sparse matrix operations
  // Use sparse matrix multiplication instead of dense matrix operations
  NumericVector result(ncol);
  for (int j = 0; j < ncol; j++) {
    result[j] = dot(r_stack.col(j), w0.col(j)); // Efficient sparse matrix-vector multiplication
  }

  return result;
}
