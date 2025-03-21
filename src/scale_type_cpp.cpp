#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector scale_type_cpp(NumericVector d,
                             std::string kernel = "gaussian",
                             Nullable<NumericVector> sigma_ = R_NilValue,
                             Nullable<NumericVector> shape_ = R_NilValue,
                             Nullable<NumericMatrix> r_stack_df = R_NilValue,
                             Nullable<int> output = R_NilValue) {

  int n = d.size();

  // Check constraints on r_stack_df and sigma
  if (sigma_.isNotNull() && r_stack_df.isNull() && output.isNull()) {
    stop("r_stack_df cannot be NULL when using sigma unless output is provided.");
  }

  NumericMatrix r_stack;
  int ncol = 0;
  if (r_stack_df.isNotNull()) {
    r_stack = NumericMatrix(r_stack_df);
    ncol = r_stack.ncol();
  }

  // Ensure sigma is provided and matches ncol
  if (sigma_.isNull()) {
    stop("sigma must be provided and must have the same length as the number of columns in r_stack_df.");
  }

  NumericVector sigma = NumericVector(sigma_);
  if (sigma.size() == 1 && ncol > 1) {
    sigma = rep(sigma[0], ncol); // Expand sigma to match ncol
  } else if (r_stack_df.isNotNull() && sigma.size() != ncol) {
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

  NumericMatrix w0(n, ncol);

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

  // Normalize weights for each column
  for (int j = 0; j < ncol; j++) {
    double sum_w = 0.0;
    for (int i = 0; i < n; i++) {
      sum_w += w0(i, j);
    }
    if (sum_w == 0) {
      warning("Sum of weights for column %d is zero. Normalization skipped.", j);
    } else {
      for (int i = 0; i < n; i++) {
        w0(i, j) /= sum_w;
      }
    }
  }

  // If output is not NULL, return the first column of normalized weights
  if (output.isNotNull()) {
    return w0(_, 0);
  }

  // Compute weighted raster values
  NumericVector result(ncol);
  for (int j = 0; j < ncol; j++) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
      sum += r_stack(i, j) * w0(i, j);
    }
    result[j] = sum;
  }

  return result;
}
