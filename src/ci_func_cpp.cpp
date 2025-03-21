#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix ci_func_cpp(NumericMatrix x,
                          int df,
                          Nullable<double> min_D = R_NilValue,
                          Nullable<CharacterVector> names = R_NilValue) {

  int n = x.nrow();
  NumericMatrix out_df(n, x.ncol() + 2); // Original matrix + 2 columns for CI

  // Copy input matrix to output
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < x.ncol(); j++) {
      out_df(i, j) = x(i, j);
    }
  }

  for (int i = 0; i < n; i++) {
    if (NumericVector::is_na(x(i, 1)) || !R_finite(x(i, 1))) {
      NumericVector ci = NumericVector::create(NA_REAL, NA_REAL);
      out_df(i, x.ncol()) = ci[0];
      out_df(i, x.ncol() + 1) = ci[1];
    } else {
      double lower = x(i, 0) - R::qt(0.975, df, 1, 0) * x(i, 1);
      double upper = x(i, 0) + R::qt(0.975, df, 1, 0) * x(i, 1);

      if (min_D.isNotNull()) {
        double min_val = as<double>(min_D);
        if (lower < min_val) lower = min_val;
        if (upper < min_val) upper = min_val;
      }

      out_df(i, x.ncol()) = lower;
      out_df(i, x.ncol() + 1) = upper;
    }
  }

  // Assign row names if provided
  if (names.isNotNull()) {
    rownames(out_df) = as<CharacterVector>(names);
  }

  return out_df;
}

