#include <Rcpp.h>
#include <algorithm>
#include "kernel_derivatives.h"

// Validate without allocating full-sized logical comparison matrices in R.
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalVector msr_kernel_validate_cpp(Rcpp::NumericMatrix v,
                                            Rcpp::NumericMatrix c,
                                            Rcpp::IntegerVector rows) {
  int n = v.nrow();
  std::vector<unsigned char> support(n, 0);
  const double *vp = v.begin(), *cp = c.begin();
  for (R_xlen_t k = 0; k < v.size(); ++k) {
    if (!R_finite(vp[k]) || !R_finite(cp[k]) || cp[k] < 0 ||
        (cp[k] == 0 && vp[k] != 0))
      return Rcpp::LogicalVector::create(false, false);
    if (cp[k] > 0) support[k % n] = 1;
  }
  for (int i : rows) if (!support[i - 1])
    return Rcpp::LogicalVector::create(true, false);
  return Rcpp::LogicalVector::create(true, true);
}

// R input validation and row selection happen once in the prepared engine.
Rcpp::List msr_kernel_eval_impl(Rcpp::List V, Rcpp::List C,
                              Rcpp::NumericVector distance,
                              Rcpp::NumericVector log_scale,
                              Rcpp::NumericVector log_shape,
                              int kernel, bool standardize, bool columnwise) {
  int p = V.size(), n = Rcpp::as<Rcpp::NumericMatrix>(V[0]).nrow();
  Rcpp::NumericMatrix raw(n, p), value(n, p), ds(n, p), dq(n, p);
  Rcpp::NumericVector center(p), spread(p);
  for (int j = 0; j < p; ++j) {
    Rcpp::NumericMatrix v = V[j], c = C[j];
    const double *vp = v.begin(), *cp = c.begin();
    auto out = columnwise ? msr::kernel_column<double>(
      n, distance.size(), [&](int i, int k) { return vp[i + n * k]; },
      [&](int i, int k) { return cp[i + n * k]; },
      [&](int k) { return distance[k]; },
      log_scale[j], log_shape[j], kernel, standardize, true, true) :
      msr::kernel_column<double>(
      n, distance.size(), [&](int i, int k) { return v(i, k); },
      [&](int i, int k) { return c(i, k); },
      [&](int k) { return distance[k]; },
      log_scale[j], log_shape[j], kernel, standardize, true, false);
    center[j] = out.center; spread[j] = out.spread;
    if (standardize &&
        *std::min_element(out.raw.begin(), out.raw.end()) ==
        *std::max_element(out.raw.begin(), out.raw.end()))
      Rcpp::stop("Undefined kernel value or derivative: constant predictor.");
    for (int i = 0; i < n; ++i) {
      if (!R_finite(out.value[i]) || !R_finite(out.scale[i]) ||
          !R_finite(out.shape[i]))
        Rcpp::stop("Undefined kernel value or derivative: check support, scale bounds, and predictor variance.");
      raw(i,j) = out.raw[i]; value(i,j) = out.value[i];
      ds(i,j) = out.scale[i]; dq(i,j) = out.shape[i];
    }
  }
  return Rcpp::List::create(Rcpp::_["raw"] = raw, Rcpp::_["value"] = value,
    Rcpp::_["d_scale"] = ds, Rcpp::_["d_shape"] = dq,
    Rcpp::_["center"] = center, Rcpp::_["spread"] = spread);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List msr_kernel_eval_cpp(Rcpp::List V, Rcpp::List C,
                              Rcpp::NumericVector distance,
                              Rcpp::NumericVector log_scale,
                              Rcpp::NumericVector log_shape,
                              int kernel, bool standardize) {
  return msr_kernel_eval_impl(V, C, distance, log_scale, log_shape,
                             kernel, standardize, false);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List msr_kernel_stream_cpp(Rcpp::List V, Rcpp::List C,
                                Rcpp::NumericVector distance,
                                Rcpp::NumericVector log_scale,
                                Rcpp::NumericVector log_shape,
                                int kernel, bool standardize) {
  return msr_kernel_eval_impl(V, C, distance, log_scale, log_shape,
                             kernel, standardize, true);
}
