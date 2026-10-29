#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Sorted finite values. Shifted Welford updates avoid cancellation from large
// absolute elevations and compute every prefix in one bounded linear pass.
// [[Rcpp::export]]
NumericVector msr_surface_sq_prefix_cpp(NumericVector value) {
  R_xlen_t n = value.size();
  NumericVector out(n, NA_REAL);
  if (!n) return out;
  long double anchor = value[0], mean = 0, m2 = 0;
  for (R_xlen_t j = 0; j < n; ++j) {
    if (!std::isfinite(value[j])) stop("Prefix values must be finite.");
    long double x = static_cast<long double>(value[j]) - anchor;
    long double delta = x - mean;
    mean += delta / (j + 1);
    m2 += delta * (x - mean);
    if (j) out[j] = static_cast<double>(m2 / j);
  }
  return out;
}
