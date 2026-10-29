#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
// Centered weighted population spread, matching the R weighted SQ definition.
// One shared weight buffer; no per-column filtered or deviation vectors.
// [[Rcpp::export]]
NumericVector msr_surface_weighted_sq_cpp(NumericVector d, NumericMatrix values,
                                         double sigma, double shape, std::string kernel) {
  if (values.nrow() != d.size()) stop("Distance and value rows must agree.");
  if (!std::isfinite(sigma) || sigma <= 0) stop("sigma must be finite and positive.");
  if (kernel != "gaussian" && kernel != "exp" && kernel != "fixed" && kernel != "expow")
    stop("Invalid kernel type.");
  NumericVector weights(d.size()), out(values.ncol(), NA_REAL);
  for (R_xlen_t i=0; i<d.size(); ++i) {
    double x=d[i];
    weights[i] = kernel == "gaussian" ? std::exp(-(x*x)/(2*sigma*sigma)) :
      kernel == "exp" ? std::exp(-x/sigma) :
      kernel == "fixed" ? (std::isnan(x) ? NA_REAL : double(x<sigma)) :
      std::exp(-std::pow(x/sigma,shape));
  }
  for (int j=0; j<values.ncol(); ++j) {
    long double sw=0, sx=0;
    for (R_xlen_t i=0; i<d.size(); ++i) {
      double v=values(i,j), w=weights[i];
      if (std::isfinite(v) && std::isfinite(w) && w>0) { sw+=w; sx+=w*v; }
    }
    double s0=static_cast<double>(sw);
    if (!std::isfinite(s0) || s0<=0) continue;
    double mu=static_cast<double>(sx)/s0;
    long double ss=0;
    for (R_xlen_t i=0; i<d.size(); ++i) {
      double v=values(i,j), w=weights[i];
      if (std::isfinite(v) && std::isfinite(w) && w>0) {
        double dev=v-mu;
        ss+=w*(dev*dev);
      }
    }
    double m2=static_cast<double>(ss)/s0;
    if (std::isnan(m2)) stop("missing value where TRUE/FALSE needed");
    out[j]=m2<=0 ? 0 : std::sqrt(m2);
  }
  return out;
}


