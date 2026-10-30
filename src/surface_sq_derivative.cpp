#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Raw centered population SQ and its log-scale derivative. Reuse two bounded
// work vectors across sites rather than allocating filtered/deviation arrays.
// [[Rcpp::export]]
NumericMatrix msr_surface_sq_derivative_cpp(List buffers, double sigma,
                                           std::string kernel) {
  if (!std::isfinite(sigma) || sigma <= 0 ||
      (kernel != "gaussian" && kernel != "exp"))
    stop("Invalid smooth surface scale or kernel.");
  NumericMatrix out(buffers.size(), 2);
  std::vector<double> weights, scores;
  for (R_xlen_t i = 0; i < buffers.size(); ++i) {
    if (i % 256 == 0) checkUserInterrupt();
    List b = buffers[i];
    NumericVector d = b["d"], v = b["v"];
    if (d.size() != v.size()) stop("Invalid surface buffer.");
    weights.resize(d.size()); scores.resize(d.size());
    long double sw = 0, sx = 0;
    double anchor = 0;
    bool anchored = false;
    for (R_xlen_t j = 0; j < d.size(); ++j) {
      if (!std::isfinite(d[j]) || !std::isfinite(v[j]))
        stop("Surface buffers must contain finite cells.");
      double r = d[j] / sigma;
      double score = kernel == "gaussian" ? r * r : r;
      double w = std::exp(kernel == "gaussian" ? -score / 2 : -score);
      weights[j] = std::isfinite(w) && w > 0 ? w : 0;
      scores[j] = score;
      if (weights[j] > 0) {
        if (!anchored) { anchor = v[j]; anchored = true; }
        sw += w; sx += w * (v[j] - anchor);
      }
    }
    double total = static_cast<double>(sw);
    if (!std::isfinite(total) || total <= 0)
      stop("No effective surface support at this scale.");
    double mean = static_cast<double>(sx) / total;
    long double ss = 0;
    for (R_xlen_t j = 0; j < d.size(); ++j) if (weights[j] > 0) {
      double dev = (v[j] - anchor) - mean;
      ss += weights[j] * (dev * dev);
    }
    double variance = static_cast<double>(ss) / total;
    double sq = std::sqrt(variance);
    long double gradient = 0;
    if (variance != 0) {
      for (R_xlen_t j = 0; j < d.size(); ++j) if (weights[j] > 0) {
        double dev = (v[j] - anchor) - mean;
        gradient += weights[j] * scores[j] * (dev * dev - variance);
      }
    }
    out(i, 0) = sq;
    out(i, 1) = variance == 0 ? 0 : static_cast<double>(gradient) / (2 * total * sq);
    if (!std::isfinite(out(i, 0)) || !std::isfinite(out(i, 1)))
      stop("Nonfinite surface values or derivatives.");
  }
  return out;
}
