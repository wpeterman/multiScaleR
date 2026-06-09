#include <Rcpp.h>

#include <cmath>
#include <utility>
#include <vector>

using namespace Rcpp;

// Compiled masked-window focal pass for the `sa` (average roughness) surface
// metric. Average roughness is the mean absolute deviation of the neighborhood
// values from the neighborhood mean. Unlike `sq` (a variance, which expands
// into neighborhood sums and is therefore computed by FFT convolution), `sa`
// has no closed-form FFT decomposition because the absolute value couples every
// cell to the output location through that window's own mean. This routine
// therefore performs an exact two-pass focal calculation in compiled code,
// which is much faster than calling an R closure per window via terra::focal().
//
// `values` is the raster as a matrix with values(i, j) = the value at raster
// row i, column j (i.e., terra::as.matrix(raster, wide = TRUE)). `kernel` is
// the circular neighborhood mask, with nonzero entries marking in-window cells
// (e.g., from .landscape_circle_kernel()); it is assumed to have odd
// dimensions so the focal cell sits at the center. The result is returned in
// the same orientation as `values`.
//
// NA and non-finite cells are dropped from each window (matching na.rm = TRUE
// behavior); the focal cell itself may be NA and the metric is still computed
// from its finite neighbors. A window with no finite cells returns NA.
// [[Rcpp::export]]
NumericMatrix surface_sa_focal_cpp(NumericMatrix values, NumericMatrix kernel) {
  const int nr = values.nrow();
  const int nc = values.ncol();
  const int kr = kernel.nrow();
  const int kc = kernel.ncol();
  const int r_off = kr / 2;
  const int c_off = kc / 2;

  // Precompute the in-window (row, col) offsets relative to the focal cell so
  // the inner loops only visit cells inside the circular mask.
  std::vector<std::pair<int, int> > offsets;
  offsets.reserve(static_cast<std::size_t>(kr) * kc);
  for (int i = 0; i < kr; ++i) {
    for (int j = 0; j < kc; ++j) {
      double w = kernel(i, j);
      if (!NumericVector::is_na(w) && w != 0.0) {
        offsets.push_back(std::make_pair(i - r_off, j - c_off));
      }
    }
  }

  NumericMatrix out(nr, nc);

  for (int r = 0; r < nr; ++r) {
    for (int c = 0; c < nc; ++c) {
      // First pass: neighborhood mean over finite cells.
      double sum = 0.0;
      int n = 0;
      for (std::size_t k = 0; k < offsets.size(); ++k) {
        int rr = r + offsets[k].first;
        int cc = c + offsets[k].second;
        if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) {
          continue;
        }
        double v = values(rr, cc);
        if (NumericVector::is_na(v) || !std::isfinite(v)) {
          continue;
        }
        sum += v;
        n += 1;
      }

      if (n == 0) {
        out(r, c) = NA_REAL;
        continue;
      }

      double mean = sum / n;

      // Second pass: mean absolute deviation from that mean.
      double ad_sum = 0.0;
      for (std::size_t k = 0; k < offsets.size(); ++k) {
        int rr = r + offsets[k].first;
        int cc = c + offsets[k].second;
        if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) {
          continue;
        }
        double v = values(rr, cc);
        if (NumericVector::is_na(v) || !std::isfinite(v)) {
          continue;
        }
        ad_sum += std::abs(v - mean);
      }

      out(r, c) = ad_sum / n;
    }
  }

  return out;
}
