#include <Rcpp.h>

#include <cmath>
#include <utility>
#include <vector>

using namespace Rcpp;

// Offset into the focal window plus its kernel weight.
struct WeightedOffset {
  int dr;
  int dc;
  double w;
};

// Compiled masked-window focal pass for the `sa` (average roughness) surface
// metric, in both unweighted and kernel-weighted form. Average roughness is the
// (weighted) mean absolute deviation of the neighborhood values from the
// (weighted) neighborhood mean. Unlike `sq` (a variance, which expands into
// neighborhood sums and is therefore computed by FFT convolution), `sa` has no
// closed-form FFT decomposition because the absolute value couples every cell
// to the output location through that window's own mean. This routine performs
// an exact two-pass focal calculation in compiled code, which is much faster
// than calling an R closure per window via terra::focal().
//
// `values` is the raster as a matrix with values(i, j) = the value at raster
// row i, column j (i.e., terra::as.matrix(raster, wide = TRUE)). `kernel`
// supplies both the window and the per-cell weights: nonzero (and non-NA)
// entries mark in-window cells and their value is the weight. A 0/1 circular
// mask (e.g., from .landscape_circle_kernel()) therefore reproduces the
// unweighted average roughness exactly, while a distance-decay weight matrix
// gives the kernel-weighted version. The kernel is assumed to have odd
// dimensions so the focal cell sits at the center. The result is returned in
// the same orientation as `values`.
//
// NA and non-finite cells are dropped from each window (matching na.rm = TRUE
// behavior); the focal cell itself may be NA and the metric is still computed
// from its finite neighbors. A window with no positive total weight returns NA.
// [[Rcpp::export]]
NumericMatrix surface_sa_focal_cpp(NumericMatrix values, NumericMatrix kernel) {
  const int nr = values.nrow();
  const int nc = values.ncol();
  const int kr = kernel.nrow();
  const int kc = kernel.ncol();
  const int r_off = kr / 2;
  const int c_off = kc / 2;

  // Precompute the in-window (row, col) offsets and weights relative to the
  // focal cell so the inner loops only visit cells inside the mask.
  std::vector<WeightedOffset> offsets;
  offsets.reserve(static_cast<std::size_t>(kr) * kc);
  for (int i = 0; i < kr; ++i) {
    for (int j = 0; j < kc; ++j) {
      double w = kernel(i, j);
      if (!NumericVector::is_na(w) && w != 0.0) {
        WeightedOffset o;
        o.dr = i - r_off;
        o.dc = j - c_off;
        o.w = w;
        offsets.push_back(o);
      }
    }
  }

  NumericMatrix out(nr, nc);

  for (int r = 0; r < nr; ++r) {
    for (int c = 0; c < nc; ++c) {
      // First pass: weighted neighborhood mean over finite cells.
      double sum_w = 0.0;
      double sum_wv = 0.0;
      for (std::size_t k = 0; k < offsets.size(); ++k) {
        int rr = r + offsets[k].dr;
        int cc = c + offsets[k].dc;
        if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) {
          continue;
        }
        double v = values(rr, cc);
        if (NumericVector::is_na(v) || !std::isfinite(v)) {
          continue;
        }
        sum_w += offsets[k].w;
        sum_wv += offsets[k].w * v;
      }

      if (sum_w <= 0.0) {
        out(r, c) = NA_REAL;
        continue;
      }

      double mean = sum_wv / sum_w;

      // Second pass: weighted mean absolute deviation from that mean.
      double ad_sum = 0.0;
      for (std::size_t k = 0; k < offsets.size(); ++k) {
        int rr = r + offsets[k].dr;
        int cc = c + offsets[k].dc;
        if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) {
          continue;
        }
        double v = values(rr, cc);
        if (NumericVector::is_na(v) || !std::isfinite(v)) {
          continue;
        }
        ad_sum += offsets[k].w * std::abs(v - mean);
      }

      out(r, c) = ad_sum / sum_w;
    }
  }

  return out;
}
