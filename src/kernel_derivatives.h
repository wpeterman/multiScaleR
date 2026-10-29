#ifndef MULTISCALER_KERNEL_DERIVATIVES_H
#define MULTISCALER_KERNEL_DERIVATIVES_H

#include <cmath>
#include <vector>

// Shared arithmetic for ordinary doubles and TMB automatic-differentiation
// types. V and C contain finite-cell sums and counts in column-major order.
// Kernel codes: 0 Gaussian, 1 exponential, 2 exponential-power.
// Each column has its own scale (and, for code 2, shape). All derivatives
// are with respect to log parameters. Hard cutoffs are intentionally excluded.
namespace msr {
template<class Type>
struct KernelColumn {
  std::vector<Type> raw, value, scale, shape;
  Type center, spread;
};

template<class Type, class Value, class Count, class Distance>
KernelColumn<Type> kernel_column(int n, int nb, Value V, Count C, Distance d,
                                Type log_scale, Type log_shape, int kernel,
                                bool standardize, bool derivatives,
                                bool columnwise = false) {
  using std::exp;
  using std::log;
  using std::sqrt;
  KernelColumn<Type> out;
  out.raw.resize(n); out.value.resize(n);
  out.scale.resize(n); out.shape.resize(n);
  std::vector<Type> w(nb), hs(nb), hq(nb);
  Type s = exp(log_scale), q = exp(log_shape);
  for (int k = 0; k < nb; ++k) {
    Type r = Type(d(k)) / s, a;
    if (kernel == 0) {
      a = r * r / Type(2);
      hs[k] = Type(2) * a; hq[k] = Type(0);
    } else if (kernel == 1) {
      a = r; hs[k] = a; hq[k] = Type(0);
    } else {
      // The d == 0 limit is exactly zero, including its shape derivative.
      Type lr = d(k) == 0 ? Type(0) : log(r);
      a = d(k) == 0 ? Type(0) : exp(q * lr);
      hs[k] = q * a; hq[k] = -q * a * lr;
    }
    w[k] = exp(-a);
  }
  out.center = Type(0);
  Type mean_ds = 0, mean_dq = 0;
  // R matrices are column-major. The optional streaming traversal keeps rows
  // contiguous while preserving each row's bin summation order. The reference
  // traversal remains the default, including for existing TMB callers.
  std::vector<Type> num_col, den_col, ns_col, ds_col, nq_col, dq_col;
  if (columnwise) {
    num_col.assign(n, Type(0)); den_col.assign(n, Type(0));
    ns_col.assign(n, Type(0)); ds_col.assign(n, Type(0));
    if (kernel == 2) {
      nq_col.assign(n, Type(0)); dq_col.assign(n, Type(0));
    }
    for (int k = 0; k < nb; ++k) {
      for (int i = 0; i < n; ++i) {
        Type vw = Type(V(i, k)) * w[k], cw = Type(C(i, k)) * w[k];
        num_col[i] += vw; den_col[i] += cw;
        if (derivatives) {
          ns_col[i] += vw * hs[k]; ds_col[i] += cw * hs[k];
          if (kernel == 2) {
            nq_col[i] += vw * hq[k]; dq_col[i] += cw * hq[k];
          }
        }
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    Type num = 0, den = 0, ns = 0, ds = 0, nq = 0, dq = 0;
    if (columnwise) {
      num = num_col[i]; den = den_col[i]; ns = ns_col[i]; ds = ds_col[i];
      if (kernel == 2) { nq = nq_col[i]; dq = dq_col[i]; }
    } else for (int k = 0; k < nb; ++k) {
      Type vw = Type(V(i, k)) * w[k], cw = Type(C(i, k)) * w[k];
      num += vw; den += cw;
      if (derivatives) {
        ns += vw * hs[k]; ds += cw * hs[k];
        nq += vw * hq[k]; dq += cw * hq[k];
      }
    }
    out.raw[i] = out.value[i] = num / den;
    out.center += out.raw[i];
    out.scale[i] = derivatives ? (ns - out.raw[i] * ds) / den : Type(0);
    out.shape[i] = derivatives ? (nq - out.raw[i] * dq) / den : Type(0);
    mean_ds += out.scale[i];
    mean_dq += out.shape[i];
  }
  out.center /= Type(n);
  mean_ds /= Type(n);
  mean_dq /= Type(n);
  out.spread = Type(1);
  if (standardize) {
    Type ss = 0, cross_s = 0, cross_q = 0;
    for (int i = 0; i < n; ++i) {
      Type c = out.raw[i] - out.center;
      ss += c * c;
      cross_s += c * (out.scale[i] - mean_ds);
      cross_q += c * (out.shape[i] - mean_dq);
    }
    out.spread = sqrt(ss / Type(n - 1));
    for (int i = 0; i < n; ++i) {
      out.value[i] = (out.raw[i] - out.center) / out.spread;
      if (derivatives) {
        out.scale[i] = (out.scale[i] - mean_ds) / out.spread -
          out.value[i] * cross_s / ss;
        out.shape[i] = (out.shape[i] - mean_dq) / out.spread -
          out.value[i] * cross_q / ss;
      }
    }
  }
  return out;
}
}
#endif
