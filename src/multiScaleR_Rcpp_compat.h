#ifndef MULTISCALER_RCPP_COMPAT_H
#define MULTISCALER_RCPP_COMPAT_H

#include <Rversion.h>

/*
 * CRAN Windows R-devel snapshots in April 2026 omitted the
 * R_NamespaceRegistry declaration from installed headers while the Rcpp
 * headers still referenced it. Declare it before Rcpp is included so package
 * compilation can proceed on that transient toolchain.
 */
#if defined(_WIN32) && R_VERSION >= R_Version(4, 6, 0)
#include <Rinternals.h>
extern SEXP R_NamespaceRegistry;
#endif

#endif
