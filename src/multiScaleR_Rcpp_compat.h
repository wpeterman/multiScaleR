#ifndef MULTISCALER_RCPP_COMPAT_H
#define MULTISCALER_RCPP_COMPAT_H

#include <Rversion.h>

/*
 * CRAN Windows R-devel snapshots in April 2026 omitted the
 * R_NamespaceRegistry declaration from installed headers while the Rcpp
 * headers still referenced it. Declare it before Rcpp is included so package
 * compilation can proceed on that transient toolchain. Do not include
 * Rinternals.h here: including it before Rcpp can expose R's compatibility
 * macros (e.g., length) to C++ standard library headers.
 */
#if defined(_WIN32) && R_VERSION >= R_Version(4, 6, 0)
struct SEXPREC;
typedef struct SEXPREC *SEXP;
extern SEXP R_NamespaceRegistry;
#endif

#endif
