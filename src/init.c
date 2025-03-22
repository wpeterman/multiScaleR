#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare function prototypes (MUST match RcppExports.cpp)
extern SEXP _multiScaleR_ci_func_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _multiScaleR_scale_type_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multiScaleR_scale_type_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);  // Add this line

// Register native routines
static const R_CallMethodDef CallEntries[] = {
  {"_multiScaleR_ci_func_cpp", (DL_FUNC) &_multiScaleR_ci_func_cpp, 4},
  {"_multiScaleR_scale_type_cpp", (DL_FUNC) &_multiScaleR_scale_type_cpp, 6},
  {"_multiScaleR_scale_type_sparse", (DL_FUNC) &_multiScaleR_scale_type_sparse, 6},  // Add this entry
  {NULL, NULL, 0}  // Sentinel to mark the end
};

// Package initialization
void R_init_multiScaleR(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
