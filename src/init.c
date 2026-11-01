#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare function prototypes (MUST match RcppExports.cpp)
extern SEXP _multiScaleR_msr_surface_sq_derivative_cpp(SEXP, SEXP, SEXP);
extern SEXP _multiScaleR_ci_func_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _multiScaleR_landscape_composition_metric_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multiScaleR_landscape_edge_metric_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multiScaleR_landscape_adjacency_metric_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multiScaleR_scale_type_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multiScaleR_scale_type_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);  // Add this line
extern SEXP _multiScaleR_surface_sa_focal_cpp(SEXP, SEXP);
extern SEXP _multiScaleR_msr_kernel_eval_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multiScaleR_msr_kernel_stream_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multiScaleR_msr_kernel_validate_cpp(SEXP, SEXP, SEXP);

extern SEXP _multiScaleR_msr_adjacency_prepare_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _multiScaleR_msr_adjacency_eval_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multiScaleR_msr_adjacency_events_prepare_cpp(SEXP);
extern SEXP _multiScaleR_msr_adjacency_events_eval_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP _multiScaleR_msr_surface_weighted_sq_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multiScaleR_msr_surface_sq_prefix_cpp(SEXP);

extern SEXP _multiScaleR_msr_choice_softmax_cpp(SEXP, SEXP);

// Register native routines
static const R_CallMethodDef CallEntries[] = {
  {"_multiScaleR_msr_choice_softmax_cpp", (DL_FUNC) &_multiScaleR_msr_choice_softmax_cpp, 2},
  {"_multiScaleR_ci_func_cpp", (DL_FUNC) &_multiScaleR_ci_func_cpp, 4},
  {"_multiScaleR_landscape_composition_metric_cpp", (DL_FUNC) &_multiScaleR_landscape_composition_metric_cpp, 9},
  {"_multiScaleR_landscape_edge_metric_cpp", (DL_FUNC) &_multiScaleR_landscape_edge_metric_cpp, 7},
  {"_multiScaleR_landscape_adjacency_metric_cpp", (DL_FUNC) &_multiScaleR_landscape_adjacency_metric_cpp, 8},
  {"_multiScaleR_scale_type_cpp", (DL_FUNC) &_multiScaleR_scale_type_cpp, 6},
  {"_multiScaleR_scale_type_sparse", (DL_FUNC) &_multiScaleR_scale_type_sparse, 6},  // Add this entry
  {"_multiScaleR_surface_sa_focal_cpp", (DL_FUNC) &_multiScaleR_surface_sa_focal_cpp, 2},
  {"_multiScaleR_msr_kernel_eval_cpp", (DL_FUNC) &_multiScaleR_msr_kernel_eval_cpp, 7},
  {"_multiScaleR_msr_kernel_stream_cpp", (DL_FUNC) &_multiScaleR_msr_kernel_stream_cpp, 7},
  {"_multiScaleR_msr_kernel_validate_cpp", (DL_FUNC) &_multiScaleR_msr_kernel_validate_cpp, 3},
  {"_multiScaleR_msr_adjacency_prepare_cpp", (DL_FUNC) &_multiScaleR_msr_adjacency_prepare_cpp, 4},
  {"_multiScaleR_msr_adjacency_eval_cpp", (DL_FUNC) &_multiScaleR_msr_adjacency_eval_cpp, 5},
  {"_multiScaleR_msr_adjacency_events_prepare_cpp", (DL_FUNC) &_multiScaleR_msr_adjacency_events_prepare_cpp, 1},
  {"_multiScaleR_msr_adjacency_events_eval_cpp", (DL_FUNC) &_multiScaleR_msr_adjacency_events_eval_cpp, 5},
  {"_multiScaleR_msr_surface_weighted_sq_cpp", (DL_FUNC) &_multiScaleR_msr_surface_weighted_sq_cpp, 5},
  {"_multiScaleR_msr_surface_sq_prefix_cpp", (DL_FUNC) &_multiScaleR_msr_surface_sq_prefix_cpp, 1},
  {"_multiScaleR_msr_surface_sq_derivative_cpp", (DL_FUNC) &_multiScaleR_msr_surface_sq_derivative_cpp, 3},
  {NULL, NULL, 0}  // Sentinel to mark the end
};

// Package initialization
void R_init_multiScaleR(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
