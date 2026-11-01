#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>

// Registered directly in init.c. No RNG scope: normalization is deterministic.
RcppExport SEXP _multiScaleR_msr_choice_softmax_cpp(SEXP etaSEXP, SEXP groupsSEXP) {
BEGIN_RCPP
  if (TYPEOF(etaSEXP) != REALSXP && TYPEOF(etaSEXP) != INTSXP)
    Rcpp::stop("Choice scores must be numeric.");
  Rcpp::NumericVector eta = Rcpp::as<Rcpp::NumericVector>(etaSEXP);
  Rcpp::List groups(groupsSEXP);
  const R_xlen_t n = eta.size();
  Rcpp::NumericVector probability(n), log_probability(n);
  std::vector<unsigned char> seen(static_cast<size_t>(n), 0);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (!R_FINITE(eta[i])) Rcpp::stop("Choice scores must be finite.");
  }
  for (R_xlen_t s = 0; s < groups.size(); ++s) {
    if (s % 1024 == 0) Rcpp::checkUserInterrupt();
    SEXP index = groups[s];
    const bool integer = TYPEOF(index) == INTSXP;
    if (!integer && TYPEOF(index) != REALSXP)
      Rcpp::stop("Choice-set row indices must be numeric integers.");
    const R_xlen_t count = Rf_xlength(index);
    if (!count) Rcpp::stop("Choice sets must not be empty.");
    // Read the retained index vector directly, without per-set coercion copies.
    const auto at = [&](R_xlen_t j) -> double {
      return integer ? static_cast<double>(INTEGER(index)[j]) : REAL(index)[j];
    };
    double shift = -std::numeric_limits<double>::infinity();
    for (R_xlen_t j = 0; j < count; ++j) {
      const double row = at(j);
      if (!R_FINITE(row) || row < 1 || row > static_cast<double>(n) ||
          std::floor(row) != row)
        Rcpp::stop("Choice-set row indices are invalid.");
      const R_xlen_t i = static_cast<R_xlen_t>(row) - 1;
      if (seen[static_cast<size_t>(i)])
        Rcpp::stop("Choice sets must partition rows without overlap.");
      seen[static_cast<size_t>(i)] = 1;
      if (eta[i] > shift) shift = eta[i];
    }
    long double total = 0;
    for (R_xlen_t j = 0; j < count; ++j) {
      const R_xlen_t i = static_cast<R_xlen_t>(at(j)) - 1;
      probability[i] = std::exp(eta[i] - shift);
      total += probability[i];
    }
    const double denominator = static_cast<double>(total);
    const double log_denominator = std::log(denominator);
    for (R_xlen_t j = 0; j < count; ++j) {
      const R_xlen_t i = static_cast<R_xlen_t>(at(j)) - 1;
      probability[i] /= denominator;
      log_probability[i] = eta[i] - shift - log_denominator;
    }
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    if (!seen[static_cast<size_t>(i)])
      Rcpp::stop("Choice sets must partition every prediction row.");
  }
  return Rcpp::List::create(Rcpp::_["probability"] = probability,
                            Rcpp::_["log_probability"] = log_probability);
END_RCPP
}
