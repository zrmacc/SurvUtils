// [[Rcpp::depends(RcppArmadillo)]]
#include "utils.h"

namespace {

// Check if value is in vector (vectorized, avoids loop).
inline bool IsInImpl(double a, const arma::colvec& b) {
  return arma::any(b == a);
}

}  // namespace

bool IsIn(double a, const arma::colvec& b) {
  return IsInImpl(a, b);
}

// Union of two time vectors; preserves order and uniqueness.
// Second parameter by value so symbol matches ABI expected by existing object files.
arma::colvec Union(const arma::colvec& a, arma::colvec b) {
  const double a0 = a(0);
  for (arma::uword i = 0; i < b.n_elem; ++i) {
    if (IsInImpl(b(i), a)) b(i) = a0;
  }
  return arma::unique(arma::join_cols(a, b));
}

// Truncate times at tau; optionally append tau if not present.
arma::colvec Truncate(const arma::colvec& time, double tau, bool add_tau) {
  arma::colvec unique_times = arma::unique(time);
  for (arma::uword i = 0; i < unique_times.n_elem; ++i) {
    if (unique_times(i) > tau) unique_times(i) = tau;
  }
  if (add_tau && !arma::any(unique_times == tau)) {
    const arma::uword n = unique_times.n_elem;
    arma::colvec out(n + 1);
    out.head(n) = unique_times;
    out(n) = tau;
    return arma::unique(out);
  }
  return arma::unique(unique_times);
}

arma::colvec AddLeadVal(const arma::colvec& x, const double value) {
  if (arma::any(x == value)) return x;
  arma::colvec out(x.n_elem + 1);
  out(0) = value;
  out.tail(x.n_elem) = x;
  return out;
}

