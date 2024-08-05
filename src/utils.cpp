// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// For debugging: Rcpp::Rcout << << std::endl; 

// ----------------------------------------------------------------------------
// Utilities.
// ----------------------------------------------------------------------------


// Check if Value is in Vector
// 
// @param a Value to search for.
// @param b Vector to search.
// @return bool.
bool IsIn(const double a, const arma::colvec &b) {
  
  for(int i=0; i<b.size(); i++) {
    if(b(i) == a) {
      return true;
    }
  }
  return false;
};


// Union
// 
// @param a Vector.
// @param b Vector.
// @return Vector.
arma::colvec Union(const arma::colvec &a, arma::colvec b) {
  
  // Sets all elements of b that are in a to a value that is known
  // to be in a. Then concatenates a and b and filters to unique values.
  
  double a0 = a(0);
  for(int i=0; i<b.size(); i++) {
    if(IsIn(b(i), a)){
      b(i) = a0;
    }
  }

  arma::colvec out = arma::unique(arma::join_cols(a, b));
  return out;
};


// Truncate
//
// @param time Vector of time points.
// @param tau Truncation time.
// @param add_tab Add truncation time if not present?
// @return Truncated vector.
arma::colvec Truncate(
  const arma::colvec &time,
  const double tau,
  const bool add_tau=false
) {
  arma::colvec unique_times = arma::unique(time);
  for (int i=0; i<unique_times.size(); i++) {
    if (unique_times(i) > tau) {
      unique_times(i) = tau;
    }
  }

  // Add truncation time, if not already present.
  arma::colvec out;
  if (add_tau && !arma::any(unique_times == tau)) {
    int n = unique_times.n_elem;
    out = arma::zeros(n + 1);
    out.subvec(0, n - 1) = unique_times;
    out(n) = tau;
  } else {
    out = unique_times;
  }

  return arma::unique(out);
};


// Add Leading Value
//
// @param x Input vector.
// @param value Value to insert at leading position.
arma::colvec AddLeadVal(const arma::colvec &x, const double value) {
  const int n = x.n_elem;
  arma::colvec out;
  if (!arma::any(x == value)) {
    out = arma::zeros(n + 1);
    out(0) = value;
    out.subvec(1, n) = x;
  } else {
    out = x;
  }
  return out;
};

