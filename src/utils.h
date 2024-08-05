// utils.h
#ifndef UTILS_H
#define UTILS_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Function declarations.
arma::colvec AddLeadVal(const arma::colvec &x, const double value);
bool IsIn(const double a, const arma::colvec &b);
arma::colvec Truncate(const arma::colvec &time, const double tau, const bool add_tau=false);
arma::colvec Union(const arma::colvec &a, arma::colvec b);

#endif // UTILS_H
