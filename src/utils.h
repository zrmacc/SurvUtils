// utils.h
#ifndef SURVUTILS_UTILS_H
#define SURVUTILS_UTILS_H

#include <RcppArmadillo.h>

arma::colvec AddLeadVal(const arma::colvec& x, const double value);
bool IsIn(double a, const arma::colvec& b);
arma::colvec Truncate(const arma::colvec& time, double tau, bool add_tau = false);
arma::colvec Union(const arma::colvec& a, arma::colvec b);

#endif
