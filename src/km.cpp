// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils.h"

// ----------------------------------------------------------------------------
// Kaplan-Meier
// ----------------------------------------------------------------------------

// Tabulate KM: time, nar, surv, haz at each evaluation time.
arma::mat KaplanMeier(
    const arma::colvec& eval_times,
    const arma::colvec& status,
    const arma::colvec& time
) {
  const arma::uword n = time.n_elem;
  arma::colvec unique_times = Union(eval_times, arma::unique(time));
  const arma::uword n_unique_time = unique_times.n_elem;

  arma::colvec censor(n_unique_time);
  arma::colvec death(n_unique_time);
  arma::colvec nar(n_unique_time);
  double current_nar = static_cast<double>(n);

  for (arma::uword i = 0; i < n_unique_time; ++i) {
    const double current_time = unique_times(i);
    const arma::uvec idx = arma::find(time == current_time);
    nar(i) = current_nar;
    censor(i) = arma::accu(status.elem(idx) == 0.0);
    death(i) = arma::accu(status.elem(idx) == 1.0);
    current_nar -= censor(i) + death(i);
  }

  const arma::colvec haz = death / nar;
  const arma::colvec surv = arma::cumprod(1.0 - haz);

  const arma::uword n_eval = eval_times.n_elem;
  arma::colvec nar_out(n_eval);
  arma::colvec haz_out(n_eval);
  arma::colvec surv_out(n_eval);
  arma::uword ptr = 0;

  for (arma::uword i = 0; i < n_unique_time && ptr < n_eval; ++i) {
    if (IsIn(unique_times(i), eval_times)) {
      nar_out(ptr) = nar(i);
      haz_out(ptr) = haz(i);
      surv_out(ptr) = surv(i);
      ++ptr;
    }
  }

  return arma::join_rows(eval_times, nar_out, surv_out, haz_out);
}


// ----------------------------------------------------------------------------
// AUC
// ----------------------------------------------------------------------------

//' Calculate Restricted Mean Survival Time
//'
//' @param status Status, coded as 0 for censoring, 1 for death.
//' @param time Observation time.
//' @param extend Extend AUC calculation if tau exceeds max(time)?
//' @param tau Truncation time.
//' @return Numeric restricted mean survival time.
//' @export
// [[Rcpp::export]]
SEXP RMST(
  const arma::colvec& status,
  const arma::colvec& time,
  bool extend = false,
  Rcpp::Nullable<double> tau = R_NilValue
) {
  arma::colvec unique_times = arma::unique(time);
  double max_time = time.max();
  double trunc_time = max_time;
  if (tau.isNotNull()) {
    trunc_time = Rcpp::as<double>(tau);
    unique_times = Truncate(unique_times, trunc_time, false);
  }

  unique_times = AddLeadVal(unique_times, 0);
  const arma::uword n_times = unique_times.n_elem;

  const arma::mat km_mat = KaplanMeier(unique_times, status, time);
  const arma::colvec surv = km_mat.col(2);

  const arma::colvec delta_t = arma::diff(unique_times);
  const arma::colvec integrand = surv.head(n_times - 1);
  double auc = arma::sum(integrand % delta_t);

  if (extend && trunc_time > max_time) {
    auc += (trunc_time - max_time) * surv.min();
  }
  return Rcpp::wrap(auc);
}


// ----------------------------------------------------------------------------
// Integrate Kaplan Meier
// ----------------------------------------------------------------------------

// Area under the KM curve between a and b (step-function integral).
double IntegrateKM(
    double a,
    double b,
    const arma::colvec& eval_times,
    const arma::colvec& surv
) {
  if (a >= b) return 0.0;

  const arma::uvec idx_low = arma::find(eval_times <= a, 1, "last");
  const arma::uvec idx_high = arma::find(eval_times <= b, 1, "last");
  const double lower_surv = idx_low.n_elem > 0 ? surv(idx_low(0)) : 1.0;
  const double upper_surv = idx_high.n_elem > 0 ? surv(idx_high(0)) : 0.0;

  const arma::uvec mid = arma::find(eval_times > a && eval_times < b);
  arma::colvec int_times = arma::join_cols(
    eval_times.elem(mid),
    arma::vec({a, b})
  );
  int_times = arma::unique(int_times);
  const arma::uword n_t = int_times.n_elem;
  arma::colvec int_vals(n_t);

  for (arma::uword i = 0; i < n_t; ++i) {
    const arma::uvec k = arma::find(eval_times == int_times(i));
    if (i == 0) int_vals(i) = lower_surv;
    else if (i == n_t - 1) int_vals(i) = upper_surv;
    else int_vals(i) = k.n_elem > 0 ? surv(k(0)) : lower_surv;
  }

  const arma::colvec dt = arma::diff(int_times);
  return arma::as_scalar(arma::sum(int_vals.head(dt.n_elem) % dt));
}


// ----------------------------------------------------------------------------
// Martingales.
// ----------------------------------------------------------------------------


// Calculate Martingales Cpp
//
// Construct a subject (row) by evaluation time (col) matrix where
// dM[i, t] = dN[i, t] - Y[i, t]dA[i, t].
// 
// @param eval_times Unique times at which to evaluate the martingale.
// @param haz Value of the hazard at each unique time.
// @param status Subject status.
// @param time Subject observation times.
// @return Matrix with subjects as rows and evaluation times as columns.
arma::mat CalcMartingaleKM(
    const arma::colvec& eval_times,
    const arma::colvec& haz,
    const arma::colvec& status,
    const arma::colvec& time
) {
  const arma::uword n = time.n_elem;
  const arma::uword n_times = eval_times.n_elem;
  arma::mat dm = arma::zeros(n, n_times);

  for (arma::uword i = 0; i < n; ++i) {
    const double subj_time = time(i);
    const double subj_status = status(i);
    for (arma::uword j = 0; j < n_times; ++j) {
      if (eval_times(j) == subj_time && subj_status == 1.0) dm(i, j) += 1.0;
      if (eval_times(j) <= subj_time) dm(i, j) -= haz(j);
      else break;
    }
  }
  return dm;
}


// ----------------------------------------------------------------------------
// Influence functions.
// ----------------------------------------------------------------------------

//' Kaplan-Meier Influence Function
//' 
//' Influence function of the Kaplan-Meier estimator at time t. Specifically,
//' \eqn{\psi_{i}(t) = -S(t)\int_{0}^{t} dM_{i}(u) / Y(u)}.
//' 
//' @param status Status, coded as 0 for censoring, 1 for death.
//' @param time Observation time.
//' @param trunc_time Truncation time.
//' @return Numeric vector of influence function values for each observation.
// [[Rcpp::export]]
SEXP InfluenceKM(
    const arma::colvec& status,
    const arma::colvec& time,
    double trunc_time
) {
  arma::colvec eval_times = Truncate(arma::unique(time), trunc_time);
  const arma::mat km_mat = KaplanMeier(eval_times, status, time);
  const arma::colvec nar = km_mat.col(1);
  const arma::colvec surv = km_mat.col(2);
  const arma::colvec haz = km_mat.col(3);

  // Survival at trunc_time: use last eval time <= trunc_time (avoids empty find).
  const arma::uvec idx_tau = arma::find(eval_times <= trunc_time, 1, "last");
  const double st = (idx_tau.n_elem > 0) ? arma::as_scalar(surv(idx_tau)) : 1.0;

  const arma::mat dm = CalcMartingaleKM(eval_times, haz, status, time);
  const arma::uword n = dm.n_rows;
  arma::colvec influence(n);

  for (arma::uword i = 0; i < n; ++i) {
    influence(i) = -st * static_cast<double>(n) *
      arma::as_scalar(arma::sum(dm.row(i).t() / nar));
  }
  return Rcpp::wrap(influence);
}


//' RMST Influence Function
//' 
//' Influence function of the restricted mean survival time at time t. Specifically,
//' \eqn{\psi_{i}(t) = -S(t)\int_{0}^{t} \mu_{\tau} dM_{i}(u) / Y(u)}.
//' 
//' @param status Status, coded as 0 for censoring, 1 for death.
//' @param time Observation time.
//' @param trunc_time Truncation time.
//' @return Numeric vector of influence function values for each observation.
// [[Rcpp::export]]
SEXP InfluenceRMST(
    const arma::colvec& status,
    const arma::colvec& time,
    double trunc_time
) {
  arma::colvec unique_times = AddLeadVal(arma::unique(time), 0.0);
  unique_times = Truncate(unique_times, trunc_time);
  const arma::uword n_subj = time.n_elem;

  const arma::mat km_mat = KaplanMeier(unique_times, status, time);
  const arma::colvec km_times = km_mat.col(0);
  const arma::colvec nar = km_mat.col(1);
  const arma::colvec surv = km_mat.col(2);
  const arma::colvec haz = km_mat.col(3);
  const arma::colvec par = nar / static_cast<double>(n_subj);

  arma::mat mart = CalcMartingaleKM(km_times, haz, status, time).t();
  const arma::uword n_times = km_times.n_elem;
  arma::colvec mu(n_times);
  for (arma::uword j = 0; j < n_times; ++j) {
    mu(j) = IntegrateKM(km_times(j), trunc_time, km_times, surv);
  }

  arma::colvec psi(n_subj);
  for (arma::uword i = 0; i < n_subj; ++i) {
    psi(i) = -arma::as_scalar(arma::sum(mu / par % mart.col(i)));
  }
  return Rcpp::wrap(psi);
}
