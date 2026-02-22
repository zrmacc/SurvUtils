// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils.h"

// ----------------------------------------------------------------------------
// Cumulative incidence
// ----------------------------------------------------------------------------

struct CICTab {
  arma::colvec time;
  arma::colvec censor;
  arma::colvec event;
  arma::colvec death;
  arma::colvec nar;
  arma::colvec death_rate;
  arma::colvec event_rate;
  arma::colvec haz;
  arma::colvec surv_init;
  arma::colvec cic_event;
  arma::colvec cic_death;
  arma::colvec var_cic_event;
  arma::colvec se_cic_event;
};


CICTab TabulateEventsCpp(
  const arma::colvec& eval_times,
  const arma::colvec& status,
  const arma::colvec& time
) {
  const arma::uword n = time.n_elem;
  arma::colvec unique_times = AddLeadVal(Union(eval_times, arma::unique(time)), 0);
  const arma::uword n_time = unique_times.n_elem;

  arma::colvec n_cens = arma::zeros<arma::colvec>(n_time);
  arma::colvec n_event = arma::zeros<arma::colvec>(n_time);
  arma::colvec n_death = arma::zeros<arma::colvec>(n_time);
  arma::colvec nar = arma::zeros<arma::colvec>(n_time);
  double current_nar = static_cast<double>(n);

  for (arma::uword i = 0; i < n_time; ++i) {
    const double t = unique_times(i);
    const arma::uvec idx = arma::find(time == t);
    nar(i) = current_nar;
    n_cens(i) = arma::accu(status.elem(idx) == 0.0);
    n_event(i) = arma::accu(status.elem(idx) == 1.0);
    n_death(i) = arma::accu(status.elem(idx) == 2.0);
    current_nar -= n_cens(i) + n_event(i) + n_death(i);
  }

  const arma::uword n_eval = eval_times.n_elem;
  arma::colvec censor_out(n_eval);
  arma::colvec event_out(n_eval);
  arma::colvec death_out(n_eval);
  arma::colvec nar_out(n_eval);
  arma::uword j = 0;

  for (arma::uword i = 0; i < n_time && j < n_eval; ++i) {
    if (IsIn(unique_times(i), eval_times)) {
      censor_out(j) = n_cens(i);
      event_out(j) = n_event(i);
      death_out(j) = n_death(i);
      nar_out(j) = nar(i);
      ++j;
    }
  }

  CICTab out;
  out.time = eval_times;
  out.censor = censor_out;
  out.event = event_out;
  out.death = death_out;
  out.nar = nar_out;
  return out;
}


// ----------------------------------------------------------------------------

//' Calculate CIC
//'
//' Estimate the cumulative incidence curve. Specifically:
//' \eqn{F_{1}(t) = P(T \leq t, \delta = 1).}
//' 
//' @param status Status, coded as 0 for censoring, 1 for an event, 2 for death.
//' @param time Observation time.
// [[Rcpp::export]]
SEXP CalcCIC(const arma::colvec& status, const arma::colvec& time) {
  arma::colvec eval_times = AddLeadVal(arma::unique(time), 0);
  CICTab out = TabulateEventsCpp(eval_times, status, time);
  const arma::uword n_time = out.time.n_elem;

  // Avoid division by zero when nar == 0 (e.g. past last event).
  arma::colvec nar_safe = arma::clamp(out.nar, 1.0, arma::datum::inf);
  out.death_rate = out.death / nar_safe;
  out.event_rate = out.event / nar_safe;
  out.haz = (out.death + out.event) / nar_safe;

  out.surv_init.set_size(n_time);
  out.surv_init(0) = 1.0;
  const arma::colvec cumprod_haz = arma::cumprod(1.0 - out.haz);
  if (n_time > 1) {
    out.surv_init.subvec(1, n_time - 1) = cumprod_haz.subvec(0, n_time - 2);
  }

  out.cic_event = arma::cumsum(out.surv_init % out.event_rate);
  out.cic_death = arma::cumsum(out.surv_init % out.death_rate);

  const arma::colvec nar2 = arma::square(nar_safe);
  arma::colvec var1 = arma::square(out.cic_event) % arma::cumsum(out.event / nar2) +
    arma::cumsum(arma::square(1.0 - out.cic_death) % out.event / nar2) -
    2.0 * out.cic_event % arma::cumsum((1.0 - out.cic_death) % out.event / nar2);
  arma::colvec var2 = arma::square(out.cic_event) % arma::cumsum(out.death / nar2) +
    arma::cumsum(arma::square(out.cic_event) % out.death / nar2) -
    2.0 * out.cic_event % arma::cumsum(out.cic_event % out.death / nar2);

  out.var_cic_event = var1 + var2;
  out.se_cic_event = arma::sqrt(out.var_cic_event);

  return Rcpp::DataFrame::create(
    Rcpp::Named("time")=out.time,
    Rcpp::Named("nar")=out.nar,
    Rcpp::Named("censor")=out.censor,
    Rcpp::Named("event")=out.event,
    Rcpp::Named("death")=out.death,
    Rcpp::Named("rate_event")=out.event_rate,
    Rcpp::Named("rate_death")=out.death_rate,
    Rcpp::Named("haz_total")=out.haz,
    Rcpp::Named("surv_init")=out.surv_init,
    Rcpp::Named("cic_event")=out.cic_event,
    Rcpp::Named("cic_death")=out.cic_death,
    Rcpp::Named("var_cic_event")=out.var_cic_event,
    Rcpp::Named("se_cic_event")=out.se_cic_event
  );
};


// ----------------------------------------------------------------------------
// Martingales.
// ----------------------------------------------------------------------------


// Calculate Martingales Cpp
//
// Construct a subject (row) by evaluation time (col) matrix where
// dMj[i, t] = dNj[i, t] - Y[i, t]dAj[i, t].
// 
// @param code Status value for the event of interest.
// @param cshaz Cause-specific hazard at each unique time.
// @param eval_times Unique times at which to evaluate the martingale.
// @param status Subject status.
// @param time Subject observation times.
// @return Matrix with subjects as rows and evaluations times as columns.
arma::mat CalcMartingaleCI(
    int code,
    const arma::colvec& cshaz,
    const arma::colvec& eval_times,
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
      if (eval_times(j) == subj_time && subj_status == static_cast<double>(code)) dm(i, j) += 1.0;
      if (eval_times(j) <= subj_time) dm(i, j) -= cshaz(j);
      else break;
    }
  }
  return dm;
}


// ----------------------------------------------------------------------------

//' CIC Influence Function
//'
//' Influence function for the cumulative incidence of the status == 1 event.
//' 
//' @param status Status, coded as 0 for censoring, 1 for an event, 2 for death.
//' @param time Observation time.
//' @param trunc_time Time at which to evaluate the influence function.
//' @return Vector of per-subject influence function evaluations.
// [[Rcpp::export]]
SEXP InfluenceCIC(
  const arma::colvec& status,
  const arma::colvec& time,
  double trunc_time
) {
  arma::colvec eval_times = Truncate(arma::unique(time), trunc_time);
  const arma::uword n = time.n_elem;
  CICTab tab = TabulateEventsCpp(eval_times, status, time);
  const arma::uword n_time = tab.time.n_elem;

  arma::colvec nar_safe = arma::clamp(tab.nar, 1.0, arma::datum::inf);
  tab.death_rate = tab.death / nar_safe;
  tab.event_rate = tab.event / nar_safe;
  tab.haz = (tab.death + tab.event) / nar_safe;

  tab.surv_init.set_size(n_time);
  tab.surv_init(0) = 1.0;
  if (n_time > 1) {
    const arma::colvec cumprod_haz = arma::cumprod(1.0 - tab.haz);
    tab.surv_init.subvec(1, n_time - 1) = cumprod_haz.subvec(0, n_time - 2);
  }

  tab.cic_event = arma::cumsum(tab.surv_init % tab.event_rate);
  tab.cic_death = arma::cumsum(tab.surv_init % tab.death_rate);

  const arma::mat dM1 = CalcMartingaleCI(1, tab.event_rate, tab.time, status, time);
  const arma::mat dM2 = CalcMartingaleCI(2, tab.death_rate, tab.time, status, time);
  const arma::mat dM = dM1 + dM2;

  // CIC at trunc_time: use last eval time <= trunc_time (avoids empty find).
  const arma::uvec idx_tau = arma::find(eval_times <= trunc_time, 1, "last");
  const double ft = (idx_tau.n_elem > 0) ? arma::as_scalar(tab.cic_event(idx_tau)) : 0.0;

  arma::colvec influence(n);
  for (arma::uword i = 0; i < n; ++i) {
    const arma::colvec dMi = dM.row(i).t();
    const arma::colvec dM1i = dM1.row(i).t();
    const double nn = static_cast<double>(n);
    influence(i) = -ft * nn * arma::as_scalar(arma::sum(dMi / nar_safe)) +
      nn * arma::as_scalar(arma::sum(tab.cic_event % dMi / nar_safe)) +
      nn * arma::as_scalar(arma::sum(tab.surv_init % dM1i / nar_safe));
  }
  return Rcpp::wrap(influence);
}

