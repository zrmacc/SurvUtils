// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils.h"

// For debugging: Rcpp::Rcout << << std::endl; 

// ----------------------------------------------------------------------------
// Kaplan Meier.
// ----------------------------------------------------------------------------

// Tabulate Kaplan Meier Cpp
//  
// Constructs a matrix with evaluation times as rows, and 4 columns:
// \itemize{
// \item{time: Evaluation times.}
// \item{nar: Number at risk.}
// \item{surv: Survival probability.}
// \item{haz: Hazard.}
// }
//
// @param eval_times Evaluation times.
// @param status Status, coded as 0 for censoring, 1 for death.
// @param time Observation time.
// @return Numeric matrix.
arma::mat KaplanMeier(
    const arma::colvec eval_times,
    const arma::colvec status,
    const arma::colvec time
){
  
  // Subjects.
  const int n = time.size();
  
  // Unique times.
  arma::colvec unique_times = Union(eval_times, arma::unique(time));
  const int n_unique_time = unique_times.size();
  
  // Censoring, death, and at risk counts.
  arma::colvec censor(n_unique_time);
  arma::colvec death(n_unique_time);
  arma::colvec nar(n_unique_time);
  
  // Loop over unique times.
  double current_nar = n;

  for(int i=0; i<n_unique_time; i++) {
    
    double current_time = unique_times(i);
    const arma::colvec current_status = status.elem(arma::find(time == current_time));
    
    nar(i) = current_nar;
    censor(i) = arma::sum(current_status == 0.0);
    death(i) = arma::sum(current_status == 1.0);
    
    // Update NAR.
    current_nar -= censor(i) + death(i);
  }
  
  // Hazard.
  const arma::colvec haz = death / nar;
  
  // Survival probability.
  arma::colvec surv = arma::cumprod(1 - haz);
  
  // Restrict to evaluation times.
  const int n_eval_time = eval_times.size();
  arma::colvec nar_out(n_eval_time);
  arma::colvec haz_out(n_eval_time);
  arma::colvec surv_out(n_eval_time);

  int pointer = 0;
  for(int i=0; i<n_unique_time; i++) {
    
    double current_time = unique_times(i);
    if(IsIn(current_time, eval_times)) {
      nar_out(pointer) = nar(i);
      haz_out(pointer) = haz(i);
      surv_out(pointer) = surv(i);
      pointer += 1;
    }
    
  }
  
  // Output.
  arma::mat out = arma::join_rows(eval_times, nar_out, surv_out, haz_out);
  return out;
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
  const arma::colvec status,
  const arma::colvec time,
  const bool extend=false,
  Rcpp::Nullable<double> tau=R_NilValue
){
  
  // Unique times.
  arma::colvec unique_times = arma::unique(time);
  double max_time = time.max();
  double trunc_time = max_time;
  if (tau.isNotNull()) {
    trunc_time = Rcpp::as<double>(tau);
    unique_times = Truncate(unique_times, trunc_time, false);
  }

  // Add leading zero.
  unique_times = AddLeadVal(unique_times, 0);
  const int n_times = unique_times.size();

  // Calculate Kaplan-Meier curve.
  const arma::mat km_mat = KaplanMeier(unique_times, status, time);
  arma::colvec km_times = km_mat.col(0);
  arma::colvec surv = km_mat.col(2);
  
  // Rcpp::Rcout << km_times << std::endl; 
  // Rcpp::Rcout << surv << std::endl; 

  // Calculate AUC.
  const arma::colvec delta_t = arma::diff(unique_times);
  const arma::colvec integrand = surv.subvec(0, n_times - 2);
  double auc = arma::sum(integrand % delta_t);

  // Project past end of Kaplan-Meier curve, if applicable. 
  if (extend && (trunc_time > max_time)) {
    auc += (trunc_time - max_time) * surv.min(); 
  }
  return Rcpp::wrap(auc);
}


// ----------------------------------------------------------------------------
// Integrate Kaplan Meier
// ----------------------------------------------------------------------------

// Integrate Kaplan Meier
//
// Integrate the Kaplan-Meier curve between two values.
// 
// @param a limit of integration.
// @param b limit of integration.
// @param eval_times Times at which the Kaplan-Meier curve is evaluated.
// @param surv Estimated survival at the evaluation times.
// @return Area under the KM curve between the type values.
double IntegrateKM(
    const double a, 
    const double b,
    const arma::colvec eval_times,    
    const arma::colvec surv
) {

  if (a == b) {
    return 0;
  }
  
  // Determine survival at a.
  arma::colvec lower_time_vec = eval_times.elem(arma::find(eval_times <= a, 1, "last"));
  double lower_time = lower_time_vec.at(0);
  arma::colvec lower_surv_vec = surv.elem(arma::find(eval_times == lower_time));
  double lower_surv = lower_surv_vec.at(0);
  
  // Determine survival at b.
  arma::colvec upper_time_vec = eval_times.elem(arma::find(eval_times <= b, 1, "last"));
  const double upper_time = upper_time_vec.at(0);
  arma::colvec upper_surv_vec = surv.elem(arma::find(eval_times == upper_time));
  const double upper_surv = upper_surv_vec.at(0);
  
  // Subset to the evaluation times strictly between a and b.
  // Append a and b to the evaluation times.
  // Append surv(a) and surv(b) to the survival values.
  const arma::colvec eval_times_subset = eval_times.elem(arma::find(
    (eval_times > a) && (eval_times < b)
  ));
  const arma::colvec surv_subset = surv.elem(arma::find(
    (eval_times > a) && (eval_times < b)
  ));
  arma::colvec int_times = arma::join_cols(eval_times_subset, arma::vec({a, b}));

  // Sort and reduce to unique values.
  int_times = arma::unique(int_times);
  const double n_times = int_times.n_elem;
  
  // Construct vector of survival values corresponding to integration times.
  arma::colvec int_values = arma::zeros(n_times);
  for (int i=0; i<n_times; i++){
    if (i == 0) {
      int_values(i) = lower_surv;
    } else if (i == n_times - 1) {
      int_values(i) = upper_surv;
    } else {
      arma::colvec lookup_surv = surv.elem(arma::find(eval_times == int_times(i)));
      int_values(i) = lookup_surv.at(0);
    }
  }
  
  // Integrate.
  const arma::colvec delta_t = arma::diff(int_times);
  const arma::colvec integrand = int_values.subvec(0, delta_t.n_elem - 1);
  const double out = arma::sum(integrand % delta_t);
  return out;
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
// @return Matrix with subjects as rows and unique times as columns.
arma::mat CalcMartingaleKM(
    const arma::colvec &eval_times,    
    const arma::colvec &haz,
    const arma::colvec &status,
    const arma::colvec &time
) {
  
  // Subjects.
  const int n = time.size();
  
  // Unique times.
  const int n_times = eval_times.size();
  
  // Create a subject by evaluation times matrix, where
  // dM[i, t] is the martingale increment for subject i at time t.
  arma::mat dm = arma::zeros(n, n_times);
  
  // Loop over subjects.
  for(int i=0; i<n; i++) {
    
    // Time and status for the focus subject.
    const double subj_time = time(i);
    const int subj_status = status(i);

    // Loop over times.
    for(int j=0; j<n_times; j++) {
      
      const double current_time = eval_times(j);
      const double current_haz = haz(j);
      
      // Add dN_{i}(t).
      if(current_time == subj_time && subj_status == 1) {
        dm(i, j) += 1;
      }
      
      // Add -Y_{i}(t)dA(t).
      if(current_time <= subj_time) {
        dm(i, j) -= current_haz;
      } else {
        break;  
      }

    } // End loop over times.
  } // End loop over subjects.
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
    const arma::colvec status,
    const arma::colvec time,
    const double trunc_time 
){

	// Evaluation times.
	arma::colvec eval_times = Truncate(arma::unique(time), trunc_time);
	arma::mat km_mat = KaplanMeier(eval_times, status, time);
	// Rcpp::Rcout << km_mat << std::endl; 

	const arma::colvec nar = km_mat.col(1); // Number at risk.
	const arma::colvec surv = km_mat.col(2); // Survival.
	const arma::colvec haz = km_mat.col(3); // Instantaneous hazard.

	// Survival at truncation time.
	double st = arma::as_scalar(surv.elem(arma::find(eval_times == trunc_time)));

	// Martingales.
	const arma::mat dm = CalcMartingaleKM(eval_times, haz, status, time);
	// Rcpp::Rcout << dm << std::endl; 

	// Influence functions.
	int n = dm.n_rows;
	arma::colvec dmi;
	arma::colvec influence = arma::zeros(n);

	for(int i=0; i<n; i++) {
		dmi = arma::trans(dm.row(i));
		influence(i) = -1 * st * n * arma::sum(dmi / nar);	
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
    const arma::colvec status,
    const arma::colvec time,
    const double trunc_time 
){

  // Unique times.
  arma::colvec unique_times = arma::unique(time);
  
  // Add leading zero.
  unique_times = AddLeadVal(unique_times, 0);
  unique_times = Truncate(unique_times, trunc_time);
  const int n_subj = time.n_elem;
  
  // Calculate Kaplan-Meier curve.
  const arma::mat km_mat = KaplanMeier(unique_times, status, time);
  const arma::colvec km_times = km_mat.col(0);
  const arma::colvec nar = km_mat.col(1);
  const arma::colvec surv = km_mat.col(2);
  const arma::colvec haz = km_mat.col(3);
  
  // Proportion at risk.
  const arma::colvec par = nar / n_subj;
  
  // Calculate martingales.
  arma::mat mart = CalcMartingaleKM(km_times, haz, status, time);
  mart = mart.t();
  
  // Calculate mu;
  const double n_times = km_times.n_elem;
  arma::colvec mu = arma::zeros(n_times);
  for (int j=0; j<n_times; j++) {
    mu(j) = IntegrateKM(km_times(j), trunc_time, km_times, surv);
  }
  
  // Calculate influence function.
  arma::colvec psi = arma::zeros(n_subj);
  for (int i=0; i<n_subj; i++) {
    psi(i) = -1 * arma::sum(mu / par % mart.col(i));
  }
  
  return Rcpp::wrap(psi);
}
