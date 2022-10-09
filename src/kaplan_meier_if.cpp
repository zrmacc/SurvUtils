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

bool IsIn(const double &a, const arma::vec &b) {
  
  for(int i=0; i<b.size(); i++) {
    if(b(i) == a) {
      return true;
    }
  }
  return false;
}


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
}


// Truncate
//
// @param time Vector of time points.
// @param tau Truncation time.
// @return Truncated vector.
arma::colvec Truncate(const arma::colvec &time, const double tau) {
  arma::colvec unique_times = arma::unique(time);
  for(int i=0; i<unique_times.size(); i++){
    if(unique_times(i) > tau){
      unique_times(i) = tau;
    }
  }
  arma::colvec out = arma::unique(unique_times);
  return out;
}


// ----------------------------------------------------------------------------
// At risk matrix.
// ----------------------------------------------------------------------------

// Tabulate At-Risk Matrix Cpp
//
// Construct a subject (row) by evaluation time (col) matrix.
//  
// @param eval_times Evaluation times.
// @param idx Unique subject index. 
// @param time Observation time.
// @return Numeric matrix.

arma::mat AtRiskMatrixCpp(
    const arma::colvec eval_times,
    const arma::colvec idx,
    const arma::colvec time
){
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
  // Create a subject by evaluation times matrix, where
  // Y[i, t] is 1 if subject i is at-risk at time t.
  const int n_times = eval_times.size();
  arma::mat y = arma::zeros(n, n_times);
  
  // Loop over subjects.
  for(int i=0; i<n; i++) {
    
    // Current subject.
    arma::colvec subj_times = time.elem(arma::find(idx == unique_idx(i)));
    
    // Subject is at risk until time > the subject's last time.
    double subj_last_time = subj_times.max();
    
    // Loop over times.
    for(int j=0; j<n_times; j++) {
      
      double current_time = eval_times(j);

      if(current_time <= subj_last_time) {
        y(i, j) = 1;
      } else {
        break;
      }
    }
  }
  return y;
}


// ----------------------------------------------------------------------------
// Kaplan Meier.
// ----------------------------------------------------------------------------

// Tabulate Kaplan Meier Cpp
//  
// Constructs a matrix with evaluation times as rows, and 4 columns:
// \itemize{
// \item{time}{Evaluation times.}
// \item{nar}{Number at risk.}
// \item{surv}{Survival probability.}
// \item{haz}{Hazard.}
// }
//
// @param eval_times Evaluation times.
// @param idx Unique subject index.
// @param status Status, coded as 0 for censoring, 1 for death.
// @param time Observation time.
// @return Numeric matrix.

arma::mat KaplanMeierCpp(
    const arma::colvec eval_times,
    const arma::colvec idx,
    const arma::colvec status,
    const arma::colvec time
){
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
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
// Martingales.
// ----------------------------------------------------------------------------


// Calculate Martingales Cpp
//
// Construct a subject (row) by evaluation time (col) matrix where
// dM[i, t] = dN[i, t] - Y[i, t]dA[i, t].
// 
// @param haz Value of the hazard at each unique time.
// @param idx Subject index.
// @param status Subject status.
// @param time Subject observation times.
// @param unique_times Unique times at which to obtain the martingale.
// @return Matrix with subjects as rows and unique times as columns.

arma::mat CalcMartingaleCpp(
    const arma::colvec haz,
    const arma::colvec idx,
    const arma::colvec status,
    const arma::colvec time,
    const arma::colvec unique_times
) {
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
  // Unique times.
  const int n_times = unique_times.size();
  
  // Create a subject by evaluation times matrix, where
  // dM[i, t] is the martingale increment for subject i at time t.
  arma::mat dm = arma::zeros(n, n_times);
  
  // Loop over subjects.
  for(int i=0; i<n; i++) {
    
    // Time and status for the focus subject.
    arma::colvec subj_times = time.elem(arma::find(idx == unique_idx(i)));
    arma::colvec subj_status = status.elem(arma::find(idx == unique_idx(i)));
    
    // Subject is at risk until time > the subject's las time.
    double subj_last_time = subj_times.max();
    
    // Subject's final status.
    arma::uvec key_final_status = arma::find(
      subj_times == subj_last_time, 1, "last");
    double final_status = arma::as_scalar(
      subj_status.elem(key_final_status)
    );
    
    // Loop over times.
    for(int j=0; j<n_times; j++) {
      
      double current_time = unique_times(j);
      double current_haz = haz(j);
      
      // Add dN_{i}(t).
      if(current_time == subj_last_time) {
        if(final_status != 0) {
          dm(i, j) += 1;
        }
      }
      
      // Add -Y_{i}(t)d\Lambda(t).
      if(current_time <= subj_last_time) {
        dm(i, j) -= current_haz;
      } else {
        break;  
      }
    }
  }
  return dm;
}


// ----------------------------------------------------------------------------
// Influence functions.
// ----------------------------------------------------------------------------

//' Influence Function R
//' 
//' Influence function of the Kaplan-Meier estimator at time t. Specifically,
//' \eqn{\psi_{i}(t) = -S(t)\int_{0}^{t} dM_{i}(u) / Y(u)}.
//'  
//' @param idx Unique subject index. 
//' @param status Status, coded as 0 for censoring, 1 for death.
//' @param time Observation time.
//' @param trunc_time Truncation time.
//' @return Data.frame.
// [[Rcpp::export]]

SEXP InfluenceKM(
    const arma::colvec idx,
    const arma::colvec status,
    const arma::colvec time,
    const float trunc_time 
){

	// Evaluation times.
	arma::colvec eval_times = Truncate(arma::unique(time), trunc_time);
	arma::mat km_mat = KaplanMeierCpp(eval_times, idx, status, time);
	// Rcpp::Rcout << km_mat << std::endl; 

	const arma::colvec y = km_mat.col(1); // Number at risk.
	const arma::colvec surv = km_mat.col(2); // Survival.
	const arma::colvec haz = km_mat.col(3); // Instantaneous hazard.

	// Survival at truncation time.
	double st = arma::as_scalar(surv.elem(arma::find(eval_times == trunc_time)));

	// Martingales.
	const arma::mat dm = CalcMartingaleCpp(haz, idx, status, time, eval_times);
	// Rcpp::Rcout << dm << std::endl; 

	// Influence functions.
	int n = dm.n_rows;
	arma::colvec dmi;
	arma::colvec influence = arma::zeros(n);

	for(int i=0; i<n; i++) {
		dmi = arma::trans(dm.row(i));
		influence(i) = -1 * st * n * arma::sum(dmi / y);	
	}

	return Rcpp::wrap(influence);
}
