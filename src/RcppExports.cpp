// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CalcCIC
SEXP CalcCIC(const arma::vec& status, const arma::vec& time);
RcppExport SEXP _SurvUtils_CalcCIC(SEXP statusSEXP, SEXP timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type time(timeSEXP);
    rcpp_result_gen = Rcpp::wrap(CalcCIC(status, time));
    return rcpp_result_gen;
END_RCPP
}
// InfluenceCIC
SEXP InfluenceCIC(const arma::colvec& status, const arma::colvec& time, const double trunc_time);
RcppExport SEXP _SurvUtils_InfluenceCIC(SEXP statusSEXP, SEXP timeSEXP, SEXP trunc_timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const double >::type trunc_time(trunc_timeSEXP);
    rcpp_result_gen = Rcpp::wrap(InfluenceCIC(status, time, trunc_time));
    return rcpp_result_gen;
END_RCPP
}
// RMST
SEXP RMST(const arma::colvec status, const arma::colvec time, const bool extend, Rcpp::Nullable<double> tau);
RcppExport SEXP _SurvUtils_RMST(SEXP statusSEXP, SEXP timeSEXP, SEXP extendSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const bool >::type extend(extendSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(RMST(status, time, extend, tau));
    return rcpp_result_gen;
END_RCPP
}
// InfluenceKM
SEXP InfluenceKM(const arma::colvec status, const arma::colvec time, const double trunc_time);
RcppExport SEXP _SurvUtils_InfluenceKM(SEXP statusSEXP, SEXP timeSEXP, SEXP trunc_timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const double >::type trunc_time(trunc_timeSEXP);
    rcpp_result_gen = Rcpp::wrap(InfluenceKM(status, time, trunc_time));
    return rcpp_result_gen;
END_RCPP
}
// InfluenceRMST
SEXP InfluenceRMST(const arma::colvec status, const arma::colvec time, const double trunc_time);
RcppExport SEXP _SurvUtils_InfluenceRMST(SEXP statusSEXP, SEXP timeSEXP, SEXP trunc_timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const double >::type trunc_time(trunc_timeSEXP);
    rcpp_result_gen = Rcpp::wrap(InfluenceRMST(status, time, trunc_time));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SurvUtils_CalcCIC", (DL_FUNC) &_SurvUtils_CalcCIC, 2},
    {"_SurvUtils_InfluenceCIC", (DL_FUNC) &_SurvUtils_InfluenceCIC, 3},
    {"_SurvUtils_RMST", (DL_FUNC) &_SurvUtils_RMST, 4},
    {"_SurvUtils_InfluenceKM", (DL_FUNC) &_SurvUtils_InfluenceKM, 3},
    {"_SurvUtils_InfluenceRMST", (DL_FUNC) &_SurvUtils_InfluenceRMST, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_SurvUtils(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
