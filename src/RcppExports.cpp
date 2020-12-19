// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// pair_dist2
double pair_dist2(arma::vec x, arma::vec y);
RcppExport SEXP _nett_pair_dist2(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(pair_dist2(x, y));
    return rcpp_result_gen;
END_RCPP
}
// pair_dist2_mat
arma::mat pair_dist2_mat(arma::mat z);
RcppExport SEXP _nett_pair_dist2_mat(SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(pair_dist2_mat(z));
    return rcpp_result_gen;
END_RCPP
}
// sample_dclvm_cpp
arma::sp_mat sample_dclvm_cpp(arma::mat z, double scale, arma::vec theta);
RcppExport SEXP _nett_sample_dclvm_cpp(SEXP zSEXP, SEXP scaleSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_dclvm_cpp(z, scale, theta));
    return rcpp_result_gen;
END_RCPP
}
// sample_dcsbm_cpp
arma::sp_mat sample_dcsbm_cpp(arma::uvec z, arma::mat Pmat, arma::vec theta);
RcppExport SEXP _nett_sample_dcsbm_cpp(SEXP zSEXP, SEXP PmatSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Pmat(PmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_dcsbm_cpp(z, Pmat, theta));
    return rcpp_result_gen;
END_RCPP
}
// sample_r1
arma::sp_mat sample_r1(arma::vec theta, double p);
RcppExport SEXP _nett_sample_r1(SEXP thetaSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_r1(theta, p));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nett_pair_dist2", (DL_FUNC) &_nett_pair_dist2, 2},
    {"_nett_pair_dist2_mat", (DL_FUNC) &_nett_pair_dist2_mat, 1},
    {"_nett_sample_dclvm_cpp", (DL_FUNC) &_nett_sample_dclvm_cpp, 3},
    {"_nett_sample_dcsbm_cpp", (DL_FUNC) &_nett_sample_dcsbm_cpp, 3},
    {"_nett_sample_r1", (DL_FUNC) &_nett_sample_r1, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_nett(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
