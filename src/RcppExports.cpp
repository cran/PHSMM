// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// nll_Rcpp
double nll_Rcpp(arma::mat allprobs, arma::mat gamma, arma::rowvec delta, int T_y);
RcppExport SEXP _PHSMM_nll_Rcpp(SEXP allprobsSEXP, SEXP gammaSEXP, SEXP deltaSEXP, SEXP T_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type allprobs(allprobsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type T_y(T_ySEXP);
    rcpp_result_gen = Rcpp::wrap(nll_Rcpp(allprobs, gamma, delta, T_y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PHSMM_nll_Rcpp", (DL_FUNC) &_PHSMM_nll_Rcpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_PHSMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
