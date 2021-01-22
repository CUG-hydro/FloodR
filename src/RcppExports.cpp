// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// leastSqrRegression
List leastSqrRegression(NumericVector x_in, NumericVector y_in, IntegerMatrix Comb, bool intercept_zero);
RcppExport SEXP _FloodR_leastSqrRegression(SEXP x_inSEXP, SEXP y_inSEXP, SEXP CombSEXP, SEXP intercept_zeroSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x_in(x_inSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y_in(y_inSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type Comb(CombSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept_zero(intercept_zeroSEXP);
    rcpp_result_gen = Rcpp::wrap(leastSqrRegression(x_in, y_in, Comb, intercept_zero));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FloodR_leastSqrRegression", (DL_FUNC) &_FloodR_leastSqrRegression, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_FloodR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}