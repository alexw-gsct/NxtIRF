// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// test_IRF
int test_IRF(std::string bam_file, std::string reference_path, std::string output_path);
RcppExport SEXP _rIRFinder_test_IRF(SEXP bam_fileSEXP, SEXP reference_pathSEXP, SEXP output_pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bam_file(bam_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_path(reference_pathSEXP);
    Rcpp::traits::input_parameter< std::string >::type output_path(output_pathSEXP);
    rcpp_result_gen = Rcpp::wrap(test_IRF(bam_file, reference_path, output_path));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _rIRFinder_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _rIRFinder_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _rIRFinder_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _rIRFinder_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rIRFinder_test_IRF", (DL_FUNC) &_rIRFinder_test_IRF, 3},
    {"_rIRFinder_rcpparma_hello_world", (DL_FUNC) &_rIRFinder_rcpparma_hello_world, 0},
    {"_rIRFinder_rcpparma_outerproduct", (DL_FUNC) &_rIRFinder_rcpparma_outerproduct, 1},
    {"_rIRFinder_rcpparma_innerproduct", (DL_FUNC) &_rIRFinder_rcpparma_innerproduct, 1},
    {"_rIRFinder_rcpparma_bothproducts", (DL_FUNC) &_rIRFinder_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_rIRFinder(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
