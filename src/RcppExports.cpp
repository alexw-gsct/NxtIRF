// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// IRF_RLE_From_Cov
List IRF_RLE_From_Cov(std::string s_in, int strand);
RcppExport SEXP _rIRFinder_IRF_RLE_From_Cov(SEXP s_inSEXP, SEXP strandSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type s_in(s_inSEXP);
    Rcpp::traits::input_parameter< int >::type strand(strandSEXP);
    rcpp_result_gen = Rcpp::wrap(IRF_RLE_From_Cov(s_in, strand));
    return rcpp_result_gen;
END_RCPP
}
// IRF_gunzip
int IRF_gunzip(std::string s_in, std::string s_out);
RcppExport SEXP _rIRFinder_IRF_gunzip(SEXP s_inSEXP, SEXP s_outSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type s_in(s_inSEXP);
    Rcpp::traits::input_parameter< std::string >::type s_out(s_outSEXP);
    rcpp_result_gen = Rcpp::wrap(IRF_gunzip(s_in, s_out));
    return rcpp_result_gen;
END_RCPP
}
// IRF_main
int IRF_main(std::string bam_file, std::string reference_file, std::string output_file);
RcppExport SEXP _rIRFinder_IRF_main(SEXP bam_fileSEXP, SEXP reference_fileSEXP, SEXP output_fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bam_file(bam_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_file(reference_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type output_file(output_fileSEXP);
    rcpp_result_gen = Rcpp::wrap(IRF_main(bam_file, reference_file, output_file));
    return rcpp_result_gen;
END_RCPP
}
// IRF_SupplyMappaRegionReads
int IRF_SupplyMappaRegionReads(std::string genome_file, std::string region_file, std::string out_fa, int read_len, int read_stride, int error_pos);
RcppExport SEXP _rIRFinder_IRF_SupplyMappaRegionReads(SEXP genome_fileSEXP, SEXP region_fileSEXP, SEXP out_faSEXP, SEXP read_lenSEXP, SEXP read_strideSEXP, SEXP error_posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type genome_file(genome_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type region_file(region_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type out_fa(out_faSEXP);
    Rcpp::traits::input_parameter< int >::type read_len(read_lenSEXP);
    Rcpp::traits::input_parameter< int >::type read_stride(read_strideSEXP);
    Rcpp::traits::input_parameter< int >::type error_pos(error_posSEXP);
    rcpp_result_gen = Rcpp::wrap(IRF_SupplyMappaRegionReads(genome_file, region_file, out_fa, read_len, read_stride, error_pos));
    return rcpp_result_gen;
END_RCPP
}
// IRF_SupplyMappaReads
int IRF_SupplyMappaReads(std::string genome_file, std::string out_fa, int read_len, int read_stride, int error_pos);
RcppExport SEXP _rIRFinder_IRF_SupplyMappaReads(SEXP genome_fileSEXP, SEXP out_faSEXP, SEXP read_lenSEXP, SEXP read_strideSEXP, SEXP error_posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type genome_file(genome_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type out_fa(out_faSEXP);
    Rcpp::traits::input_parameter< int >::type read_len(read_lenSEXP);
    Rcpp::traits::input_parameter< int >::type read_stride(read_strideSEXP);
    Rcpp::traits::input_parameter< int >::type error_pos(error_posSEXP);
    rcpp_result_gen = Rcpp::wrap(IRF_SupplyMappaReads(genome_file, out_fa, read_len, read_stride, error_pos));
    return rcpp_result_gen;
END_RCPP
}
// IRF_genmap
int IRF_genmap(std::string bam_file, std::string output_path);
RcppExport SEXP _rIRFinder_IRF_genmap(SEXP bam_fileSEXP, SEXP output_pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bam_file(bam_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type output_path(output_pathSEXP);
    rcpp_result_gen = Rcpp::wrap(IRF_genmap(bam_file, output_path));
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
    {"_rIRFinder_IRF_RLE_From_Cov", (DL_FUNC) &_rIRFinder_IRF_RLE_From_Cov, 2},
    {"_rIRFinder_IRF_gunzip", (DL_FUNC) &_rIRFinder_IRF_gunzip, 2},
    {"_rIRFinder_IRF_main", (DL_FUNC) &_rIRFinder_IRF_main, 3},
    {"_rIRFinder_IRF_SupplyMappaRegionReads", (DL_FUNC) &_rIRFinder_IRF_SupplyMappaRegionReads, 6},
    {"_rIRFinder_IRF_SupplyMappaReads", (DL_FUNC) &_rIRFinder_IRF_SupplyMappaReads, 5},
    {"_rIRFinder_IRF_genmap", (DL_FUNC) &_rIRFinder_IRF_genmap, 2},
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
