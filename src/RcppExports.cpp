// Generated manually for postNet vendored G4 search.

#include <Rcpp.h>

Rcpp::List find_gquadruplexes_cpp(SEXP subject, int min_score);

RcppExport SEXP _postNet_find_gquadruplexes_cpp(SEXP subjectSEXP, SEXP min_scoreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type subject(subjectSEXP);
    Rcpp::traits::input_parameter< int >::type min_score(min_scoreSEXP);
    rcpp_result_gen = Rcpp::wrap(find_gquadruplexes_cpp(subject, min_score));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_postNet_find_gquadruplexes_cpp", (DL_FUNC) &_postNet_find_gquadruplexes_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_postNet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
