// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// kcountDNA
NumericMatrix kcountDNA(List x, int k);
RcppExport SEXP aphid_kcountDNA(SEXP xSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(kcountDNA(x, k));
    return rcpp_result_gen;
END_RCPP
}
// kdist
NumericMatrix kdist(NumericMatrix x, IntegerVector from, IntegerVector to, IntegerVector seqlengths, int k);
RcppExport SEXP aphid_kdist(SEXP xSEXP, SEXP fromSEXP, SEXP toSEXP, SEXP seqlengthsSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type seqlengths(seqlengthsSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(kdist(x, from, to, seqlengths, k));
    return rcpp_result_gen;
END_RCPP
}
// map
LogicalVector map(NumericMatrix ecs, LogicalMatrix notgaps, List pseudocounts, NumericVector seqweights, NumericVector qe, double lambda);
RcppExport SEXP aphid_map(SEXP ecsSEXP, SEXP notgapsSEXP, SEXP pseudocountsSEXP, SEXP seqweightsSEXP, SEXP qeSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ecs(ecsSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type notgaps(notgapsSEXP);
    Rcpp::traits::input_parameter< List >::type pseudocounts(pseudocountsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type seqweights(seqweightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qe(qeSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(map(ecs, notgaps, pseudocounts, seqweights, qe, lambda));
    return rcpp_result_gen;
END_RCPP
}
// acount
IntegerVector acount(IntegerVector x, int arity);
RcppExport SEXP aphid_acount(SEXP xSEXP, SEXP aritySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type arity(aritySEXP);
    rcpp_result_gen = Rcpp::wrap(acount(x, arity));
    return rcpp_result_gen;
END_RCPP
}
// ecount
IntegerVector ecount(IntegerVector states, int statearity, IntegerVector residues, int resarity);
RcppExport SEXP aphid_ecount(SEXP statesSEXP, SEXP statearitySEXP, SEXP residuesSEXP, SEXP resaritySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type states(statesSEXP);
    Rcpp::traits::input_parameter< int >::type statearity(statearitySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type residues(residuesSEXP);
    Rcpp::traits::input_parameter< int >::type resarity(resaritySEXP);
    rcpp_result_gen = Rcpp::wrap(ecount(states, statearity, residues, resarity));
    return rcpp_result_gen;
END_RCPP
}
// atab
NumericMatrix atab(IntegerMatrix x, NumericVector seqweights);
RcppExport SEXP aphid_atab(SEXP xSEXP, SEXP seqweightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type seqweights(seqweightsSEXP);
    rcpp_result_gen = Rcpp::wrap(atab(x, seqweights));
    return rcpp_result_gen;
END_RCPP
}
// logsum
double logsum(NumericVector x);
RcppExport SEXP aphid_logsum(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logsum(x));
    return rcpp_result_gen;
END_RCPP
}
// whichmax
int whichmax(NumericVector x, int start);
RcppExport SEXP aphid_whichmax(SEXP xSEXP, SEXP startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type start(startSEXP);
    rcpp_result_gen = Rcpp::wrap(whichmax(x, start));
    return rcpp_result_gen;
END_RCPP
}
// probDNA
double probDNA(int x, NumericVector probs);
RcppExport SEXP aphid_probDNA(SEXP xSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(probDNA(x, probs));
    return rcpp_result_gen;
END_RCPP
}
// probAA
double probAA(int x, NumericVector probs);
RcppExport SEXP aphid_probAA(SEXP xSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(probAA(x, probs));
    return rcpp_result_gen;
END_RCPP
}
// ViterbiD
List ViterbiD(IntegerVector x, IntegerVector y, int type, double d, double e, NumericMatrix S, IntegerVector windowspace, double offset);
RcppExport SEXP aphid_ViterbiD(SEXP xSEXP, SEXP ySEXP, SEXP typeSEXP, SEXP dSEXP, SEXP eSEXP, SEXP SSEXP, SEXP windowspaceSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type windowspace(windowspaceSEXP);
    Rcpp::traits::input_parameter< double >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(ViterbiD(x, y, type, d, e, S, windowspace, offset));
    return rcpp_result_gen;
END_RCPP
}
// ViterbiH
List ViterbiH(IntegerVector y, NumericMatrix A, NumericMatrix E, bool DNA, bool AA);
RcppExport SEXP aphid_ViterbiH(SEXP ySEXP, SEXP ASEXP, SEXP ESEXP, SEXP DNASEXP, SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type E(ESEXP);
    Rcpp::traits::input_parameter< bool >::type DNA(DNASEXP);
    Rcpp::traits::input_parameter< bool >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(ViterbiH(y, A, E, DNA, AA));
    return rcpp_result_gen;
END_RCPP
}
// ViterbiP
List ViterbiP(IntegerVector y, NumericMatrix A, NumericMatrix E, NumericVector qe, NumericVector qey, int type, IntegerVector windowspace, double offset, bool DI, bool ID, bool DNA, bool AA);
RcppExport SEXP aphid_ViterbiP(SEXP ySEXP, SEXP ASEXP, SEXP ESEXP, SEXP qeSEXP, SEXP qeySEXP, SEXP typeSEXP, SEXP windowspaceSEXP, SEXP offsetSEXP, SEXP DISEXP, SEXP IDSEXP, SEXP DNASEXP, SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type E(ESEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qe(qeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qey(qeySEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type windowspace(windowspaceSEXP);
    Rcpp::traits::input_parameter< double >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< bool >::type DI(DISEXP);
    Rcpp::traits::input_parameter< bool >::type ID(IDSEXP);
    Rcpp::traits::input_parameter< bool >::type DNA(DNASEXP);
    Rcpp::traits::input_parameter< bool >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(ViterbiP(y, A, E, qe, qey, type, windowspace, offset, DI, ID, DNA, AA));
    return rcpp_result_gen;
END_RCPP
}
// ViterbiPP
List ViterbiPP(NumericMatrix Ax, NumericMatrix Ay, NumericMatrix Ex, NumericMatrix Ey, NumericVector qe, int type, IntegerVector windowspace, double offset);
RcppExport SEXP aphid_ViterbiPP(SEXP AxSEXP, SEXP AySEXP, SEXP ExSEXP, SEXP EySEXP, SEXP qeSEXP, SEXP typeSEXP, SEXP windowspaceSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Ax(AxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ay(AySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ex(ExSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ey(EySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qe(qeSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type windowspace(windowspaceSEXP);
    Rcpp::traits::input_parameter< double >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(ViterbiPP(Ax, Ay, Ex, Ey, qe, type, windowspace, offset));
    return rcpp_result_gen;
END_RCPP
}
// forwardH
List forwardH(IntegerVector y, NumericMatrix A, NumericMatrix E, bool DNA, bool AA);
RcppExport SEXP aphid_forwardH(SEXP ySEXP, SEXP ASEXP, SEXP ESEXP, SEXP DNASEXP, SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type E(ESEXP);
    Rcpp::traits::input_parameter< bool >::type DNA(DNASEXP);
    Rcpp::traits::input_parameter< bool >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(forwardH(y, A, E, DNA, AA));
    return rcpp_result_gen;
END_RCPP
}
// forwardP
List forwardP(IntegerVector y, NumericMatrix A, NumericMatrix E, NumericVector qe, NumericVector qey, int type, IntegerVector windowspace, bool DI, bool ID, bool DNA, bool AA);
RcppExport SEXP aphid_forwardP(SEXP ySEXP, SEXP ASEXP, SEXP ESEXP, SEXP qeSEXP, SEXP qeySEXP, SEXP typeSEXP, SEXP windowspaceSEXP, SEXP DISEXP, SEXP IDSEXP, SEXP DNASEXP, SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type E(ESEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qe(qeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qey(qeySEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type windowspace(windowspaceSEXP);
    Rcpp::traits::input_parameter< bool >::type DI(DISEXP);
    Rcpp::traits::input_parameter< bool >::type ID(IDSEXP);
    Rcpp::traits::input_parameter< bool >::type DNA(DNASEXP);
    Rcpp::traits::input_parameter< bool >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(forwardP(y, A, E, qe, qey, type, windowspace, DI, ID, DNA, AA));
    return rcpp_result_gen;
END_RCPP
}
// backwardH
List backwardH(IntegerVector y, NumericMatrix A, NumericMatrix E, bool DNA, bool AA);
RcppExport SEXP aphid_backwardH(SEXP ySEXP, SEXP ASEXP, SEXP ESEXP, SEXP DNASEXP, SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type E(ESEXP);
    Rcpp::traits::input_parameter< bool >::type DNA(DNASEXP);
    Rcpp::traits::input_parameter< bool >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(backwardH(y, A, E, DNA, AA));
    return rcpp_result_gen;
END_RCPP
}
// backwardP
List backwardP(IntegerVector y, NumericMatrix A, NumericMatrix E, NumericVector qe, NumericVector qey, int type, IntegerVector windowspace, bool DI, bool ID, bool DNA, bool AA);
RcppExport SEXP aphid_backwardP(SEXP ySEXP, SEXP ASEXP, SEXP ESEXP, SEXP qeSEXP, SEXP qeySEXP, SEXP typeSEXP, SEXP windowspaceSEXP, SEXP DISEXP, SEXP IDSEXP, SEXP DNASEXP, SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type E(ESEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qe(qeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qey(qeySEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type windowspace(windowspaceSEXP);
    Rcpp::traits::input_parameter< bool >::type DI(DISEXP);
    Rcpp::traits::input_parameter< bool >::type ID(IDSEXP);
    Rcpp::traits::input_parameter< bool >::type DNA(DNASEXP);
    Rcpp::traits::input_parameter< bool >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(backwardP(y, A, E, qe, qey, type, windowspace, DI, ID, DNA, AA));
    return rcpp_result_gen;
END_RCPP
}
