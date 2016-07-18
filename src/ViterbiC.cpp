#include <Rcpp.h>
using namespace Rcpp;
//' Find the optimal path through a HMM or PHMM
//'
//' \code{ViterbiC} finds the optimal path of a sequence through a HMM
//' or PHMM and returns its log-odds score.
//'
//' @param x an object of class \code{HMM} or \code{PHMM}, or a character vector.
//' @param y a character vector consisting of residues emitted by the
//' HMM or PHMM.
//' @param logspace logical argument indicating whether the emission
//' and transmission probabilities for the model(s) are logged.
//' @name ViterbiC
//' @export
// [[Rcpp::export(name = "ViterbiC.HMM")]]
List ViterbiHMM(List x, CharacterVector y, bool logspace = false){
  NumericVector s = clone(VECTOR_ELT(x, 0));
  NumericMatrix A = clone(VECTOR_ELT(x, 1));
  NumericMatrix E = clone(VECTOR_ELT(x, 2));
  List names = E.attr("dimnames");
  CharacterVector states = VECTOR_ELT(names, 0);
  CharacterVector residues = VECTOR_ELT(names, 1);
  int nstates = E.nrow();
  int nrolls = y.size();
  int nresidues = residues.size();
  if(!logspace){
    for(int i = 0; i < nstates; i++){
      s[i] = log(s[i]);
      A(i, _) = log(A(i, _));
      E(i, _) = log(E(i, _));
    }
  }
  IntegerVector yind(nrolls);
  //IntegerVector rseq = seq_along(residues); // this starts at 1
  //IntegerVector rseq = Range(0, nresidues - 1);
  for(int i = 0; i < nrolls; i++){
    yind[i] = match(CharacterVector::create(y[i]), residues)[0] - 1;
  }
  NumericMatrix V(nstates, nrolls);
  IntegerVector prerolls = Range(1, nrolls);
  CharacterVector rolls = as<CharacterVector>(prerolls);
  List vnames = List::create(Named("state") = states, Named("roll") = rolls);
  V.attr("dimnames") = vnames;
  V(_, 0) = E(_, yind[0]) + s;
  IntegerMatrix rawP(nstates, nrolls);
  //CharacterMatrix P(nstates, nrolls);
  //P.attr("dimnames") = vnames;
  NumericMatrix tmp(nstates, nstates);
  NumericVector colmaxs(nstates);
  IntegerVector maxstates(nstates);
  for(int i = 1; i < nrolls; i++){
    for(int j = 0; j < nstates; j++){
      for(int k = 0; k < nstates; k++){
        tmp(j, k) = V(j, i - 1) + A(j, k);
      }
    }
    for(int l = 0; l < nstates; l++){
      int maxstate = which_max(tmp(_, l));
      colmaxs[l] = tmp(_, l)[maxstate];
      maxstates[l] = maxstate;
    }
    V(_, i) = E(_, yind[i]) + colmaxs;
    rawP(_, i) = maxstates;
    //P(_, i) = as<CharacterVector>(states[maxstates]);
    checkUserInterrupt(); // could speed up by checking less reqularly
  }
  int maxstate = which_max(V(_, nrolls - 1));
  double score = V(_, nrolls - 1)[maxstate];
  IntegerVector rawpath(nrolls);
  CharacterVector path(nrolls);
  rawpath[nrolls - 1] = maxstate;
  path[nrolls - 1] = states[maxstate];
  for(int i = nrolls - 2; i >= 0; i--){
    rawpath[i] = rawP(maxstate, i + 1);
    path[i] = states[rawpath[i]];
    maxstate = rawpath[i];
  }
  List out = List::create(//Named("V") = V,
    //Named("P") = P,
    Named("score") = score,
    Named("path") = path);
  out.attr("class") = "Viterbi";
  return out;
}

//' @rdname ViterbiC
//' @export
// [[Rcpp::export(name = "ViterbiC.PHMM")]]
List ViterbiPHMM(List x, CharacterVector y, bool logspace = false){
  Rcout << "sorry Viterbi method is not available for PHMMs yet";
  return 0;
}
