#include <Rcpp.h>
using namespace Rcpp;
//' Full log probability of sequence given model.
//'
//' Implementation of the forward algorithm to caluclate the full (log) probability
//' of a sequence through a given a HMM or profile HMM.
//' @param x an object of class \code{HMM} or \code{PHMM}.
//'
// [[Rcpp::export]]
List forwardC(List x, CharacterVector y, bool logspace = false) {
  NumericMatrix A = clone(VECTOR_ELT(x, 0));
  NumericMatrix E = clone(VECTOR_ELT(x, 1));
  List names = E.attr("dimnames");
  CharacterVector states = VECTOR_ELT(names, 0);
  CharacterVector residues = VECTOR_ELT(names, 1);
  int nrolls = y.size();
  int nstates = states.size();
  int nres = residues.size();
  if(!logspace){
    for(int i = 0; i < nstates; i++){
      A(i, _) = log(A(i, _));
      E(i, _) = log(E(i, _));
    }
    A(nstates, _) = log(A(nstates, _));
  }
  if(nrolls == 0){
    List out = List::create(
      Named("score") = A(1, 1),
      Named("path") = NULL);
    out.attr("class") = "forward";
    return out;
  }
  IntegerVector yind(nrolls);
  for(int i = 0; i < nrolls; i++){
    yind[i] = match(CharacterVector::create(y[i]), residues)[0] - 1;
  }
  NumericMatrix R(nstates, nrolls);
  IntegerVector prerolls = Range(1, nrolls);
  CharacterVector rolls = as<CharacterVector>(prerolls);
  List rnames = List::create(Named("state") = states, Named("roll") = rolls);
  R.attr("dimnames") = rnames;
  NumericVector s = A(0, _);
  NumericVector e = A(_, 0);
  R(_, 0) = E(_, yind[0]) + s[seq(1, nstates)];
  NumericMatrix tmp(nstates, nstates);
  NumericVector tmpk = tmp(_, 0);
  NumericVector coltotals(nstates);
  NumericVector newcol(nstates);
  double res = 0;
  for(int i = 1; i < nrolls; i++){
    for(int k = 0; k < nstates; k++){
      for(int l = 0; l < nstates; l++){
        tmp(k, l) = R(k, i - 1) + A(k + 1, l + 1);
      }
    }
    for(int k = 0; k < nstates; k++){
      tmpk = tmp(_, k);
      res = tmpk[0];
      for(int l = 1; l < nstates; l++){
        if(res == -INFINITY) res = tmpk[l];
        else if(tmpk[l] == -INFINITY) res += 0;
        else if(res > tmpk[l]) res = res + log1p(exp(tmpk[l] - res));
        else res = tmpk[l] + log1p(exp(res - tmpk[l]));
      }
      coltotals[k] = res;
    }
    newcol = E(_, yind[i]) + coltotals;
    //coltotals = rep(0, nstates);
    R(_, i) = newcol;
  }
  NumericVector ak0(nstates);
  if(any(e[seq(1, nstates)] != rep(-INFINITY, nstates - 1)).is_true()) ak0 = e[seq(1, nstates)];
  tmpk = R(_, nrolls - 1) + ak0;
  res = tmpk[0];
  for(int l = 1; l < nstates; l++){
    if(res == -INFINITY) res = tmpk[l];
    else if(tmpk[l] == -INFINITY) res += 0;
    else if(res > tmpk[l]) res = res + log1p(exp(tmpk[l] - res));
    else res = tmpk[l] + log1p(exp(res - tmpk[l]));
  }
  List out = List::create(Named("score") = res, Named("array") = R);
  out.attr("class") = "forward";
  return out;
}
