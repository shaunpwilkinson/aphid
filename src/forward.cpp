#include <Rcpp.h>
using namespace Rcpp;
//' Full log probability of sequence given model.
//'
//' Implementation of the forward algorithm to caluclate the full (log) probability
//' of a sequence through a given a HMM or profile HMM.
//' @param x an object of class \code{HMM} or \code{PHMM}.
//'
// [[Rcpp::export]]
double logsum(NumericVector x) {
  int n = x.size();
  double res = x[0];
  if(n == 1) return(res);
  for(int i = 1; i < n; i++){
    if(res == -INFINITY) res = x[i];
    else if(x[i] == -INFINITY) res += 0;
    else if(res > x[i]) res = res + log1p(exp(x[i] - res));
    else res = x[i] + log1p(exp(res - x[i]));
  }
  return(res);
}

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
    List out = List::create(Named("score") = A(1, 1), Named("path") = NULL);
    out.attr("class") = "forward";
    return out;
  }
  IntegerVector yind(nrolls);
  CharacterVector yi(1);
  for(int i = 0; i < nrolls; i++){
    yi = y[i];
    yind[i] = match(yi, residues)[0] - 1;
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
  NumericVector coltotals(nstates);
  for(int i = 1; i < nrolls; i++){
    for(int k = 0; k < nstates; k++){
      for(int l = 0; l < nstates; l++){
        tmp(k, l) = R(k, i - 1) + A(k + 1, l + 1);
      }
    }
    for(int k = 0; k < nstates; k++) coltotals[k] = logsum(tmp(_, k));
    R(_, i) =  E(_, yind[i]) + coltotals;
  }
  NumericVector ak0(nstates);
  if(any(e[seq(1, nstates)] != rep(-INFINITY, nstates)).is_true()) ak0 = e[seq(1, nstates)];
  double res = logsum(R(_, nrolls - 1) + ak0);
  List out = List::create(Named("score") = res, Named("array") = R);
  out.attr("class") = "forward";
  return out;
}




// [[Rcpp::export]]
List backwardC(List x, CharacterVector y, bool logspace = false) {
  NumericMatrix A = clone(VECTOR_ELT(x, 0));
  NumericMatrix E = clone(VECTOR_ELT(x, 1));
  List names = E.attr("dimnames");
  CharacterVector states = VECTOR_ELT(names, 0);
  CharacterVector residues = VECTOR_ELT(names, 1);
  int nrolls = y.size();
  int nstates = states.size(); // does not include begin/end state
  int nres = residues.size();
  if(!logspace){
    for(int i = 0; i < nstates; i++){
      A(i, _) = log(A(i, _));
      E(i, _) = log(E(i, _));
    }
    A(nstates, _) = log(A(nstates, _));
  }
  if(nrolls == 0){
    List out = List::create(Named("score") = A(1, 1), Named("path") = NULL);
    out.attr("class") = "forward";
    return out;
  }
  IntegerVector yind(nrolls);
  CharacterVector yi(1);
  for(int i = 0; i < nrolls; i++){
    yi = y[i];
    yind[i] = match(yi, residues)[0] - 1;
  }
  NumericMatrix R(nstates, nrolls);
  //IntegerVector prerolls = Range(nrolls, 1); // doesn't work in reverse
  IntegerVector prerolls(nrolls);
  prerolls[0] = nrolls;
  for(int i = 1; i < nrolls; i++) prerolls[i] = prerolls[i - 1] - 1;
  CharacterVector rolls = as<CharacterVector>(prerolls);
  List rnames = List::create(Named("state") = states, Named("roll") = rolls);
  R.attr("dimnames") = rnames;
  NumericVector s = A(0, _);
  NumericVector e = A(_, 0);
  NumericVector ak0(nstates);
  if(any(e[seq(1, nstates)] != rep(-INFINITY, nstates)).is_true()) ak0 = e[seq(1, nstates)];
  R(_, nrolls - 1) = ak0;
  NumericMatrix tmp(nstates, nstates);
  NumericVector rowtotals(nstates);
  for(int i = nrolls - 1; i > 0; i--){
    for(int k = 0; k < nstates; k++){
      for(int l = 0; l < nstates; l++){
        tmp(k, l) =  A(k + 1, l + 1) + E(l, yind[i]) + R(l, i);
      }
    }
    for(int k = 0; k < nstates; k++) rowtotals[k] = logsum(tmp(k, _));
    R(_, i - 1) =  rowtotals;
  }
  NumericVector logprobs(nstates);
  for(int l = 0; l < nstates; l++) logprobs[l] = A(0, l + 1) + E(l, yind[0]) + R(l, 0);
  double res = logsum(logprobs);
  List out = List::create(Named("score") = res, Named("array") = R);
  out.attr("class") = "backward";
  return out;
}
