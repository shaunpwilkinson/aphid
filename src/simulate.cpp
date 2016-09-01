#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List returnmod(List x) {
  NumericVector A = clone(VECTOR_ELT(x, 0));
  NumericVector E = clone(VECTOR_ELT(x, 1));
  NumericVector qa = clone(VECTOR_ELT(x, 2));
  NumericVector qe = clone(VECTOR_ELT(x, 3));
  NumericVector dimA = A.attr("dim");
  int tosize = dimA[0] * dimA[1];
  IntegerVector todims = IntegerVector::create(dimA[0], dimA[1]);
  IntegerVector toDindices = seq(0, (tosize - 1));
  IntegerVector toMindices = seq(tosize, (2 * tosize - 1));
  IntegerVector toIindices = seq(2 * tosize, (3 * tosize - 1));
  NumericVector toD = A[toDindices];
  NumericVector toM = A[toMindices];
  NumericVector toI = A[toIindices];
  toD.attr("dim") = todims;
  toM.attr("dim") = todims;
  toI.attr("dim") = todims;
  List outls(5);
  int position = 0;
  for(int i = 0; i < 5; i++){
    outls[i] = qa;
  }
  return outls;
}

