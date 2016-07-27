#include <Rcpp.h>
using namespace Rcpp;
//' Sum of logged values
//'
//' \code{logsum} takes a vector of logged values and returns the log of the sum of their
//' exponentials.
//'
//' @param x a numeric vector of logged values (usually probabilities) to be summed over.
//' @return a numerically stable calculation of \code{log(sum(exp(x)))}.
//' @export
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


