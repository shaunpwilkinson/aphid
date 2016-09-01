#include <Rcpp.h>
using namespace Rcpp;


// this is actually slower than the r version for small numbers of residues
// IntegerVector tab(CharacterVector v, CharacterVector residues){
//   IntegerVector out(residues.size());
//   for(int i = 0; i < residues.size(); i++){
//     for(int j = 0; j < v.size(); j++){
//       if(v[j] == residues[i]) out[i]++;
//     }
//   }
//   out.attr("names") = residues;
//   return(out);
// }





// elements of x can be 0, 1, 2, or NA
// number of modules should include the begin and end states
// outputs a 9 row count (integer) matrix with ncol = nmodules - 1 (doesn't include end)

// [[Rcpp::export]]
IntegerMatrix tab9C(IntegerMatrix x, int modules) {
  IntegerMatrix out(9, modules - 1);
  CharacterVector transtype = CharacterVector::create("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II");
  IntegerVector tmp = seq(0, modules - 2);
  CharacterVector frommod = as<CharacterVector>(tmp);
  List outnames = List::create(Named("type") = transtype, Named("module") = frommod);
  out.attr("dimnames") = outnames;
  int j = 0;
  int k = 1;
  int module = 0; // from module
  for(int i = 0; i < x.nrow(); i++){
    while(j < x.ncol() - 1){
      while(IntegerMatrix::is_na(x(i, k))){
        k++;
      }
      if(x(i, j) == 0 & x(i, k) == 0) out(0, module)++;
      else if(x(i, j) == 0 & x(i, k) == 1) out(1, module)++;
      else if(x(i, j) == 0 & x(i, k) == 2) out(2, module)++;
      else if(x(i, j) == 1 & x(i, k) == 0) out(3, module)++;
      else if(x(i, j) == 1 & x(i, k) == 1) out(4, module)++;
      else if(x(i, j) == 1 & x(i, k) == 2) out(5, module)++;
      else if(x(i, j) == 2 & x(i, k) == 0) out(6, module)++;
      else if(x(i, j) == 2 & x(i, k) == 1) out(7, module)++;
      else if(x(i, j) == 2 & x(i, k) == 2) out(8, module)++;
      //
      if(x(i, k) != 2) module++;
      j = k;
      k++;
      checkUserInterrupt();
    }
    j = 0;
    k = 1;
    module = 0;
  }
  return out;
}
