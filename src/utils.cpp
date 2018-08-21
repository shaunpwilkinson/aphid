#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(name = ".acount")]]
IntegerVector acount(IntegerVector x, int arity){
  // x is an integer vector in a certain coding scheme (starting at 0)
  // arity is an integer representing the numbering system of the input vector,
  // 2 for binary, 3 for ternary, 4 for quarternary, etc.
  // returns the transition frequencies as an integer vector
  // int newarity = pow(arity, 2);
  int newarity = arity * arity;
  IntegerVector guide = seq(0, newarity - 1);
  guide.attr("dim") = IntegerVector::create(arity, arity);
  IntegerVector out(newarity);
  for(int i = 1; i < x.size(); i++) out[guide(x[i - 1], x[i])]++;
  out.attr("dim") = IntegerVector::create(arity, arity);
  return(out);
}

// [[Rcpp::export(name = ".ecount")]]
IntegerVector ecount(IntegerVector states, int statearity,
                            IntegerVector residues, int resarity){
  int newarity = statearity * resarity;
  IntegerVector guide = seq(0, newarity - 1);
  guide.attr("dim") = IntegerVector::create(statearity, resarity);
  IntegerVector out(newarity);
  for(int i = 0; i < states.size(); i++) out[guide(states[i], residues[i])]++;
  out.attr("dim") = IntegerVector::create(statearity, resarity);
  return(out);
}

// [[Rcpp::export(name = ".atab")]]
NumericMatrix atab(IntegerMatrix x, NumericVector seqweights){
  // tabulates the 9 transition types found in PHMMs
  // x is a ternary matrix with ncol = phmm length + 2
  // length of seqweights should be same as nrow(x)
  // elements of x can be 0, 1, 2, or NA
  // number of modules should include the begin and end states
  // outputs a 9 row count (integer) matrix with ncol =
  // nmodules - 1 (doesn't include end)
  if(x.nrow() != seqweights.size()){
    throw Rcpp::exception("length of seqweights vector should equal number of sequences");
  }
  int modules = 0;
  for(int i = 0; i < x.ncol(); i++) if((x(0, i) == 0) | (x(0, i) == 1)) modules++;
  // modules is the total number of modules in the model including begin and end states
  NumericMatrix out(9, modules - 1);
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
      while(IntegerMatrix::is_na(x(i, k))) k++;
      if((x(i, j) == 0) & (x(i, k) == 0)) out(0, module) = out(0, module) + seqweights[i];
      else if((x(i, j) == 0) & (x(i, k) == 1)) out(1, module) = out(1, module) + seqweights[i];
      else if((x(i, j) == 0) & (x(i, k) == 2)) out(2, module) = out(2, module) + seqweights[i];
      else if((x(i, j) == 1) & (x(i, k) == 0)) out(3, module) = out(3, module) + seqweights[i];
      else if((x(i, j) == 1) & (x(i, k) == 1)) out(4, module) = out(4, module) + seqweights[i];
      else if((x(i, j) == 1) & (x(i, k) == 2)) out(5, module) = out(5, module) + seqweights[i];
      else if((x(i, j) == 2) & (x(i, k) == 0)) out(6, module) = out(6, module) + seqweights[i];
      else if((x(i, j) == 2) & (x(i, k) == 1)) out(7, module) = out(7, module) + seqweights[i];
      else if((x(i, j) == 2) & (x(i, k) == 2)) out(8, module) = out(8, module) + seqweights[i];
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
