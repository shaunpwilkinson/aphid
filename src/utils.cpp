#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix progression(IntegerVector path, IntegerVector start){
  IntegerMatrix res(2, path.size());
  res(_, 0) = start;
  int xpos = start[0];
  int ypos = start[1];
  for(int i = 0; i < path.size() - 1; i ++){
    if(path[i] == 1){
      xpos++;
      ypos++;
    }else if(path[i] == 0){
      xpos++;
    }else if(path[i] == 2){
      ypos++;
    }else throw Rcpp::exception("path contains unknown elements");
    res (0, i + 1) = xpos;
    res (1, i + 1) = ypos;
  }
  return(res);
}

// [[Rcpp::export]]
IntegerMatrix progression2(IntegerVector path, IntegerVector start){
  IntegerMatrix res(2, path.size());
  res(_, 0) = start;
  int xpos = start[0];
  int ypos = start[1];
  for(int i = 0; i < path.size() - 1; i ++){
    if(path[i] == 2){
      xpos++;
      ypos++;
    }else if(path[i] < 2){
      xpos++;
    }else if(path[i] > 2){
      ypos++;
    }
    res (0, i + 1) = xpos;
    res (1, i + 1) = ypos;
  }
  return(res);
}


//' Count transition frequencies.
//'
//' Summation of transition frequencies in an integer vector.
//'
//' @param x an integer vector.
//' @param arity an integer representing the numbering system of the input vector,
//' 2 for binary, 3 for ternary, etc.
//'
// [[Rcpp::export]]
IntegerVector transitioncount(IntegerVector x, int arity){
  int newarity = pow(arity, 2);
  IntegerVector guide = seq(0, newarity - 1);
  guide.attr("dim") = IntegerVector::create(arity, arity);
  IntegerVector out(newarity);
  for(int i = 1; i < x.size(); i++) out[guide(x[i - 1], x[i])]++;
  out.attr("dim") = IntegerVector::create(arity, arity);
  return(out);
}


// [[Rcpp::export]]
IntegerVector emissioncount(IntegerVector states, int statearity,
                            IntegerVector residues, int resarity){
  int newarity = statearity * resarity;
  IntegerVector guide = seq(0, newarity - 1);
  guide.attr("dim") = IntegerVector::create(statearity, resarity);
  IntegerVector out(newarity);
  for(int i = 0; i < states.size(); i++) out[guide(states[i], residues[i])]++;
  out.attr("dim") = IntegerVector::create(statearity, resarity);
  return(out);
}


// elements of x can be 0, 1, 2, or NA
// number of modules should include the begin and end states
// outputs a 9 row count (integer) matrix with ncol = nmodules - 1 (doesn't include end)

// [[Rcpp::export]]
NumericMatrix tab9C(IntegerMatrix x, NumericVector seqweights){
  // x is a ternary matrix with ncol = phmm length + 2
  // length of seqweights should be same as nrow(x)
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
      if(x(i, j) == 0 & x(i, k) == 0) out(0, module) = out(0, module) + seqweights[i];
      else if(x(i, j) == 0 & x(i, k) == 1) out(1, module) = out(1, module) + seqweights[i];
      else if(x(i, j) == 0 & x(i, k) == 2) out(2, module) = out(2, module) + seqweights[i];
      else if(x(i, j) == 1 & x(i, k) == 0) out(3, module)= out(3, module) + seqweights[i];
      else if(x(i, j) == 1 & x(i, k) == 1) out(4, module)= out(4, module) + seqweights[i];
      else if(x(i, j) == 1 & x(i, k) == 2) out(5, module)= out(5, module) + seqweights[i];
      else if(x(i, j) == 2 & x(i, k) == 0) out(6, module)= out(6, module) + seqweights[i];
      else if(x(i, j) == 2 & x(i, k) == 1) out(7, module)= out(7, module) + seqweights[i];
      else if(x(i, j) == 2 & x(i, k) == 2) out(8, module)= out(8, module) + seqweights[i];
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





