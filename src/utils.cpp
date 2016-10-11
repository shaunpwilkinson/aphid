#include <Rcpp.h>
using namespace Rcpp;


// this is actually slower than the R version for small numbers of residues
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


// [[Rcpp::export]]
double DNAprobC(RawVector a, NumericVector probs){
  // a is a raw byte in format of Paradis (2007)
  // probs is a 4-element numeric vector of probabilities for the set {a,c,g,t}
  if(probs.size() != 4){
    throw Rcpp::exception("probs argument must be a numeric vector of length 4");
  }
  RawVector tmp = RawVector::create(4, 55, 0, 136, 72, 199, 40, 24,
                                    224, 176, 208, 112, 240, 8, 160,
                                    144, 96, 80);
  if(a[0] <= tmp[0]){
    throw Rcpp::exception("Input sequence contains gaps or unknown characters");
  }else{
    if((a[0] & tmp[1]) == tmp[2]){ // is purine?
      if(a[0] == tmp[3]){
        return(probs[0]); // A
      } else if(a[0] == tmp[4]){
        return(probs[2]); // G
      } else{
        return(log((exp(probs[0]) + exp(probs[2]))/2)); // R (A or G)
      }
    } else if((a[0] & tmp[5]) == tmp[2]){ // is pyrimidine?
      if(a[0] == tmp[6]){
        return(probs[1]); // C
      } else if(a[0] == tmp[7]){
        return(probs[3]); // T
      } else{
        return(log((exp(probs[1]) + exp(probs[3]))/2)); // Y (C or T)
      }
    }else if(a[0] == tmp[14]){
      return(log((exp(probs[0]) + exp(probs[1]))/2)); // M (A or C)
    }else if(a[0] == tmp[15]){
      return(log((exp(probs[0]) + exp(probs[3]))/2)); // W (A or T)
    }else if(a[0] == tmp[16]){
      return(log((exp(probs[1]) + exp(probs[2]))/2)); // S (G or C)
    }else if(a[0] == tmp[17]){
      return(log((exp(probs[2]) + exp(probs[3]))/2)); // K (G or T)
    }else if(a[0] == tmp[8]){
      return(log((exp(probs[0]) + exp(probs[1]) + exp(probs[2]))/3));
      // V (A or C or G)
    } else if(a[0] == tmp[9]){
      return(log((exp(probs[0]) + exp(probs[1]) + exp(probs[3]))/3));
      // H (A or C or T)
    } else if(a[0] == tmp[10]){
      return(log((exp(probs[0]) + exp(probs[2]) + exp(probs[3]))/3));
      // D (A or G or T)
    } else if(a[0] == tmp[11]){
      return(log((exp(probs[1]) + exp(probs[2]) + exp(probs[3]))/3));
      // B (C or G or T)
    } else if(a[0] == tmp[12]){
      return(log((exp(probs[0]) + exp(probs[1]) +
             exp(probs[2]) + exp(probs[3]))/4));   // N
    } else {
      throw Rcpp::exception("Invalid byte");
    }
  }
}

