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
//' @param numbersystem an integer representing the numbering system of the input vector,
//' 2 for binary, 3 for ternary, etc.
//'
// [[Rcpp::export]]
IntegerVector transitioncount(IntegerVector x, int numbersystem){
  int newnumsys = pow(numbersystem, 2);
  IntegerVector guide = seq(0, newnumsys - 1);
  guide.attr("dim") = IntegerVector::create(numbersystem, numbersystem);
  IntegerVector out(newnumsys);
  for(int i = 1; i < x.size(); i++) out[guide(x[i], x[i - 1])]++;
  out.attr("dim") = IntegerVector::create(numbersystem, numbersystem);
  return(out);
}


// [[Rcpp::export]]
IntegerVector emissioncount(IntegerVector states, int statenumbersystem,
                            IntegerVector residues, int resnumbersystem){
  int newnumsys = statenumbersystem * resnumbersystem;
  IntegerVector guide = seq(0, newnumsys - 1);
  guide.attr("dim") = IntegerVector::create(statenumbersystem, resnumbersystem);
  IntegerVector out(newnumsys);
  for(int i = 0; i < states.size(); i++) out[guide(states[i], residues[i])]++;
  out.attr("dim") = IntegerVector::create(statenumbersystem, resnumbersystem);
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


// [[Rcpp::export]]
IntegerVector DNA2pentadecimal(RawVector x){
  // a is a vector of raw bytes in format of Paradis (2007)
  // return order A, T, G, C, S, W, R, Y, K, M, B, V, H, D, N (same as NUC4.4)
  // many times faster than following R function:
  // fun <- function(v, residues){
  //   vcoded <- integer(length(v))
  //   for(i in seq_along(residues)) vcoded[v == residues[i]] <- i
  //     vcoded - 1
  // }
  IntegerVector res(x.size());
  RawVector tmp = RawVector::create(4, 55, 0, 136, 72, 199, 40, 24,
                                    224, 176, 208, 112, 240, 8, 160,
                                    144, 96, 80, 2);
  for(int i = 0; i < x.size(); i++){
    if((x[i] & tmp[1]) == tmp[2]){ // is purine?
      if(x[i] == tmp[3]){
        res[i] = 0; // A
      } else if(x[i] == tmp[4]){
        res[i] = 2; // G
      } else{
        res[i] = 6; // R (A or G)
      }
    } else if((x[i] & tmp[5]) == tmp[2]){ // is pyrimidine?
      if(x[i] == tmp[6]){
        res[i] = 3; // C
      } else if(x[i] == tmp[7]){
        res[i] = 1; // T
      } else{
        res[i] = 7; // Y (C or T)
      }
    }else if(x[i] == tmp[14]){
      res[i] = 9; // M (A or C)
    }else if(x[i] == tmp[15]){
      res[i] = 5; // W (A or T)
    }else if(x[i] == tmp[16]){
      res[i] = 4; // S (G or C)
    }else if(x[i] == tmp[17]){
      res[i] = 8; // K (G or T)
    }else if(x[i] == tmp[8]){
      res[i] = 11;
      // V (A or C or G)
    } else if(x[i] == tmp[9]){
      res[i] = 12;
      // H (A or C or T)
    } else if(x[i] == tmp[10]){
      res[i] = 13;
      // D (A or G or T)
    } else if(x[i] == tmp[11]){
      res[i] = 10;
      // B (C or G or T)
    } else if(x[i] == tmp[12]){
      res[i] = 14;   // N
    } else if(x[i] == tmp[0] | x[i] == tmp[18]){
      res[i] = NA_INTEGER;
    }else{
      throw Rcpp::exception("Invalid byte");
    }
  }
  return(res);
}


// // [[Rcpp::export]]
// IntegerVector DNA2pentadecimal(RawVector x){
//   // a is a raw byte in format of Paradis (2007)
//   IntegerVector res(1);
//   RawVector tmp = RawVector::create(4, 55, 0, 136, 72, 199, 40, 24,
//                                     224, 176, 208, 112, 240, 8, 160,
//                                     144, 96, 80);
//   if(x[0] <= tmp[0]){
//     throw Rcpp::exception("Input sequence contains gaps or unknown characters");
//   }else{
//     //NUC4.4 order
//     if((x[0] & tmp[1]) == tmp[2]){ // is purine?
//       if(x[0] == tmp[3]){
//         res[0] = 0; // A
//       } else if(x[0] == tmp[4]){
//         res[0] = 2; // G
//       } else{
//         res[0] = 6; // R (A or G)
//       }
//     } else if((x[0] & tmp[5]) == tmp[2]){ // is pyrimidine?
//       if(x[0] == tmp[6]){
//         res[0] = 3; // C
//       } else if(x[0] == tmp[7]){
//         res[0] = 1; // T
//       } else{
//         res[0] = 7; // Y (C or T)
//       }
//     }else if(x[0] == tmp[14]){
//       res[0] = 9; // M (A or C)
//     }else if(x[0] == tmp[15]){
//       res[0] = 5; // W (A or T)
//     }else if(x[0] == tmp[16]){
//       res[0] = 4; // S (G or C)
//     }else if(x[0] == tmp[17]){
//       res[0] = 8; // K (G or T)
//     }else if(x[0] == tmp[8]){
//       res[0] = 11;
//       // V (A or C or G)
//     } else if(x[0] == tmp[9]){
//       res[0] = 12;
//       // H (A or C or T)
//     } else if(x[0] == tmp[10]){
//       res[0] = 13;
//       // D (A or G or T)
//     } else if(x[0] == tmp[11]){
//       res[0] = 10;
//       // B (C or G or T)
//     } else if(x[0] == tmp[12]){
//       res[0] = 14;   // N
//     } else {
//       throw Rcpp::exception("Invalid byte");
//     }
//   }
//   return(res);
// }



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





