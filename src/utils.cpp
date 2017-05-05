#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(name = ".acount")]]
IntegerVector acount(IntegerVector x, int arity){
  // x is an integer vector in a certain coding scheme (starting at 0)
  // arity is an integer representing the numbering system of the input vector,
  // 2 for binary, 3 for ternary, 4 for quarternary, etc.
  // returns the transition frequencies as an integer vector
  int newarity = pow(arity, 2);
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
      else if((x(i, j) == 1) & (x(i, k) == 0)) out(3, module)= out(3, module) + seqweights[i];
      else if((x(i, j) == 1) & (x(i, k) == 1)) out(4, module)= out(4, module) + seqweights[i];
      else if((x(i, j) == 1) & (x(i, k) == 2)) out(5, module)= out(5, module) + seqweights[i];
      else if((x(i, j) == 2) & (x(i, k) == 0)) out(6, module)= out(6, module) + seqweights[i];
      else if((x(i, j) == 2) & (x(i, k) == 1)) out(7, module)= out(7, module) + seqweights[i];
      else if((x(i, j) == 2) & (x(i, k) == 2)) out(8, module)= out(8, module) + seqweights[i];
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

// used for fragmenting a sequence into matches and inserts, given a path
// [[Rcpp::export(name = ".fragR")]]
List fragR(RawVector x, IntegerVector path, int l, RawVector gap){
  // x is a sequence vector of mode "character"
  // path is a trinary integer vector, must end with 1
  // l is PHMM model size (saves on computation time)
  // returns fragmented sequence list of length 2l + 1
  int outlen = 2 * l + 1; // length of returned list, includes
  //one element for each emitting match state and one for each insert state
  List out = List(outlen);
  int seqcounter = 0; // keeps track of current position in the sequence
  int pathcounter = 0; // keeps track of current position in the path
  int insstart = 0; // keeps track of position of beginning of current insert
  int outcounter = 0; // keeps track of current position in the returned list
  RawVector nothing(0);
  bool frominsert = false;
  while(outcounter < outlen){
    if(path[pathcounter] == 1){// match state
      if(!frominsert) out[outcounter] = nothing;
      outcounter ++; // advances to match state in 'out'
      if(outcounter == outlen) break;
      out[outcounter] = RawVector::create(x[seqcounter]);
      outcounter ++; // advances to the next insert state
      pathcounter++;
      seqcounter++;
      frominsert = false;
    }else if(path[pathcounter] == 0){ //delete state
      if(!frominsert) out[outcounter] = nothing;
      outcounter ++; // advances to match state in 'out'
      if(outcounter == outlen) break;
      out[outcounter] = gap;
      outcounter ++; // advances to the next insert state
      pathcounter++;
      // seqcounter++; // no advance for seqcounter
      frominsert = false;
    }else if(path[pathcounter] == 2){// insert state
      insstart = seqcounter;
      pathcounter++;
      seqcounter++;
      while(path[pathcounter] == 2){
        pathcounter++;
        seqcounter++;
        // no advance in outcounter
      }
      out[outcounter] = x[seq(insstart, seqcounter - 1)];
      frominsert = true;
    }else throw Rcpp::exception("Invalid integer in path");
    //checkUserInterrupt();
  }
  return(out);
}



// [[Rcpp::export(name = ".fragC")]]
List fragC(CharacterVector x, IntegerVector path, int l, CharacterVector gap){
  // x is a sequence vector of mode "character"
  // path is a trinary integer vector, must end with 1
  // l is PHMM model size (saves on computation time)
  // returns fragmented sequence list of length 2l + 1
  int outlen = 2 * l + 1; // length of returned list, includes
  //one element for each emitting match state and one for each insert state
  List out = List(outlen);
  int seqcounter = 0; // keeps track of current position in the sequence
  int pathcounter = 0; // keeps track of current position in the path
  int insstart = 0; // keeps track of position of beginning of current insert
  int outcounter = 0; // keeps track of current position in the returned list
  // CharacterVector gap = CharacterVector::create(gap);
  CharacterVector nothing(0);
  bool frominsert = false;
  while(outcounter < outlen){
    if(path[pathcounter] == 1){// match state
      if(!frominsert) out[outcounter] = nothing;
      outcounter ++; // advances to match state in 'out'
      if(outcounter == outlen) break;
      out[outcounter] = CharacterVector::create(x[seqcounter]);
      outcounter ++; // advances to the next insert state
      pathcounter++;
      seqcounter++;
      frominsert = false;
    }else if(path[pathcounter] == 0){ //delete state
      if(!frominsert) out[outcounter] = nothing;
      outcounter ++; // advances to match state in 'out'
      if(outcounter == outlen) break;
      out[outcounter] = gap;
      outcounter ++; // advances to the next insert state
      pathcounter++;
      // seqcounter++; // no advance for seqcounter
      frominsert = false;
    }else if(path[pathcounter] == 2){// insert state
      insstart = seqcounter;
      pathcounter++;
      seqcounter++;
      while(path[pathcounter] == 2){
        pathcounter++;
        seqcounter++;
        // no advance in outcounter
      }
      out[outcounter] = x[seq(insstart, seqcounter - 1)];
      frominsert = true;
    }else throw Rcpp::exception("Invalid integer in path");
    // checkUserInterrupt();
  }
  return(out);
}


