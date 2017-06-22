#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export(name = ".map")]]
LogicalVector map(NumericMatrix ecs, LogicalMatrix notgaps, List pseudocounts,
                   NumericVector seqweights,
                   NumericVector qe, double lambda = 0){
  NumericVector Apscs = clone(VECTOR_ELT(pseudocounts, 0));
  NumericVector Epscs = clone(VECTOR_ELT(pseudocounts, 1));
  NumericMatrix ecs2 = clone(ecs);
  int L = ecs.ncol();
  int nres = ecs.nrow();
  int n = notgaps.nrow();
  NumericVector S(L + 2);
  IntegerVector sigma(L + 2);
  LogicalVector Scutter(L + 2);
  LogicalVector res(L + 2);
  res[L + 1] = true;
  IntegerVector columnnumbers = Range(0, L + 1);
  CharacterVector resnames = as<CharacterVector>(columnnumbers);
  res.attr("names") = resnames;
  LogicalVector rescutter(L + 2);
  rescutter[0] = true;
  rescutter[L + 1] = true;
  rescutter = !rescutter;
  NumericVector tmp(nres);
  NumericVector M(L + 2);
  for(int j = 0; j < L; j++){
    tmp = ecs(_, j) + Epscs;
    M[j + 1] = sum(ecs(_, j) * log(tmp/sum(tmp))); // could poss speed this up a bit
  }
  NumericVector ecsj(nres);
  NumericVector tcsij(9);
  NumericVector icsj(n);
  LogicalVector zeroinserts(n);
  LogicalMatrix tcounts(n, 9);
  NumericVector cxy(9);
  NumericVector axy(9);
  NumericVector denoms(3);
  IntegerVector denomindex = IntegerVector::create(0,0,0,1,1,1,2,2,2);
  for(int j = 1; j < L + 2; j++){
    NumericVector tau(j);
    NumericVector iota(j);
    NumericVector ecsij = clone(ecsj);
    NumericVector icsij = clone(icsj);
    for(int i = 0; i < j; i++){
      if(i < j - 1){
        iota[i] = sum(ecsij * qe);
        ecsij = ecsij - ecs(_, i);
        zeroinserts = icsij < 0.00001;
        tcounts(_, 0) = !notgaps(_, i) & !notgaps(_, j) & zeroinserts;
        tcounts(_, 1) = !notgaps(_, i) & notgaps(_, j) & zeroinserts;
        tcounts(_, 2) = !notgaps(_, i) & !zeroinserts;
        tcounts(_, 3) = notgaps(_, i) & !notgaps(_, j) & zeroinserts;
        tcounts(_, 4) = notgaps(_, i) & notgaps(_, j) & zeroinserts;
        tcounts(_, 5) = notgaps(_, i) & !zeroinserts;
        tcounts(_, 6) = !notgaps(_, j) & !zeroinserts;
        tcounts(_, 7) = notgaps(_, j) & !zeroinserts;
        //tcounts column 8 not required as can be calculated separately
        for(int k = 0; k < 8; k++){
          NumericVector tmp = seqweights[tcounts(_, k)];
          tcsij[k] = sum(tmp);
        }
        tcsij[8] = sum(icsij) - (tcsij[2] + tcsij[5]); //II
        NumericVector tmp(n);
        tmp[notgaps(_, i + 1)] = 1;
        icsij = icsij - seqweights * tmp;
        // ends up as vec of zeros at end of each i cycle
      }else{
        tcounts(_, 0) = !notgaps(_, i) & !notgaps(_, j); //DD
        tcounts(_, 1) = !notgaps(_, i) & notgaps(_, j); //DM
        tcounts(_, 3) = notgaps(_, i) & !notgaps(_, j); //MD
        tcounts(_, 4) = notgaps(_, i) & notgaps(_, j); //MM
        NumericVector DD = seqweights[tcounts(_, 0)];
        NumericVector DM = seqweights[tcounts(_, 1)];
        NumericVector MD = seqweights[tcounts(_, 3)];
        NumericVector MM = seqweights[tcounts(_, 4)];
        tcsij[0] = sum(DD);
        tcsij[1] = sum(DM);
        tcsij[2] = 0;
        tcsij[3] = sum(MD);
        tcsij[4] = sum(MM);
        tcsij[5] = 0;
        tcsij[6] = 0;
        tcsij[7] = 0;
        tcsij[8] = 0;
      }
      cxy = tcsij + Apscs;
      denoms[0] = cxy[0] + cxy[1] + cxy[2];
      denoms[1] = cxy[3] + cxy[4] + cxy[5];
      denoms[2] = cxy[6] + cxy[7] + cxy[8];
      for(int k = 0; k < 9; k ++) axy[k] = cxy[k]/denoms[denomindex[k]];
      tau[i] = sum(cxy * log(axy));
    }
    if(j < L + 1){
      ecsj = ecsj + ecs(_, j - 1);
      NumericVector tmp(n);
      tmp[notgaps(_, j)] = 1;
      icsj = icsj + seqweights * tmp;
    }
    Scutter[j - 1] = true;
    NumericVector Stmp = S[Scutter];
    NumericVector tmp = Stmp + tau + iota + M[j] + lambda;
    sigma[j] = which_max(tmp) + 1;
    S[j] = tmp[sigma[j] - 1];
    checkUserInterrupt();
  }
  int j = sigma[L + 1];
  while(j > 0){
    res[j] = true;
    j = sigma[j] - 1;
    checkUserInterrupt();
  }
  LogicalVector out = res[rescutter];
  return(out);
}

