#include <Rcpp.h>
using namespace Rcpp;
//' Optimal path of sequence through model.
//'
//' \code{ViterbiC} finds the optimal path of a sequence through a HMM
//' or PHMM and returns its log-odds score.
//'
//' @param x an object of class \code{HMM} or \code{PHMM}, or a character vector.
//' @param y a character vector consisting of residues emitted by the
//' HMM or PHMM.
//' @param logspace logical argument indicating whether the emission
//' and transmission probabilities for the model(s) are logged.
//' @name ViterbiC
//' @export
// [[Rcpp::export(name = "ViterbiC.default")]]
List Viterbidefault(IntegerVector x, IntegerVector y,
                    int type, double d, double e, NumericMatrix S,
                    LogicalMatrix itertab, double offset){

  int n = x.size() + 1;
  int m = y.size() + 1;
  double sij;
  NumericVector IXcdt(3);
  IXcdt[2] = -INFINITY;
  NumericVector MMcdt(4);
  if(type == 2) MMcdt[3] = 0; else MMcdt[3] = -INFINITY;
  NumericVector IYcdt(3);
  IYcdt[0] = -INFINITY;
  int IXmax;
  int MMmax;
  int IYmax;
  NumericVector neginf = NumericVector::create(-INFINITY);
  NumericVector tmp = rep_len(neginf, n * m);
  tmp.attr("dim") = IntegerVector::create(n, m);
  NumericMatrix MIX = as<NumericMatrix>(tmp);// x aligns to gap in y
  NumericMatrix MMM = clone(MIX);// match-match
  NumericMatrix MIY = clone(MIX); // y aligns to gap in x
  IntegerMatrix PIX(n, m); // pointer matrices
  IntegerMatrix PMM = clone(PIX);
  IntegerMatrix PIY = clone(PIX);

  // initialize scoring matrices
  MMM(0, 0) = 0;
  if(type == 0){
    MIX(1, 0) = -d;
    MIY(0, 1) = -d;
    for(int i = 2; i < n; i++) MIX(i, 0) = MIX(i - 1, 0) - e;
    for(int j = 2; j < m; j++) MIY(0, j) = MIY(0, j - 1) - e;
  }else{
    for(int i = 1; i < n; i++) MIX(i, 0) = 0;
    for(int j = 1; j < m; j++) MIY(0, j) = 0;
  }
  for(int i = 1; i < n; i++) PIX(i, 0) = 0; // 0 represente IX
  for(int j = 1; j < m; j++) PIY(0, j) = 2; // 2 represents IY

  // recursion
  for(int i = 1; i < n; i++){
    for(int j = 1; j < m; j++){
      if(itertab(i, j)){
        sij = S(x[i - 1], y[j - 1]) + offset;
        IXcdt[0] = MIX(i - 1, j) - e;
        IXcdt[1] = MMM(i - 1, j) - (d + e);
        // IXcdt[2] = -INFINITY;
        MMcdt[0] = MIX(i - 1, j - 1) + sij;
        MMcdt[1] = MMM(i - 1, j - 1) + sij;
        MMcdt[2] = MIY(i - 1, j - 1) + sij;
        // IYcdt[0] = -INFINITY;
        IYcdt[1] = MMM(i, j - 1) - (d + e);
        IYcdt[2] = MIY(i, j - 1) - e;
        IXmax = which_max(IXcdt);
        MMmax = which_max(MMcdt);
        IYmax = which_max(IYcdt);
        MIX(i, j) = IXcdt[IXmax];
        MMM(i, j) = MMcdt[MMmax];
        MIY(i, j) = IYcdt[IYmax];
        PIX(i, j) = IXmax;
        PMM(i, j) = MMmax;
        PIY(i, j) = IYmax;
      }
    }
  }
  // traceback
  IntegerVector path(m + n); // if being picky should be m + n - 2 but doesnt matter
  LogicalVector keeppath(m + n);
  int counter = m + n - 1;
  int tbr;
  int tbc;
  int tbm;
  double score;
  IntegerVector startposition(2);
  if(type == 0){
    tbr = n - 1;
    tbc = m - 1;
    NumericVector brc = NumericVector::create(MIX(tbr, tbc), MMM(tbr, tbc), MIY(tbr, tbc));
    tbm = which_max(brc);
    score = brc[tbm];
    while(tbr > 0 | tbc > 0){
      keeppath[counter] = true;
      path[counter] = tbm;
      if(tbm == 0) {
        tbm = PIX(tbr, tbc);
        tbr--;
      }else if(tbm == 1){
        tbm = PMM(tbr, tbc);
        tbr--;
        tbc--;
      }else if(tbm == 2){
        tbm = PIY(tbr, tbc);
        tbc--;
      }else throw Rcpp::exception("error 1");
      counter--;
      checkUserInterrupt();
    }
  }else if(type == 1){
    tbr = n - 1;
    tbc = m - 1;
    tbm = 1;
    //tbm = 1; // MMM matrix
    score = MMM(tbr, tbc);
    for(int i = n - 2; i >= 0; i--){
      if(MMM(i, m - 1) > score){
        tbr = i;
        score = MMM(i, m - 1);
      }
    }
    for(int j = m - 1; j >= 0; j--){
      if(MMM(n - 1, j) > score){
        tbc = j;
        tbr = n - 1;
        score = MMM(n - 1, j);
      }
    }//could have equal option too using R::runif
    if(tbr < n - 1){
      for(int i = tbr + 1; i < n; i++){
        path[counter] = 0;
        keeppath[counter] = true;
        counter--;
      }
    }else if(tbc < m - 1){
      for(int j = tbc + 1; j < m; j++){
        path[counter] = 2;
        keeppath[counter] = true;
        counter--;
      }
    }
    while(tbr > 0 & tbc > 0){
      keeppath[counter] = true;
      path[counter] = tbm;
      if(tbm == 0) {
        tbm = PIX(tbr, tbc);
        tbr--;
      }else if(tbm == 1){
        tbm = PMM(tbr, tbc);
        tbr--;
        tbc--;
      }else if(tbm == 2){
        tbm = PIY(tbr, tbc);
        tbc--;
      }else throw Rcpp::exception("error 1");
      counter--;
      checkUserInterrupt();
    }
    if(tbr > 0){
      for(int i = 0; i < tbr; i++){
        path[counter] = 0;
        keeppath[counter] = true;
        counter--;
      }
    }else if(tbc > 0){
      for(int j = 0; j < tbc; j++){
        path[counter] = 2;
        keeppath[counter] = true;
        counter--;
      }
    }
  }else{
    PMM(0, 0) = 3;
    int maxind = which_max(MMM);
    score = MMM[maxind];
    tbr = maxind % (n); // remainder
    tbc = maxind/(n); // quotient
    tbm = 1;
    bool advance = PMM(tbr, tbc) != 3;
    while(advance){
      keeppath[counter] = true;
      path[counter] = tbm;
      if(tbm == 0) {
        tbm = PIX(tbr, tbc);
        tbr--;
      }else if(tbm == 1){
        tbm = PMM(tbr, tbc);
        tbr--;
        tbc--;
      }else if(tbm == 2){
        tbm = PIY(tbr, tbc);
        tbc--;
      }else throw Rcpp::exception("error 1");
      counter--;
      startposition[0] = tbr;
      startposition[1] = tbc;
      if(tbm == 0) {
        advance = PIX(tbr, tbc) != 3;
      }else if(tbm == 1) {
        advance = PMM(tbr, tbc) != 3;
      }else if(tbm == 2){
        advance = PIY(tbr, tbc) != 3;
      } else throw Rcpp::exception("error 2");
      checkUserInterrupt();
    }
  }
  IntegerVector truncpath = path[keeppath];
  //Rcout << any(duplicated(IntegerVector::create(1,2,2,3))).is_true();
  //Rcout << R::runif(-0.0001, 0.0001);
  startposition = startposition + 1; // R indexing style
  List res = List::create(Named("score") = score,
                          Named("path") = truncpath,
                          Named("start") = startposition);
  res.attr("class") = "Viterbi";
  return(res);
}


//' @name ViterbiC
//' @export
// [[Rcpp::export(name = "ViterbiC.HMM")]]
List ViterbiHMM(List x, CharacterVector y, bool logspace = false){
  //NumericVector s = clone(VECTOR_ELT(x, 0));
  NumericMatrix A = clone(VECTOR_ELT(x, 0));
  NumericMatrix E = clone(VECTOR_ELT(x, 1));
  List names = E.attr("dimnames");
  CharacterVector states = VECTOR_ELT(names, 0);
  CharacterVector residues = VECTOR_ELT(names, 1);
  int nstates = E.nrow();
  int nrolls = y.size();
  int nresidues = residues.size();
  if(!logspace){
    for(int i = 0; i < nstates; i++){
      A(i, _) = log(A(i, _));
      E(i, _) = log(E(i, _));
    }
    A(nstates, _) = log(A(nstates, _));
  }
  IntegerVector yind(nrolls); // to code residues
  //IntegerVector rseq = seq_along(residues); // this starts at 1
  //IntegerVector rseq = Range(0, nresidues - 1);
  for(int i = 0; i < nrolls; i++){
    yind[i] = match(CharacterVector::create(y[i]), residues)[0] - 1; // this is really slow
  }
  NumericMatrix V(nstates, nrolls);
  IntegerVector prerolls = Range(1, nrolls);
  CharacterVector rolls = as<CharacterVector>(prerolls);
  List vnames = List::create(Named("state") = states, Named("roll") = rolls);
  V.attr("dimnames") = vnames;
  NumericVector s = A(0, _);
  NumericVector e = A(_, 0);
  V(_, 0) = E(_, yind[0]) + s[seq(1, nstates)];
  IntegerMatrix rawP(nstates, nrolls);
  //CharacterMatrix P(nstates, nrolls);
  //P.attr("dimnames") = vnames;
  NumericMatrix tmp(nstates, nstates);
  NumericVector colmaxs(nstates);
  IntegerVector maxstates(nstates);
  for(int i = 1; i < nrolls; i++){
    for(int j = 0; j < nstates; j++){
      for(int k = 0; k < nstates; k++){
        tmp(j, k) = V(j, i - 1) + A(j + 1, k + 1);
      }
    }
    for(int l = 0; l < nstates; l++){
      int maxstate = which_max(tmp(_, l));
      colmaxs[l] = tmp(_, l)[maxstate];
      maxstates[l] = maxstate;
    }
    V(_, i) = E(_, yind[i]) + colmaxs;
    rawP(_, i) = maxstates;
    //P(_, i) = as<CharacterVector>(states[maxstates]);
    checkUserInterrupt(); //could speed up slightly by checking less frequently
  }
  NumericVector ak0(nstates);
  bool allinfinite = all(e[seq(1, nstates)] == rep(-INFINITY, nstates)).is_true();

  if(!allinfinite){ak0 += e[seq(1, nstates)];}

  int maxstate = which_max(V(_, nrolls - 1) + ak0);
  double score = V(_, nrolls - 1)[maxstate] + ak0[maxstate];
  IntegerVector rawpath(nrolls);
  CharacterVector path(nrolls);
  rawpath[nrolls - 1] = maxstate;
  path[nrolls - 1] = states[maxstate];
  for(int i = nrolls - 2; i >= 0; i--){
    rawpath[i] = rawP(maxstate, i + 1);
    path[i] = states[rawpath[i]];
    maxstate = rawpath[i];
  }
  List out = List::create(//Named("V") = V,
    //Named("P") = P,
    Named("score") = score,
    Named("path") = path);
  out.attr("class") = "Viterbi";
  return out;
}

//' @rdname ViterbiC
//' @export
// [[Rcpp::export(name = "ViterbiC.PHMM")]]
List ViterbiPHMM(List x, CharacterVector y, bool logspace = false){
  Rcout << "sorry Viterbi method is not available for PHMMs yet\n";
  return 0;
}
