#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double DNAprobC2(int x, NumericVector probs){
  // a is an integer vector with length 1 and arity = 15, representing DNA with ambiguities
  // format of Paradis (2007). Eg output of a = DNA2pentadecimal(x.DNAbin)[1]
  // probs is a 4-element numeric vector of probabilities for the set {a,c,g,t}
  if(probs.size() != 4){
    throw Rcpp::exception("probs argument must be a numeric vector of length 4");
  }
  if(x < 4){
    return(probs[x]); //knownbase
  }else if(x < 10){
    if(x < 7){
      if(x == 4){
        return((probs[2] + probs[3])/2);
      }else if(x == 5){
        return((probs[0] + probs[1])/2);
      }else{
        return((probs[0] + probs[2])/2);
      }
    }else{
      if(x == 7){
        return((probs[1] + probs[3])/2);
      }else if(x == 8){
        return((probs[1] + probs[2])/2);
      }else{
        return((probs[0] + probs[3])/2);
      }
    }
  }else if(x < 14){
    if(x == 10){
      return((probs[1] + probs[2] + probs[3])/3); //B
    }else if(x == 11){
      return((probs[0] + probs[2] + probs[3])/3); //V
    }else if(x == 12){
      return((probs[0] + probs[1] + probs[3])/3); //H
    }else{
      return((probs[0] + probs[1] + probs[2])/3); //D
    }
  }else if(x == 14){
    return((probs[0] + probs[1] + probs[2] + probs[3])/4);
  }else throw Rcpp::exception("expected integers between 0 and 14");
  return(0);
}



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
// [[Rcpp::export]]
List Viterbi_default(IntegerVector x, IntegerVector y,
                    int type, double d, double e, NumericMatrix S,
                    IntegerVector windowspace, double offset){

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
      if(j - i >= windowspace[0] & j - i <= windowspace[1]){
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
    checkUserInterrupt();
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
// [[Rcpp::export]]
List Viterbi_HMM(IntegerVector y, NumericMatrix A, NumericMatrix E){
  List names = E.attr("dimnames");
  CharacterVector states = VECTOR_ELT(names, 0);
  CharacterVector residues = VECTOR_ELT(names, 1);
  int nstates = E.nrow();
  int nrolls = y.size();
  int nresidues = residues.size();
  NumericMatrix V(nstates, nrolls);
  IntegerVector prerolls = Range(1, nrolls);
  CharacterVector rolls = as<CharacterVector>(prerolls);
  List vnames = List::create(Named("state") = states, Named("roll") = rolls);
  V.attr("dimnames") = vnames;
  NumericVector s = A(0, _);
  NumericVector e = A(_, 0);
  V(_, 0) = E(_, y[0]) + s[seq(1, nstates)];
  IntegerMatrix P(nstates, nrolls);
  //CharacterMatrix P(nstates, nrolls);
  P.attr("dimnames") = vnames;
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
    V(_, i) = E(_, y[i]) + colmaxs;
    P(_, i) = maxstates;
    //P(_, i) = as<CharacterVector>(states[maxstates]);
    checkUserInterrupt(); //could speed up slightly by checking less frequently
  }
  NumericVector ak0(nstates);
  bool allinfinite = all(e[seq(1, nstates)] == rep(-INFINITY, nstates)).is_true();
  if(!allinfinite){ak0 += e[seq(1, nstates)];}
  int maxstate = which_max(V(_, nrolls - 1) + ak0);
  double score = V(_, nrolls - 1)[maxstate] + ak0[maxstate];
  IntegerVector path(nrolls);
  //CharacterVector charpath(nrolls);
  path[nrolls - 1] = maxstate;
  //charpath[nrolls - 1] = states[maxstate];
  for(int i = nrolls - 2; i >= 0; i--){
    path[i] = P(maxstate, i + 1);
    //charpath[i] = states[path[i]];
    maxstate = path[i];
  }
  List out = List::create(
    Named("score") = score,
    Named("path") = path,
    Named("V") = V,
    Named("P") = P);
  out.attr("class") = "Viterbi";
  return out;
}

//' @rdname ViterbiC
//' @export
// [[Rcpp::export]]
List Viterbi_PHMM(IntegerVector y, NumericMatrix A, NumericMatrix E, NumericVector qe,
                  NumericVector qey, int type, IntegerVector windowspace, double offset,
                  bool DI, bool ID, bool DNA){
  int n = E.ncol() + 1;
  int m = y.size() + 1;
  double sij;
  NumericVector Dcdt(3);
  Dcdt[2] = -INFINITY;
  NumericVector Mcdt(4);
  if(type == 2) Mcdt[3] = 0; else Mcdt[3] = -INFINITY;
  NumericVector Icdt(3);
  Icdt[0] = -INFINITY;
  int Dmax;
  int Mmax;
  int Imax;
  NumericVector neginf = NumericVector::create(-INFINITY);
  NumericVector tmp = rep_len(neginf, n * m);
  tmp.attr("dim") = IntegerVector::create(n, m);
  NumericMatrix Dmatrix = as<NumericMatrix>(tmp);// x aligns to gap in y
  NumericMatrix Mmatrix = clone(Dmatrix);// match-match
  NumericMatrix Imatrix = clone(Dmatrix); // y aligns to gap in x
  IntegerMatrix Dpointer(n, m); // pointer matrices
  IntegerMatrix Mpointer = clone(Dpointer);
  IntegerMatrix Ipointer = clone(Dpointer);

  // initialize scoring matrices
  Mmatrix(0, 0) = 0;
  if(type == 0){
    Dmatrix(1, 0) = A(3, 0); // Match -> Delete transition prob for position 0
    Imatrix(0, 1) = A(5, 0) + qey[0];
    // V[-1, 1, 1] <- cumsum(c(0, A["DD", 2:(n - 1)])) + A["MD", 1]
    for(int i = 2; i < n; i++) Dmatrix(i, 0) = Dmatrix(i - 1, 0) + A(0, i - 1);
    //V[1, 2, "I"] <- A["MI", 1] + qey[1]
    //for(j in 3:m) V[1, j, "I"] <- V[1, j - 1, "I"] + A["II", 1] + qey[j - 1]
    for(int j = 2; j < m; j++) Imatrix(0, j) = Imatrix(0, j - 1) + A(8, 0) + qey[j - 1];
  }else{
    for(int i = 1; i < n; i++) Dmatrix(i, 0) = 0;
    for(int j = 1; j < m; j++) Imatrix(0, j) = 0; // qey??
  }
  Dpointer(1, 0) = 1;
  Ipointer(0, 1) = 1;
  for(int i = 2; i < n; i++) Dpointer(i, 0) = 0; // 0 represents IX
  for(int j = 2; j < m; j++) Ipointer(0, j) = 2; // 2 represents IY

  // recursion
  for(int i = 1; i < n; i++){
    for(int j = 1; j < m; j++){
      if(j - i >= windowspace[0] & j - i <= windowspace[1]){
        if(DNA){
          sij = DNAprobC2(y[j - 1], E(_, i - 1)) + offset;
        }else{
          sij = E(y[j - 1], i - 1) + offset;
        }
        Dcdt[0] = Dmatrix(i - 1, j) + A(0, i - 1); //DD
        Dcdt[1] = Mmatrix(i - 1, j) + A(3, i - 1); //MD
        if(ID) Dcdt[2] = Imatrix(i - 1, j) + A(6, i - 1);  //ID
        Mcdt[0] = Dmatrix(i - 1, j - 1) + A(1, i - 1) + sij; //DM
        Mcdt[1] = Mmatrix(i - 1, j - 1) + A(4, i - 1) + sij; //MM
        Mcdt[2] = Imatrix(i - 1, j - 1) + A(7, i - 1) + sij; //IM
        // Mcdt[3] is either 0 or -inf, depending on type
        if(DI) Icdt[0] = Dmatrix(i, j - 1) + A(2, i); //DI
        Icdt[1] = Mmatrix(i, j - 1) + A(5, i); //MI
        Icdt[2] = Imatrix(i, j - 1) +A (8, i); //II
        Dmax = which_max(Dcdt);
        Mmax = which_max(Mcdt);
        Imax = which_max(Icdt);
        Dmatrix(i, j) = Dcdt[Dmax];
        Mmatrix(i, j) = Mcdt[Mmax];
        Imatrix(i, j) = Icdt[Imax] + qey[j - 1];
        Dpointer(i, j) = Dmax;
        Mpointer(i, j) = Mmax;
        Ipointer(i, j) = Imax;
      }
    }
    checkUserInterrupt();
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
    // LLcdt <- c(V[n, m, "D"] + A["DM", n],
    //            V[n, m, "M"] + A["MM", n],
    //                            V[n, m, "I"] + A["IM", n])
    NumericVector LLcdt = NumericVector::create(Dmatrix(n - 1, m - 1) + A(1, n - 1), //DM
                                                Mmatrix(n - 1, m - 1) + A(4, n - 1), //MM
                                                Dmatrix(n - 1, m - 1) + A(7, n - 1)); //IM

    tbm = which_max(LLcdt);
    score = LLcdt[tbm];
    tbr = n - 1; // traceback row
    tbc = m - 1; // traceback column
    // brc = bottom right corner
    //NumericVector brc = NumericVector::create(Dmatrix(tbr, tbc), Mmatrix(tbr, tbc), Imatrix(tbr, tbc));
    //tbm = which_max(brc);
    //score = brc[tbm];
    while(tbr > 0 | tbc > 0){
      keeppath[counter] = true;
      path[counter] = tbm;
      if(tbm == 0) {
        tbm = Dpointer(tbr, tbc);
        tbr--;
      }else if(tbm == 1){
        tbm = Mpointer(tbr, tbc);
        tbr--;
        tbc--;
      }else if(tbm == 2){
        tbm = Ipointer(tbr, tbc);
        tbc--;
      }else throw Rcpp::exception("error 1");
      counter--;
      checkUserInterrupt();
    }
  }else if(type == 1){
    tbr = n - 1;
    tbc = m - 1;
    tbm = 1;
    //tbm = 1; // Mmatrix
    score = Mmatrix(tbr, tbc);
    for(int i = n - 2; i >= 0; i--){
      if(Mmatrix(i, m - 1) > score){
        tbr = i;
        score = Mmatrix(i, m - 1);
      }
    }
    for(int j = m - 1; j >= 0; j--){
      if(Mmatrix(n - 1, j) > score){
        tbc = j;
        tbr = n - 1;
        score = Mmatrix(n - 1, j);
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
        tbm = Dpointer(tbr, tbc);
        tbr--;
      }else if(tbm == 1){
        tbm = Mpointer(tbr, tbc);
        tbr--;
        tbc--;
      }else if(tbm == 2){
        tbm = Ipointer(tbr, tbc);
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
    Mpointer(0, 0) = 3;
    int maxind = which_max(Mmatrix);
    score = Mmatrix[maxind];
    tbr = maxind % (n); // remainder
    tbc = maxind/(n); // quotient
    tbm = 1;
    bool advance = Mpointer(tbr, tbc) != 3;
    while(advance){
      keeppath[counter] = true;
      path[counter] = tbm;
      if(tbm == 0) {
        tbm = Dpointer(tbr, tbc);
        tbr--;
      }else if(tbm == 1){
        tbm = Mpointer(tbr, tbc);
        tbr--;
        tbc--;
      }else if(tbm == 2){
        tbm = Ipointer(tbr, tbc);
        tbc--;
      }else throw Rcpp::exception("error 1");
      counter--;
      startposition[0] = tbr;
      startposition[1] = tbc;
      if(tbm == 0) {
        advance = Dpointer(tbr, tbc) != 3;
      }else if(tbm == 1) {
        advance = Mpointer(tbr, tbc) != 3;
      }else if(tbm == 2){
        advance = Ipointer(tbr, tbc) != 3;
      } else throw Rcpp::exception("error 2");
      checkUserInterrupt();
    }
  }
  IntegerVector truncpath = path[keeppath];
  //Rcout << any(duplicated(IntegerVector::create(1,2,2,3))).is_true();
  //Rcout << R::runif(-0.0001, 0.0001);
  startposition = startposition + 1; // R's indexing style, starts at 1 instead of 0
  //Dmatrix.attr("dim") = NULL;
  //NumericVector V = NumericVector::create(as<NumericVector>(Dmatrix));
  // NumericVector V = NumericVector::create(as<NumericVector>(Dmatrix),
  //                                         as<NumericVector>(Mmatrix),
  //                                         as<NumericVector>(Imatrix));
  List res = List::create(Named("score") = score,
                          Named("path") = truncpath,
                          Named("start") = startposition,
                          Named("Dmatrix") = Dmatrix,
                          Named("Mmatrix") = Mmatrix,
                          Named("Imatrix") = Imatrix,
                          Named("Dpointer") = Dpointer,
                          Named("Mpointer") = Mpointer,
                          Named("Ipointer") = Ipointer);
  //List res = List::create(Named("Mmat") = Dpointer);
  res.attr("class") = "Viterbi";
  return(res);
}
