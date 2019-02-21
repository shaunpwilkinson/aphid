#include <Rcpp.h>
using namespace Rcpp;

//' Sum of logged probabilities.
//'
//' \code{"logsum"} takes a vector of logged probabilities (neagtive values)
//'   and returns its sum.
//'
//' @param x a numeric vector of logged probabilities.
//' @return returns a single numeric value representing the logged sum of the
//'   values in the input vector.
//' @details This is a simple compiled function that exponentiates the values
//' in the input vector, finds their sum, and returns the log of that value.
//' @author Shaun Wilkinson
//'
// [[Rcpp::export]]
double logsum(NumericVector x){
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


// Find the position of the maximum value.
//
// Returns the position (index) of the maximum value in a numeric or integer vector,
//   with ties broken at random.
//
// x a numeric or integer vector.
// start integer indicating the index of the first value
//   of the vector. Defaults to 1 (consistent the indexing system used by R,
//   as opposed to that used by C and C++ in which the first value is given
//   an index of 0).
// returns an integer representing the position (index) of the maximum
//   value in the input vector.
// This is a simple function similar to which.max (base R) but
//   with random breaking of ties. Unlike which.max, logical vectors are
//   not accepted.
//
// [[Rcpp::export(name = ".whichmax")]]
int whichmax(NumericVector x, int start = 1){
  int maxind = 0;
  bool ties = false;
  IntegerVector backups(x.size());
  int counter = 0;
  for(int i = 1; i < x.size(); i++){
    if(x[i] > x[maxind]){
      maxind = i;
      if(ties){
        for(int j = 0; j < counter; j++) backups[j] = 0;
        counter = 0;
        ties = false;
      }
    }else if(x[i] == x[maxind]){
      backups[counter] = i;
      counter++;
      ties = true;
    }
    //checkUserInterrupt();
  }
  if(ties){
    backups[counter] = maxind;
    double rando = R::runif(0, 1);
    double counternum = counter;
    double increment = 1/(counternum + 1);
    double ceiling = increment;
    for(int i = 0; i <= counter; i++){
      if(rando < ceiling){
        return(backups[i] + start);
      }else{
        ceiling += increment;
      }
    }
  }
  return(maxind + start);
}

// [[Rcpp::export(name = ".probDNA")]]
double probDNA(int x, NumericVector probs){
  // x is an integer vector with length 1 and arity = 15, representing DNA with ambiguities
  // order A, T, G, C, S, W, R, Y, K, M, B, V, H, D, N (same as NUC.4.4)
  // format of Paradis (2007).
  // probs is a 4-element numeric vector of probabilities for the set {A,T,G,C} (in that order)
  if(probs.size() != 4){
    throw Rcpp::exception("probs argument must be a numeric vector of length 4");
  }
  if(x < 4){
    return(probs[x]); //knownbase
  }else if(x < 10){
    if(x < 7){
      if(x == 4){
        return(log((exp(probs[2]) + exp(probs[3]))/2));
      }else if(x == 5){
        return(log((exp(probs[0]) + exp(probs[1]))/2));
      }else{
        return(log((exp(probs[0]) + exp(probs[2]))/2));
      }
    }else{
      if(x == 7){
        return(log((exp(probs[1]) + exp(probs[3]))/2));
      }else if(x == 8){
        return(log((exp(probs[1]) + exp(probs[2]))/2));
      }else{
        return(log((exp(probs[0]) + exp(probs[3]))/2));
      }
    }
  }else if(x < 14){
    if(x == 10){
      return(log((exp(probs[1]) + exp(probs[2]) + exp(probs[3]))/3)); //B
    }else if(x == 11){
      return(log((exp(probs[0]) + exp(probs[2]) + exp(probs[3]))/3)); //V
    }else if(x == 12){
      return(log((exp(probs[0]) + exp(probs[1]) + exp(probs[3]))/3)); //H
    }else{
      return(log((exp(probs[0]) + exp(probs[1]) + exp(probs[2]))/3)); //D
    }
  }else if(x == 14){
    return(log(0.25));
  }else throw Rcpp::exception("expected integers between 0 and 14");
  return(0);
}

// [[Rcpp::export(name = ".probAA")]]
double probAA(int x, NumericVector probs){
  // a is an integer vector with length 1 and arity = 27,
  // representing AA with ambiguities
  // order  ACDEFGHIKLMNPQRSTVWY X BJZ OU *
  // format of Paradis (2007).
  // probs is a 20-element numeric vector of probabilities for the set
  // {ACDEFGHIKLMNPQRSTVWY} (in that order)
  if(probs.size() != 20){
    throw Rcpp::exception("probs argument must be a numeric vector of length 20");
  }
  if(x < 20){
    return(probs[x]); //knownbase
  }else if(x == 20){
    return(log(0.05));
  }else if(x < 24){
    if(x == 21){
      return(log((exp(probs[2]) + exp(probs[11]))/2)); //B (D or N)
    }else if(x == 22){
      return(log((exp(probs[7]) + exp(probs[9]))/2)); //J (I or L)
    }else{
      return(log((exp(probs[3]) + exp(probs[13]))/2)); //Z (E or Q)
    }
  }else if(x == 24){
    return(probs[8]); //O switch to K
  }else if(x == 25){
    return(probs[1]); //U switch to C
  }else if(x == 26){
    return(min(probs)); //* check this
  }else throw Rcpp::exception("expected integers between 0 and 26");
  return(0);
}

// [[Rcpp::export(name = ".ViterbiD")]]
List ViterbiD(IntegerVector x, IntegerVector y,
                    int type, double d, double e, NumericMatrix S,
                    IntegerVector windowspace, double offset = 0){
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
  if(type == 2){
    for(int i = 0; i < n; i++){
      PIX(i, 0) = 3;
      PMM(i, 0) = 3;// bit over the top but ok
      PIY(i, 0) = 3;
    }
    for(int j = 0; j < m; j++){
      PIX(0, j) = 3;
      PMM(0, j) = 3;
      PIY(0, j) = 3;
    }
  }else{
    for(int i = 1; i < n; i++) PIX(i, 0) = 0; // 0 represente IX
    for(int j = 1; j < m; j++) PIY(0, j) = 2; // 2 represents IY
  }
  // recursion
  for(int i = 1; i < n; i++){
    for(int j = 1; j < m; j++){
      if(((j - i) >= windowspace[0]) & ((j - i) <= windowspace[1])){
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
        IXmax = whichmax(IXcdt, 0);
        MMmax = whichmax(MMcdt, 0);
        IYmax = whichmax(IYcdt, 0);
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
    tbm = whichmax(brc, 0);
    score = brc[tbm];
    while((tbr > 0) | (tbc > 0)){
      if((tbr < 0) | (tbc < 0)) throw Rcpp::exception("Error 2");
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
      }else throw Rcpp::exception("Error 1");

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
    while((tbr > 0) & (tbc > 0)){
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
    int maxind = whichmax(MMM, 0);
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
      }else throw Rcpp::exception("Error 3");
      counter--;
      startposition[0] = tbr;
      startposition[1] = tbc;
      if(tbm == 0) {
        advance = PIX(tbr, tbc) != 3;
      }else if(tbm == 1) {
        advance = PMM(tbr, tbc) != 3;
      }else if(tbm == 2){
        advance = PIY(tbr, tbc) != 3;
      }else throw Rcpp::exception("Error 4");
      checkUserInterrupt();
    }
  }
  IntegerVector truncpath = path[keeppath];
  //Rcout << R::runif(-0.0001, 0.0001);
  startposition = startposition + 1; // convert to R's indexing style
  List res = List::create(Named("score") = score,
                          Named("path") = truncpath,
                          Named("start") = startposition,
                          Named("IXmatrix") = MIX,
                          Named("MMmatrix") = MMM,
                          Named("IYmatrix") = MIY,
                          Named("IXpointer") = PIX,
                          Named("MMpointer") = PMM,
                          Named("IYpointer") = PIY);

  res.attr("class") = "DPA";
  return(res);
}

// [[Rcpp::export(name = ".ViterbiH")]]
List ViterbiH(IntegerVector y, NumericMatrix A, NumericMatrix E,
                 bool DNA = false, bool AA = false){
  List names = E.attr("dimnames");
  CharacterVector states = VECTOR_ELT(names, 0);
  //CharacterVector residues = VECTOR_ELT(names, 1);
  int nstates = E.nrow();
  int nrolls = y.size();
  //int nresidues = residues.size();
  NumericMatrix V(nstates, nrolls);
  IntegerVector prerolls = Range(1, nrolls);
  CharacterVector rolls = as<CharacterVector>(prerolls);
  List vnames = List::create(Named("state") = states, Named("position") = rolls);
  V.attr("dimnames") = vnames;
  NumericVector s = A(0, _);
  NumericVector e = A(_, 0);
  LogicalVector tf(nstates + 1, true);
  //for(int k = 1; k < nstates; k++) tf[k] = true;
  tf[0] = false; // used to remove start/end row/col from A later
  IntegerMatrix P(nstates, nrolls);
  //CharacterMatrix P(nstates, nrolls);
  P.attr("dimnames") = vnames;
  NumericMatrix tmp(nstates, nstates);
  NumericVector colmaxs(nstates);
  IntegerVector maxstates(nstates);
  if(DNA){
    for(int l = 0; l < nstates; l++) V(l, 0) = probDNA(y[0], E(l,_)) + A(0, l + 1);
    for(int i = 1; i < nrolls; i++){
      for(int j = 0; j < nstates; j++){
        for(int k = 0; k < nstates; k++){
          tmp(j, k) = V(j, i - 1) + A(j + 1, k + 1);
        }
      }
      for(int l = 0; l < nstates; l++){
        P(l, i) = whichmax(tmp(_, l), 0); // 0-based indexing for cpp
        V(l, i) = probDNA(y[i], E(l,_)) + tmp(P(l, i), l);
      }
      checkUserInterrupt();
    }
  }else if(AA){
    for(int l = 0; l < nstates; l++) V(l, 0) = probAA(y[0], E(l,_)) + A(0, l + 1);
    for(int i = 1; i < nrolls; i++){
      for(int j = 0; j < nstates; j++){
        for(int k = 0; k < nstates; k++){
          tmp(j, k) = V(j, i - 1) + A(j + 1, k + 1);
        }
      }
      for(int l = 0; l < nstates; l++){
        P(l, i) = whichmax(tmp(_, l), 0); // 0-based indexing for cpp
        V(l, i) = probAA(y[i], E(l,_)) + tmp(P(l, i), l);
      }
      checkUserInterrupt();
    }
  }else{
    V(_, 0) = E(_, y[0]) + s[seq(1, nstates)];
    for(int i = 1; i < nrolls; i++){
      for(int j = 0; j < nstates; j++){
        for(int k = 0; k < nstates; k++){
          tmp(j, k) = V(j, i - 1) + A(j + 1, k + 1);
        }
      }
      for(int l = 0; l < nstates; l++){
        int maxstate = whichmax(tmp(_, l), 0);
        colmaxs[l] = tmp(_, l)[maxstate];
        maxstates[l] = maxstate;
      }
      V(_, i) = E(_, y[i]) + colmaxs;
      P(_, i) = maxstates;
      checkUserInterrupt();
    }
  }
  NumericVector ak0(nstates, 0.0);
  bool allinfinite = all(e[seq(1, nstates)] == rep(-INFINITY, nstates)).is_true();
  if(!allinfinite) ak0 += e[seq(1, nstates)];
  int maxstate = whichmax(V(_, nrolls - 1) + ak0, 0);
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
    Named("array") = V,
    Named("pointer") = P);
  out.attr("class") = "DPA";
  return out;
}

// [[Rcpp::export(name = ".ViterbiP")]]
List ViterbiP(IntegerVector y, NumericMatrix A, NumericMatrix E, NumericVector qe,
                  NumericVector qey, int type, IntegerVector windowspace, double offset = 0,
                  bool DI = false, bool ID = false, bool DNA = false, bool AA = false){
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
      if(((j - i) >= windowspace[0]) & ((j - i) <= windowspace[1])){
        if(DNA){
          sij = probDNA(y[j - 1], E(_, i - 1)) + offset;
        }else if(AA){
          sij = probAA(y[j - 1], E(_, i - 1)) + offset;
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
        Icdt[2] = Imatrix(i, j - 1) + A(8, i); //II
        Dmax = whichmax(Dcdt, 0);
        Mmax = whichmax(Mcdt, 0);
        Imax = whichmax(Icdt, 0);
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
  int counter = m + n - 1; // start at end and work backwards
  // if(type < 2){//include transition to end state if not local alignment
  //   path[counter] = 1;
  //   keeppath[counter] = true;
  //   counter--;
  // }
  int tbr;
  int tbc;
  int tbm;
  double score;
  IntegerVector startposition(2);
  if(type == 0){
    // LLcdt <- c(V[n, m, "D"] + A["DM", n],
    //            V[n, m, "M"] + A["MM", n],
    //                            V[n, m, "I"] + A["IM", n])
    NumericVector LLcdt = NumericVector::create(
      Dmatrix(n - 1, m - 1) + A(1, n - 1), //DM
      Mmatrix(n - 1, m - 1) + A(4, n - 1), //MM
      Imatrix(n - 1, m - 1) + A(7, n - 1)); //IM
    tbm = whichmax(LLcdt, 0);
    score = LLcdt[tbm];
    tbr = n - 1; // traceback row
    tbc = m - 1; // traceback column
    // brc = bottom right corner
    //NumericVector brc = NumericVector::create(Dmatrix(tbr, tbc), Mmatrix(tbr, tbc), Imatrix(tbr, tbc));
    //tbm = whichmax(brc);
    //score = brc[tbm];
    while((tbr > 0) | (tbc > 0)){
      if((tbr < 0) | (tbc < 0)) throw Rcpp::exception("Error 2");
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
      }else throw Rcpp::exception("Error 1");
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
    while((tbr > 0) & (tbc > 0)){
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
      }else throw Rcpp::exception("Error 1");
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
    int maxind = whichmax(Mmatrix, 0);
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
      }else throw Rcpp::exception("Error 5");
      counter--;
      startposition[0] = tbr;
      startposition[1] = tbc;
      if(tbm == 0) {
        advance = Dpointer(tbr, tbc) != 3;
      }else if(tbm == 1) {
        advance = Mpointer(tbr, tbc) != 3;
      }else if(tbm == 2){
        advance = Ipointer(tbr, tbc) != 3;
      } else throw Rcpp::exception("Error 6");
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
  res.attr("class") = "DPA";
  return(res);
}

// [[Rcpp::export(name = ".ViterbiPP")]]
List ViterbiPP(NumericMatrix Ax, NumericMatrix Ay,
                NumericMatrix Ex, NumericMatrix Ey,
                NumericVector qe, int type, IntegerVector windowspace,
                double offset = 0){
  int n = Ex.ncol() + 1;
  int m = Ey.ncol() + 1;
  NumericMatrix Saa(n - 1, m - 1);
  for(int i = 0; i < n - 1; i++){
    for (int j = 0; j < m - 1; j++){
      Saa(i, j) = logsum(Ex(_, i) + Ey(_, j) - qe);
    }
  }
  double sij;
  NumericVector MIcdt(2);
  NumericVector DGcdt(2);
  NumericVector MMcdt(6);
  if(type == 2) MMcdt[5] = 0; else MMcdt[5] = -INFINITY;
  NumericVector GDcdt(2);
  NumericVector IMcdt(2);
  IntegerVector MIind = IntegerVector::create(2, 0);
  IntegerVector DGind = IntegerVector::create(2, 1);
  IntegerVector MMind = IntegerVector::create(0, 1, 2, 3, 4, 5);
  IntegerVector GDind = IntegerVector::create(2, 3);
  IntegerVector IMind = IntegerVector::create(2, 4);
  int MImax;
  int DGmax;
  int MMmax;
  int GDmax;
  int IMmax;
  NumericVector neginf = NumericVector::create(-INFINITY);
  NumericVector tmp = rep_len(neginf, n * m);
  tmp.attr("dim") = IntegerVector::create(n, m);
  NumericMatrix MImatrix = as<NumericMatrix>(tmp);
  NumericMatrix DGmatrix = clone(MImatrix);
  NumericMatrix MMmatrix = clone(MImatrix);
  NumericMatrix GDmatrix = clone(MImatrix);
  NumericMatrix IMmatrix = clone(MImatrix);
  IntegerMatrix MIpointer(n, m); // pointer matrices
  IntegerMatrix DGpointer = clone(MIpointer);
  IntegerMatrix MMpointer = clone(MIpointer);
  IntegerMatrix GDpointer = clone(MIpointer);
  IntegerMatrix IMpointer = clone(MIpointer);
  // initialize scoring matrices
  MMmatrix(0, 0) = 0;
  if(type == 0){
    MImatrix(1, 0) = Ax(4, 0) + Ay(5, 0); // MM + MI
    DGmatrix(1, 0) = Ax(3, 0); //MD
    MIpointer(1, 0) = 2;
    DGpointer(1, 0) = 2;
    for(int i = 2; i < n; i++){
      MImatrix(i, 0) = MImatrix(i - 1, 0) + Ax(4, i - 1) + Ay(8, 0); // MM + II
      DGmatrix(i, 0) = DGmatrix(i - 1, 0) + Ax(0, i - 1); // DD
      MIpointer(i, 0) = 0; // added apr 4 17
      DGpointer(i, 0) = 1;// added apr 4 17
    }
    MMmatrix(0, 0) = 0;
    GDmatrix(0, 1) = Ay(3, 0); //MD
    IMmatrix(0, 1) = Ax(5, 0) + Ay(4, 0); //MI + MM
    GDpointer(0, 1) = 2;
    IMpointer(0, 1) = 2;
    for(int j = 2; j < m; j++){
      GDmatrix(0, j) = GDmatrix(0, j - 1) + Ay(0, j - 1); // DD
      IMmatrix(0, j) = IMmatrix(0, j - 1) + Ay(4, j - 1) + Ax(8, 0); // MM + II
      GDpointer(0, j) = 3;// added apr 4 17
      IMpointer(0, j) = 4;// added apr 4 17
    }
  }else{
    for(int i = 1; i < n; i++) {
      MMmatrix(i, 0) = 0;
      MIpointer(i, 0) = 0;
      DGpointer(i, 0) = 1;
    }
    MIpointer(1, 0) = 2;
    DGpointer(1, 0) = 2;
    MMmatrix(0, 0) = 0;
    for(int j = 1; j < m; j++) {
      MMmatrix(0, j) = 0;
      GDpointer(0, j) = 3;
      IMpointer(0, j) = 4;
    }
    GDpointer(0, 1) = 2;
    IMpointer(0, 1) = 2;
  }
  // recursion step
  for(int i = 1; i < n; i++){
    for(int j = 1; j < m; j++){
      if(((j - i) >= windowspace[0]) & ((j - i) <= windowspace[1])){
        sij = Saa(i - 1, j - 1) + offset;
        MIcdt[0] = MMmatrix(i - 1, j) + Ax(4, i - 1) + Ay(5, j); //MM + MI
        MIcdt[1] = MImatrix(i - 1, j) + Ax(4, i - 1) + Ay(8, j); //MM + II
        DGcdt[0] = MMmatrix(i - 1, j) + Ax(3, i - 1); //MD
        DGcdt[1] = DGmatrix(i - 1, j) + Ax(0, i - 1); //DD
        MMcdt[0] = MImatrix(i - 1, j - 1) + Ax(4, i - 1) + Ay(7, j - 1) + sij; //MM + IM
        MMcdt[1] = DGmatrix(i - 1, j - 1) + Ax(1, i - 1) + Ay(4, j - 1) + sij; //DM + MM
        MMcdt[2] = MMmatrix(i - 1, j - 1) + Ax(4, i - 1) + Ay(4, j - 1) + sij; //MM + MM
        MMcdt[3] = GDmatrix(i - 1, j - 1) + Ax(4, i - 1) + Ay(1, j - 1) + sij; //MM + DM
        MMcdt[4] = IMmatrix(i - 1, j - 1) + Ax(7, i - 1) + Ay(4, j - 1) + sij; //IM + MM
        GDcdt[0] = MMmatrix(i, j - 1) + Ay(3, j - 1); //MD
        GDcdt[1] = GDmatrix(i, j - 1) + Ay(0, j - 1); //DD
        IMcdt[0] = MMmatrix(i, j - 1) + Ax(5, i) + Ay(4, j - 1); //MI + MM
        IMcdt[1] = IMmatrix(i, j - 1) + Ax(8, i) + Ay(4, j - 1); //II + MM
        MImax = whichmax(MIcdt, 0);
        DGmax = whichmax(DGcdt, 0);
        MMmax = whichmax(MMcdt, 0);
        GDmax = whichmax(GDcdt, 0);
        IMmax = whichmax(IMcdt, 0);
        MImatrix(i, j) = MIcdt[MImax];
        DGmatrix(i, j) = DGcdt[DGmax];
        MMmatrix(i, j) = MMcdt[MMmax];
        GDmatrix(i, j) = GDcdt[GDmax];
        IMmatrix(i, j) = IMcdt[IMmax];
        MIpointer(i, j) = MIind[MImax];
        DGpointer(i, j) = DGind[DGmax];
        MMpointer(i, j) = MMind[MMmax];
        GDpointer(i, j) = GDind[GDmax];
        IMpointer(i, j) = IMind[MImax];
      }
    }
    checkUserInterrupt();
  }
  // traceback
  IntegerVector path(m + n); // if being picky should be m + n - 2 but doesnt matter
  LogicalVector keeppath(m + n);
  int counter = m + n - 1;
  // if(type < 2){//include matchy-matchy for end states if alignment type is not local
  //   path[counter] = 2;
  //   keeppath[counter] = true;
  //   counter--;
  // }
  int tbr;
  int tbc;
  int tbm;
  double score;
  IntegerVector startposition(2);
  if(type == 0){
    NumericVector LLcdt = NumericVector::create(
      MImatrix(n - 1, m - 1) + Ax(4, n - 1) + Ay(7, m - 1), //MM + IM
      DGmatrix(n - 1, m - 1) + Ax(1, n - 1) + Ay(4, m - 1), //DM + MM
      MMmatrix(n - 1, m - 1) + Ax(4, n - 1) + Ay(4, m - 1), //MM + MM
      GDmatrix(n - 1, m - 1) + Ax(4, n - 1) + Ay(1, m - 1), //MM + DM
      IMmatrix(n - 1, m - 1) + Ax(7, n - 1) + Ay(4, m - 1));
    tbm = whichmax(LLcdt, 0);
    score = LLcdt[tbm];
    tbr = n - 1; // traceback row
    tbc = m - 1; // traceback column
    while((tbr > 0) | (tbc > 0)){
      if((tbr < 0) | (tbc < 0)) throw Rcpp::exception("Error 2");
      keeppath[counter] = true;
      path[counter] = tbm;
      if(tbm < 2){
        if(tbm == 0) tbm = MIpointer(tbr, tbc); else tbm = DGpointer(tbr, tbc);
        tbr--;
      }else if(tbm == 2){
        tbm = MMpointer(tbr, tbc);
        tbr--;
        tbc--;
      }else{
        if(tbm == 3) tbm = GDpointer(tbr, tbc); else tbm = IMpointer(tbr, tbc);
        tbc--;
      }
      counter--;
      checkUserInterrupt();
    }
  }else if(type == 1){
    tbr = n - 1;
    tbc = m - 1;
    tbm = 2;
    score = MMmatrix(tbr, tbc);
    for(int i = n - 2; i >= 0; i--){
      if(MMmatrix(i, m - 1) > score){
        tbr = i;
        score = MMmatrix(i, m - 1);
      }
    }
    for(int j = m - 1; j >= 0; j--){
      if(MMmatrix(n - 1, j) > score){
        tbc = j;
        tbr = n - 1;
        score = MMmatrix(n - 1, j);
      }
    }//could have equal option too using R::runif
    if(tbr < n - 1){
      for(int i = tbr + 1; i < n; i++){
        path[counter] = 1;
        keeppath[counter] = true;
        counter--;
      }
    }else if(tbc < m - 1){
      for(int j = tbc + 1; j < m; j++){
        path[counter] = 3;
        keeppath[counter] = true;
        counter--;
      }
    }
    while((tbr > 0) & (tbc > 0)){
      keeppath[counter] = true;
      path[counter] = tbm;
      if(tbm < 2){
        if(tbm == 0) tbm = MIpointer(tbr, tbc); else tbm = DGpointer(tbr, tbc);
        tbr--;
      }else if(tbm == 2){
        tbm = MMpointer(tbr, tbc);
        tbr--;
        tbc--;
      }else{
        if(tbm == 3) tbm = GDpointer(tbr, tbc); else tbm = IMpointer(tbr, tbc);
        tbc--;
      }
      counter--;
      checkUserInterrupt();
    }
    if(tbr > 0){
      for(int i = 0; i < tbr; i++){
        path[counter] = 1;
        keeppath[counter] = true;
        counter--;
      }
    }else if(tbc > 0){
      for(int j = 0; j < tbc; j++){
        path[counter] = 3;
        keeppath[counter] = true;
        counter--;
      }
    }
  }else{
    MMpointer(0, 0) = 5;
    int maxind = whichmax(MMmatrix, 0);
    score = MMmatrix[maxind];
    tbr = maxind % (n); // remainder
    tbc = maxind/(n); // quotient
    tbm = 2;
    bool advance = MMpointer(tbr, tbc) != 5;
    while(advance){
      keeppath[counter] = true;
      path[counter] = tbm;
      if(tbm < 2){
        if(tbm == 0){
          tbm = MIpointer(tbr, tbc);
          tbr--;
          advance = MIpointer(tbr, tbc) != 5;
        }else{
          tbm = DGpointer(tbr, tbc);
          tbr--;
          advance = DGpointer(tbr, tbc) != 5;
        }
      }else if(tbm == 2){
        tbm = MMpointer(tbr, tbc);
        tbr--;
        tbc--;
        advance = MMpointer(tbr, tbc) != 5;
      }else{
        if(tbm == 3){
          tbm = GDpointer(tbr, tbc);
          tbc--;
          advance = GDpointer(tbr, tbc) != 5;
        }else{
          tbm = IMpointer(tbr, tbc);
          tbc--;
          advance = IMpointer(tbr, tbc) != 5;
        }
      }
      counter--;
      startposition[0] = tbr;
      startposition[1] = tbc;
      checkUserInterrupt();
    }
  }
  IntegerVector truncpath = path[keeppath];
  startposition = startposition + 1; // R's indexing style, starts at 1 instead of 0
  List res = List::create(Named("score") = score,
                          Named("path") = truncpath,
                          Named("start") = startposition,
                          Named("MImatrix") = MImatrix,
                          Named("DGmatrix") = DGmatrix,
                          Named("MMmatrix") = MMmatrix,
                          Named("GDmatrix") = GDmatrix,
                          Named("IMmatrix") = IMmatrix,
                          Named("MIpointer") = MIpointer,
                          Named("DGpointer") = DGpointer,
                          Named("MMpointer") = MMpointer,
                          Named("GDpointer") = GDpointer,
                          Named("IMpointer") = IMpointer,
                          Named("Saa") = Saa);
  res.attr("class") = "DPA";
  return(res);
}

// [[Rcpp::export(name = ".forwardH")]]
List forwardH(IntegerVector y, NumericMatrix A, NumericMatrix E,
                 bool DNA = false, bool AA = false){
  int nrolls = y.size();
  IntegerVector Edim = E.attr("dim");
  int nstates = Edim[0];
  //int nres = Edim[1];
  if(nrolls == 0){
    List out = List::create(Named("score") = A(0, 0), Named("array") = NULL);
    out.attr("class") = "DPA";
    return out;
  }
  NumericMatrix R(nstates, nrolls);
  NumericVector s = A(0, _);
  NumericVector e = A(_, 0);
  NumericMatrix tmp(nstates, nstates);
  NumericVector coltotals(nstates);
  if(DNA){
    for(int l = 0; l < nstates; l++) R(l, 0) = probDNA(y[0], E(l,_)) + A(0, l + 1);
    for(int i = 1; i < nrolls; i++){
      for(int j = 0; j < nstates; j++){
        for(int k = 0; k < nstates; k++){
          tmp(j, k) = R(j, i - 1) + A(j + 1, k + 1);
        }
      }
      for(int l = 0; l < nstates; l++){
        R(l, i) = probDNA(y[i], E(l,_)) + logsum(tmp(_, l));
      }
    }
  }else if(AA){
    for(int l = 0; l < nstates; l++) R(l, 0) = probAA(y[0], E(l,_)) + A(0, l + 1);
    for(int i = 1; i < nrolls; i++){
      for(int j = 0; j < nstates; j++){
        for(int k = 0; k < nstates; k++){
          tmp(j, k) = R(j, i - 1) + A(j + 1, k + 1);
        }
      }
      for(int l = 0; l < nstates; l++){
        R(l, i) = probAA(y[i], E(l,_)) + logsum(tmp(_, l));
      }
    }
  }else{
    R(_, 0) = E(_, y[0]) + s[seq(1, nstates)];
    for(int i = 1; i < nrolls; i++){
      for(int k = 0; k < nstates; k++){
        for(int l = 0; l < nstates; l++){
          tmp(k, l) = R(k, i - 1) + A(k + 1, l + 1);
        }
      }
      for(int k = 0; k < nstates; k++) coltotals[k] = logsum(tmp(_, k));
      R(_, i) =  E(_, y[i]) + coltotals;
    }
  }
  NumericVector ak0(nstates);
  if(any(e[seq(1, nstates)] != rep(-INFINITY, nstates)).is_true()) ak0 = e[seq(1, nstates)];
  double res = logsum(R(_, nrolls - 1) + ak0);
  bool odds = false;
  List out = List::create(Named("score") = res, Named("array") = R, Named("odds") = odds);
  out.attr("class") = "DPA";
  return out;
}

// [[Rcpp::export(name = ".forwardP")]]
List forwardP(IntegerVector y, NumericMatrix A, NumericMatrix E, NumericVector qe,
                  NumericVector qey, IntegerVector windowspace,
                  bool DI = false, bool ID = false, bool DNA = false, bool AA = false){
  int n = E.ncol() + 1;
  int m = y.size() + 1;
  double sij;
  NumericVector Dcdt(3);
  Dcdt[2] = -INFINITY;
  NumericVector Mcdt(3);
  NumericVector Icdt(3);
  Icdt[0] = -INFINITY;
  NumericVector neginf = NumericVector::create(-INFINITY);
  NumericVector tmp = rep_len(neginf, n * m);
  tmp.attr("dim") = IntegerVector::create(n, m);
  NumericMatrix Dmatrix = as<NumericMatrix>(tmp);// x aligns to gap in y
  NumericMatrix Mmatrix = clone(Dmatrix);// match-match
  NumericMatrix Imatrix = clone(Dmatrix); // y aligns to gap in x
  // initialize scoring matrices
  Mmatrix(0, 0) = 0;
  Dmatrix(1, 0) = A(3, 0); // Match -> Delete transition prob for position 0
  Imatrix(0, 1) = A(5, 0) + qey[0];
  for(int i = 2; i < n; i++) Dmatrix(i, 0) = Dmatrix(i - 1, 0) + A(0, i - 1);
  for(int j = 2; j < m; j++) Imatrix(0, j) = Imatrix(0, j - 1) + A(8, 0) + qey[j - 1];
  // recursion
  for(int i = 1; i < n; i++){
    for(int j = 1; j < m; j++){
      if(((j - i) >= windowspace[0]) & ((j - i) <= windowspace[1])){
        if(DNA){
          sij = probDNA(y[j - 1], E(_, i - 1));
        }else if(AA){
          sij = probAA(y[j - 1], E(_, i - 1));
        }else{
          sij = E(y[j - 1], i - 1);
        }
        Dcdt[0] = Dmatrix(i - 1, j) + A(0, i - 1); //DD
        Dcdt[1] = Mmatrix(i - 1, j) + A(3, i - 1); //MD
        if(ID) Dcdt[2] = Imatrix(i - 1, j) + A(6, i - 1);  //ID
        Mcdt[0] = Dmatrix(i - 1, j - 1) + A(1, i - 1); //DM
        Mcdt[1] = Mmatrix(i - 1, j - 1) + A(4, i - 1); //MM
        Mcdt[2] = Imatrix(i - 1, j - 1) + A(7, i - 1); //IM
        if(DI) Icdt[0] = Dmatrix(i, j - 1) + A(2, i); //DI
        Icdt[1] = Mmatrix(i, j - 1) + A(5, i); //MI
        Icdt[2] = Imatrix(i, j - 1) + A(8, i); //II
        Dmatrix(i, j) = logsum(Dcdt);
        Mmatrix(i, j) = logsum(Mcdt) + sij;
        Imatrix(i, j) = logsum(Icdt) + qey[j - 1];
      }
    }
    checkUserInterrupt();
  }
  double score;
  NumericVector LLcdt = NumericVector::create(Dmatrix(n - 1, m - 1) + A(1, n - 1), //DM
                                              Mmatrix(n - 1, m - 1) + A(4, n - 1), //MM
                                              Imatrix(n - 1, m - 1) + A(7, n - 1)); //IM
  score = logsum(LLcdt);
  bool odds = all(qey == 0).is_true();
  List res = List::create(Named("score") = score,
                          Named("odds") = odds,
                          Named("Dmatrix") = Dmatrix,
                          Named("Mmatrix") = Mmatrix,
                          Named("Imatrix") = Imatrix);
  res.attr("class") = "DPA";
  return(res);
}

// [[Rcpp::export(name = ".backwardH")]]
List backwardH(IntegerVector y, NumericMatrix A, NumericMatrix E,
                  bool DNA = false, bool AA = false) {
  int nrolls = y.size();
  IntegerVector Edim = E.attr("dim");
  int nstates = Edim[0]; // does not include begin/end state
  //int nres = Edim[1];
  if(nrolls == 0){
    List out = List::create(Named("score") = A(0, 0), Named("array") = NULL);
    out.attr("class") = "DPA";
    return out;
  }
  NumericMatrix R(nstates, nrolls);
  NumericVector s = A(0, _);
  NumericVector e = A(_, 0);
  NumericVector ak0(nstates);
  if(any(e[seq(1, nstates)] != rep(-INFINITY, nstates)).is_true()) ak0 = e[seq(1, nstates)];
  R(_, nrolls - 1) = ak0;
  NumericMatrix tmp(nstates, nstates);
  NumericVector logprobs(nstates);
  if(DNA){
    for(int i = nrolls - 1; i > 0; i--){
      for(int k = 0; k < nstates; k++){
        for(int l = 0; l < nstates; l++){
          tmp(k, l) =  A(k + 1, l + 1) + probDNA(y[i], E(l,_)) + R(l, i);
        }
      }
      for(int k = 0; k < nstates; k++) R(k, i - 1) = logsum(tmp(k, _));
    }
    for(int l = 0; l < nstates; l++) {
      logprobs[l] = A(0, l + 1) + probDNA(y[0], E(l,_)) + R(l, 0);
    }
  }else if(AA){
    for(int i = nrolls - 1; i > 0; i--){
      for(int k = 0; k < nstates; k++){
        for(int l = 0; l < nstates; l++){
          tmp(k, l) =  A(k + 1, l + 1) + probAA(y[i], E(l,_)) + R(l, i);
        }
      }
      for(int k = 0; k < nstates; k++) R(k, i - 1) = logsum(tmp(k, _));
    }
    for(int l = 0; l < nstates; l++) {
      logprobs[l] = A(0, l + 1) + probAA(y[0], E(l,_)) + R(l, 0);
    }
  }else{
    for(int i = nrolls - 1; i > 0; i--){
      for(int k = 0; k < nstates; k++){
        for(int l = 0; l < nstates; l++){
          tmp(k, l) =  A(k + 1, l + 1) + E(l, y[i]) + R(l, i);
        }
      }
      for(int k = 0; k < nstates; k++) R(k, i - 1) = logsum(tmp(k, _));
    }
    for(int l = 0; l < nstates; l++) {
      logprobs[l] = A(0, l + 1) + E(l, y[0]) + R(l, 0);
    }
  }
  double res = logsum(logprobs);
  bool odds = false;
  List out = List::create(Named("score") = res, Named("array") = R, Named("odds") = odds);
  out.attr("class") = "DPA";
  return out;
}

// [[Rcpp::export(name = ".backwardP")]]
List backwardP(IntegerVector y, NumericMatrix A, NumericMatrix E, NumericVector qe,
                  NumericVector qey, IntegerVector windowspace,
                  bool DI = false, bool ID = false, bool DNA = false, bool AA = false){
  int n = E.ncol() + 1;
  int m = y.size() + 1;
  double sij;
  NumericVector Dcdt(3);
  Dcdt[2] = -INFINITY;
  NumericVector Mcdt(3);
  NumericVector Icdt(3);
  Icdt[0] = -INFINITY;
  NumericVector neginf = NumericVector::create(-INFINITY);
  NumericVector tmp = rep_len(neginf, n * m);
  tmp.attr("dim") = IntegerVector::create(n, m);
  NumericMatrix Dmatrix = as<NumericMatrix>(tmp);// x aligns to gap in y
  NumericMatrix Mmatrix = clone(Dmatrix);// match-match
  NumericMatrix Imatrix = clone(Dmatrix); // y aligns to gap in x
  // initialize scoring matrices
  Dmatrix(n - 1, m - 1) = A(1, n - 1);
  Mmatrix(n - 1, m - 1) = A(4, n - 1);
  Imatrix(n - 1, m - 1) = A(7, n - 1);
  for(int i = n - 1; i > 0; i--) {
    Dmatrix(i - 1, m - 1) = Dmatrix(i, m - 1) + A(0, i - 1); //DD
    Mmatrix(i - 1, m - 1) = Dmatrix(i, m - 1) + A(3, i - 1); //MD
    if(ID) Imatrix(i - 1, m - 1) = Dmatrix(i, m - 1) + A(6, i - 1); //ID
  }
  for(int j = m - 1; j > 0; j--){
    if(DI) Dmatrix(n - 1, j - 1) = Imatrix(n - 1, j) + A(2, n - 1) + qey[j - 1]; //DI
    Mmatrix(n - 1, j - 1) = Imatrix(n - 1, j) + A(5, n - 1) + qey[j - 1]; //MI
    Imatrix(n - 1, j - 1) = Imatrix(n - 1, j) + A(8, n - 1) + qey[j - 1]; //II
  }
  // recursion
  for(int i = n - 2; i >= 0; i--){
    for(int j = m - 2; j >= 0; j--){
      if(((j - i) >= windowspace[0]) & ((j - i) <= windowspace[1])){
        if(DNA){
          sij = probDNA(y[j], E(_, i));
        }else if(AA){
          sij = probAA(y[j], E(_, i));
        }else{
          sij = E(y[j], i);
        }
        Dcdt[0] = Dmatrix(i + 1, j) + A(0, i); //DD
        Dcdt[1] = Mmatrix(i + 1, j + 1) + A(1, i) + sij; //DM
        if(DI) Dcdt[2] = Imatrix(i, j + 1) + A(2, i) + qey[j];  //DI
        Mcdt[0] = Dmatrix(i + 1, j) + A(3, i); //MD
        Mcdt[1] = Mmatrix(i + 1, j + 1) + A(4, i) + sij; //MM
        Mcdt[2] = Imatrix(i, j + 1) + A(5, i) + qey[j]; //IM
        if(DI) Icdt[0] = Dmatrix(i + 1, j) + A(6, i); //DI
        Icdt[1] = Mmatrix(i + 1, j + 1) + A(7, i) + sij; //MI
        Icdt[2] = Imatrix(i, j + 1) + A(8, i) + qey[j]; //II
        Dmatrix(i, j) = logsum(Dcdt);
        Mmatrix(i, j) = logsum(Mcdt);
        Imatrix(i, j) = logsum(Icdt);
      }
    }
    checkUserInterrupt();
  }
  double score;
  score = Mmatrix(0, 0);
  bool odds = all(qey == 0).is_true();
  List res = List::create(Named("score") = score,
                          Named("odds") = odds,
                          Named("Dmatrix") = Dmatrix,
                          Named("Mmatrix") = Mmatrix,
                          Named("Imatrix") = Imatrix);
  res.attr("class") = "DPA";
  return(res);
}
