#include <Rcpp.h>
using namespace Rcpp;

// x is a list of sequences, can be a DNAbin object
// k is the tuple size (non-negative integer)

// [[Rcpp::export(name = "kdist.default")]]
NumericVector kdist(List x, int k = 5, bool asmatrix = false) {
  CharacterVector nms = x.attr("names");
  int nseq = x.size();
  IntegerVector seqlens(nseq);
  int mergedseqsize = 0;
  for(int i = 0; i < nseq; i++){
    seqlens[i] = Rf_length(x[i]);
    mergedseqsize += seqlens[i];
  }
  CharacterVector mergedseq(mergedseqsize);
  int counter = 0;
  for(int i = 0; i < nseq; i++){
    CharacterVector seqi = VECTOR_ELT(x, i);
    for(int j = 0; j < seqlens[i]; j++){
      mergedseq[counter] = seqi[j];
      counter++;
    }
  }
  CharacterVector uniques = unique(mergedseq);
  int alphasize = uniques.size();
  double astpk = pow(alphasize, k); // alphasize to power k
  NumericMatrix tuplemat(nseq, astpk);
  IntegerVector rsak(k); //reversed seq along k
  counter = 0;
  for(int l = k - 1; l >= 0; l--) {
    rsak[l] = pow(alphasize, counter);
    counter++;
  }
  int decimalindex;
  IntegerVector srtindex(1);
  for(int i = 0; i < nseq; i++){
    CharacterVector seqi = VECTOR_ELT(x, i);
    for(int j = 0; j < seqlens[i] - k; j++){
      decimalindex = 0;
      for(int l = 0; l < k; l++){
        counter = 0;
        while(counter < alphasize){
          if(seqi[j + l] == uniques[counter]){
            srtindex[0] = counter;
            break;
          }else{
            counter++;
          }
        }
        decimalindex += srtindex[0] * rsak[l];
      }
      tuplemat(i, decimalindex)++;
    }
  }
  double idenom;
  double jdenom;
  double SXY;
  NumericMatrix resm(nseq, nseq);
  NumericVector resv((pow(nseq, 2) - nseq)/2);
  if(asmatrix){
    resm.attr("dimnames") = List::create(nms, nms);
  }else{
    resv.attr("Size") = nseq;
    resv.attr("Labels") = nms;
    resv.attr("Diag") = false;
    resv.attr("Upper") = false;
    // res.attr("call") = match_call();
    resv.attr("method") = "YZ08";
    resv.attr("class") = "dist";
    counter = 0;
  }
  for(int i = 0; i < nseq; i++){
    //RawVector seqi = VECTOR_ELT(x, i);
    //idenom = seqi.size() - k + 1;
    for(int j = 1; j < nseq; j++){
      if(i < j){
        //RawVector seqj = VECTOR_ELT(x, j);
        //jdenom = seqj.size() - k + 1;
        // denom = std::min(seqi.size(), seqj.size()) - k + 1;
        // Rcout << denom;
        SXY = 0;
        for(int l = 0; l < astpk; l++){
          SXY += pow((tuplemat(i, l)/(seqlens[i] - k + 1)) -
            (tuplemat(j, l)/(seqlens[i] - k + 1)), 2);
        }
        if(asmatrix) {
          resm(i, j) = SXY;
          resm(j, i) = SXY;
        } else{
          resv[counter] = SXY;
          counter += 1;
        }
      }
    }
  }
  if(asmatrix){
    return(resm);
  }else{
    return(resv);
  }
}











// [[Rcpp::export(name = "kdist.DNAbin")]]
NumericVector kdistDNA(List x, int k = 5, bool asmatrix = false) {
  CharacterVector nms = x.attr("names");
  int nseq = x.size();
  double fourtopowerk = pow(4, k);
  NumericMatrix tuplemat(nseq, fourtopowerk);
  RawVector tmp = RawVector::create(4, 55, 0, 136, 72, 199, 40, 24,
                                    224, 176, 208, 112, 240, 8, 160,
                                    144, 96, 80);
  IntegerVector rsak(k); //reversed seq along k
  int counter = 0;
  for(int l = k - 1; l >= 0; l--) {
    rsak[l] = pow(4, counter);
    counter++;
  }
  IntegerVector seqlens(nseq);
  RawVector kmer(k);
  LogicalVector knownbase(k);
  LogicalVector isN(k);
  int decimalindex;
  IntegerMatrix vecs(4, k);
  IntegerVector veclengths(k);
  IntegerVector curseq(k); // current kmer sequence
  IntegerVector curseqind(k);
  double ncombos = 1;
  double kNs = 1/fourtopowerk;
  double seqcontrib; // contribution of current kmer
  int curk; // current x position in current kmer (ie wheel of combolock)
  bool advance;
  IntegerVector test(1);
  for(int i = 0; i < nseq; i++){
    RawVector seqi = VECTOR_ELT(x, i);
    seqlens[i] = seqi.size();
    for(int j = 0; j < seqlens[i] - k; j++){
      for(int l = 0; l < k; l++){
        kmer[l] = seqi[j + l];
        knownbase[l] = (seqi[j + l] & tmp[13]) == tmp[13];
        isN[l] = seqi[j + l] == tmp[12];
      }

      //if(j == 0) {Rcout << kmer; Rcout << "_";}

      if(all(knownbase).is_true()){
        decimalindex = 0;
        // a = 0, c = 1, g = 2, t = 3
        for(int l = 0; l < k; l++){
          //if(kmer[l] == tmp[3]){kmernum[l] = 0;}else... not needed as equal to 0
          if(kmer[l] == tmp[6]){
            decimalindex += rsak[l];
          }else if(kmer[l] == tmp[4]){
            decimalindex += 2 * rsak[l];
          }else if(kmer[l] == tmp[7]){
            decimalindex += 3 * rsak[l];
          }
        }
        tuplemat(i, decimalindex)++;
      }else if(all(isN).is_true()){
        for(int m = 0; m < fourtopowerk; m++) tuplemat(i, m) += kNs;
      }else{
        ncombos = 1;
        decimalindex = 0;
        for(int l = 0; l < k; l++){
          if((kmer[l] & tmp[1]) == tmp[2]){ // is purine?
            if(kmer[l] == tmp[3]){
              vecs(0, l) = 0;
              veclengths[l] = 1;
              //return(probs[0]); // A
            } else if(kmer[l] == tmp[4]){
              vecs(0, l) = 2;
              decimalindex += 2 * rsak[l];
              veclengths[l] = 1;
              //return(probs[2]); // G
            } else{
              vecs(0, l) = 0;
              vecs(1, l) = 2;
              veclengths[l] = 2;
              ncombos *= 2;
              //return(log((exp(probs[0]) + exp(probs[2]))/2)); // A or G
            }
          } else if((kmer[l] & tmp[5]) == tmp[2]){ // is pyrimidine?
            if(kmer[l] == tmp[6]){
              vecs(0, l) = 1;
              decimalindex += rsak[l];
              veclengths[l] = 1;
              //return(probs[1]); // C
            } else if(kmer[l] == tmp[7]){
              vecs(0, l) = 3;
              decimalindex += 3 * rsak[l];
              veclengths[l] = 1;
              //return(probs[3]); // T
            } else{
              vecs(0, l) = 1;
              decimalindex += rsak[l];
              vecs(1, l) = 3;
              veclengths[l] = 2;
              ncombos *= 2;
              // C or T
            }
          }else if(kmer[l] == tmp[14]){
            vecs(0, l) = 0;
            vecs(1, l) = 1;
            veclengths[l] = 2;
            ncombos *= 2;
            // M (A or C)
          }else if(kmer[l] == tmp[15]){
            vecs(0, l) = 0;
            vecs(1, l) = 3;
            veclengths[l] = 2;
            ncombos *= 2;
            // W (A or T)
          }else if(kmer[l] == tmp[16]){
            vecs(0, l) = 1;
            decimalindex += rsak[l];
            vecs(1, l) = 2;
            veclengths[l] = 2;
            ncombos *= 2;
            // S (G or C)
          }else if(kmer[l] == tmp[17]){
            vecs(0, l) = 2;
            decimalindex += 2 * rsak[l];
            vecs(1, l) = 3;
            veclengths[l] = 2;
            ncombos *= 2;
            // K (G or T)
          }else if(kmer[l] == tmp[8]){
            vecs(0, l) = 0;
            vecs(1, l) = 1;
            vecs(2, l) = 2;
            veclengths[l] = 3;
            ncombos *= 3;
            // V (A or C or G)
          } else if(kmer[l] == tmp[9]){
            vecs(0, l) = 0;
            vecs(1, l) = 1;
            vecs(2, l) = 3;
            veclengths[l] = 3;
            ncombos *= 3;
            // H (A or C or T)
          } else if(kmer[l] == tmp[10]){
            vecs(0, l) = 0;
            vecs(1, l) = 2;
            vecs(2, l) = 3;
            veclengths[l] = 3;
            ncombos *= 3;
            // D (A or G or T)
          } else if(kmer[l] == tmp[11]){
            vecs(0, l) = 1;
            decimalindex += rsak[l];
            vecs(1, l) = 2;
            vecs(2, l) = 3;
            veclengths[l] = 3;
            ncombos *= 3;
            // B (C or G or T)
          } else if(kmer[l] == tmp[12]){
            vecs(0, l) = 0;
            vecs(1, l) = 1;
            vecs(2, l) = 2;
            vecs(3, l) = 3;
            veclengths[l] = 4;
            ncombos *= 4;
            // N
          } else {
            throw Rcpp::exception("Invalid byte");
          }
          curseq[l] = vecs(0, l);
          curseqind[l] = 0;
        }
        seqcontrib = 1/ncombos;
        tuplemat(i, decimalindex) += seqcontrib;
        for(int m = 0; m < ncombos - 1; m++){
          decimalindex = 0;
          curk = k - 1;
          advance = true;
          while(advance){
            if(curseqind[curk] < veclengths[curk] - 1){
              curseqind[curk]++;
              if(curseqind[curk] > test[0]){
                test[0] = curseqind[curk];
              }
              if(curk < k - 1) for(int n = curk + 1; n < k; n++) curseqind[n] = 0;
              for(int l = 0; l < k; l++) decimalindex += vecs(curseqind[l], l) * rsak[l];
              tuplemat(i, decimalindex) += seqcontrib;
              advance = false;
              checkUserInterrupt();////
            }else{
              curk -= 1;
            }
          }
        }
      }
    }
  }
  double idenom;
  double jdenom;
  double SXY;
  NumericMatrix resm(nseq, nseq);
  NumericVector resv((pow(nseq, 2) - nseq)/2);
  if(asmatrix){
    resm.attr("dimnames") = List::create(nms, nms);
  }else{
    resv.attr("Size") = nseq;
    resv.attr("Labels") = nms;
    resv.attr("Diag") = false;
    resv.attr("Upper") = false;
    // res.attr("call") = match_call();
    resv.attr("method") = "YZ08";
    resv.attr("class") = "dist";
    counter = 0;
  }
  for(int i = 0; i < nseq; i++){
    //RawVector seqi = VECTOR_ELT(x, i);
    //idenom = seqi.size() - k + 1;
    for(int j = 1; j < nseq; j++){
      if(i < j){
        //RawVector seqj = VECTOR_ELT(x, j);
        //jdenom = seqj.size() - k + 1;
        // denom = std::min(seqi.size(), seqj.size()) - k + 1;
        // Rcout << denom;
        SXY = 0;
        for(int l = 0; l < fourtopowerk; l++){
          SXY += pow((tuplemat(i, l)/(seqlens[i] - k + 1)) -
            (tuplemat(j, l)/(seqlens[i] - k + 1)), 2);
        }
        if(asmatrix) {
          resm(i, j) = SXY;
          resm(j, i) = SXY;
        } else{
          resv[counter] = SXY;
          counter += 1;
        }
      }
    }
  }
  if(asmatrix){
    return(resm);
  }else{
    return(resv);
  }
}


// int denom;
// double F;
// for(int i = 1; i < nseq; i++){
//   RawVector seqi = VECTOR_ELT(x, i);
//   for(int j = 0; j < nseq; j++){
//     if(i > j){
//       RawVector seqj = VECTOR_ELT(x, j);
//       if(distance == YZ08)
//         denom = std::min(seqi.size(), seqj.size()) - k + 1;
//       // Rcout << denom;
//       F = 0;
//       for(int l = 0; l < fourtopowerk; l++){
//         F += std::min(tuplemat(i, l), tuplemat(j, l))/denom;
//       }
//       res(i, j) = log(0.1 + F);
//     }
//   }
// }
// return(res);
// }
