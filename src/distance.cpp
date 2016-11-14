#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix kcount_DNA(List x, int k = 5) {
  CharacterVector nms = x.attr("names");
  int nseq = x.size();
  double fourtopowerk = pow(4, k);
  NumericMatrix tuplemat(nseq, fourtopowerk);
  tuplemat.attr("dimnames") = List::create(nms, CharacterVector(fourtopowerk));
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
  return(tuplemat);
}


