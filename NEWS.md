# aphid 1.3.0 (2018-10-09)

* Enabled 'near-enough' alignment matching between Viterbi iterations in `train.PHMM` via `limit` parameter

* Prevented memory overflow in 'train.PHMM' caused by md5 hash attmpt on very large alignments 

* Included faster maximum entropy  sequence weighting option (Henikoff & Henikoff 1994)

--------------------------------------------------------------------------------

# aphid 1.2.0 (2018-08-14)

* Enabled multithread option for Baum Welch training of profile HMMs.

* Fixed checksum parse issue in readPHMM 

--------------------------------------------------------------------------------

# aphid 1.1.0

* Added posterior probabilities for PHMMs.

* derivePHMM.list seed selection changed from max length to most common length.

* Fixed issue in 'align' causing distance matrix error if seqweights not 
set to NULL (Thanks to Ben Margetts for bug report).

* Updated dependency list to include kmer instead of phylogram.

--------------------------------------------------------------------------------

# aphid 1.0.1

* Internalized function erroneously exported in initial CRAN release.

--------------------------------------------------------------------------------

# aphid 1.0.0

* Released on CRAN 2017-06-24.
