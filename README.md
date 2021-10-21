-----------------------------
Summary
-----------------------------

GCM+eHMMs is a new multiple sequence alignment (MSA) tool that combines techniques from [UPP](https://github.com/smirarab/sepp) and [MAGUS](https://github.com/vlasmirnov/MAGUS).

GCM+eHMMs aims to solve alignment problems when input sequences contain fragments or show high levels of sequence length hetereogeneity. 


Ideally, GCM+eHMMs wants to solve the following problem:
- Input: a set of unaligned sequences `S`
- Output an alignment `A` on `S`

**(Oct 21st 2021) Current implementation of GCM+eHMMs only supports the following problem. We are working on making a full pipeline for GCM+eHMMs:**
- Input: a set of unaligned full-length sequences and a set of unaligned fragmentary sequences.
- Intermediate output:
        * a backbone alignment/tree on at most 1,000 full-length sequences
        * a set of HMMs 
