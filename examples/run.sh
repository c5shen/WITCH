#!/bin/bash

# WITCH settings: -w 1 (using weights)
#                 -k 4 (using 4 HMMs to align each query sequence)
#                 -d witch_output (output directory)
#                 -o aligned.txt (alignment file, will be in output directory)

## 1) input: unaligned sequences only
#python3 ../witch.py -i data/unaligned_all.txt \
#        -w 1 -k 4 -d witch_output -o aligned.txt

# 2) input: a backbone alignment (but no backbone tree),
#           and a set of query sequences
python3 ../witch.py -b data/backbone.aln.fasta \
        -q data/unaligned_frag.txt \
        -w 1 -k 4 -d witch_output -o aligned.txt

# 3) input: an ensemble of HMMs generated from a backbone alignment/tree,
#           the backbone alignment, and a set of query sequences 
#python3 ../witch.py -b data/backbone.aln.fasta \
#        -q data/unaligned_frag.txt \
#        -p data/hmmdir \
#        -w 1 -k 4 -d witch_output -o aligned.txt
