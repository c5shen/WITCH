#!/bin/bash

# weighted, k=4
#python3 ../gcm+eHMMs.py \
#        -p ~/tallis/playground/eHMMs_group/benchmark_results/upp/A\=10/high_frag/1000M1/R0/temp_pasta_upp_align/ \
#        -b ../../R0/backbone.aln.fasta -q ../../R0/unaligned_frag.txt \
#        -w -k 4 -o 1000M1_hf_R0_weighted/k=4/ 

# weighted, k=4, normalized
python3 ../gcm+eHMMs.py \
        -p ~/tallis/playground/eHMMs_group/benchmark_results/upp/A\=10/high_frag/1000M1/R0/temp_pasta_upp_align/ \
        -b ../../R0/backbone.aln.fasta -q ../../R0/unaligned_frag.txt \
        -w -k 4 --weight-adjust normalize -o 1000M1_hf_R0_weighted_normalize/k=4/ 
