#!/bin/bash

# Scenario A - unaligned sequences only
#python3 ../witch.py -i data/unaligned_all.txt -d scenarioA_output -o aligned.txt

## Scenario B - unaligned sequences only; using bit scores;
##              using 10 HMMs to align a sequence
#python3 ../witch.py -i data/unaligned_all.txt -d scenarioB_output -o aligned.txt -w 0 -k 10

## 3) Scenario C - backbone alignment available; backbone tree missing;
##                 query sequences available; also saving weights to local
#python3 ../witch.py --num-cpus -1 -b data/backbone.aln.fasta \
#    -q data/unaligned_frag.txt -d scenarioC_output -o aligned.txt \
#    --save-weight 1

# 4) Scenario D - backbone alignment available; backbone tree available;
#                 query sequences available; saving weights to local;
#                 also saving working (temp) files to examine
python3 ../witch.py --num-cpus -1 -b data/backbone.aln.fasta \
    -e data/backbone.tre \
    -q data/unaligned_frag.txt -d scenarioD_output -o aligned.txt \
    --save-weight 1 --keeptemp --keepgcmtemp
