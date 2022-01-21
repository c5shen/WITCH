#!/bin/bash
cd ../

## Scenario A - unaligned sequences only
#python3 witch.py -i examples/data/unaligned_all.txt -d scenarioA_output -o aligned.txt

## Scenario B - unaligned sequences only; using bit scores;
##              using 10 HMMs to align a sequence
#python3 witch.py -i examples/data/unaligned_all.txt -d scenarioB_output -o aligned.txt -w 0 -k 10

## 3) Scenario C - backbone alignment available; backbone tree missing;
##                 query sequences available 
#python3 witch.py -b examples/data/backbone.aln.fasta -q examples/data/unaligned_frag.txt -d scenarioC_output -o aligned.txt
