#!/bin/bash

scenario=4
if [[ $1 != "" ]]; then
    scenario=$1
fi

if [[ $scenario == 1 ]]; then
    # Scenario A - unaligned sequences only; running WITCH for the first
    #              time and bypassing the prompt to set up the configuration
    #              file
    python3 ../witch.py -y -i data/unaligned_all.fasta -d scenarioA_output \
        -o aligned.fasta
elif [[ $scenario == 2 ]]; then
    # Scenario B - unaligned sequences only; using bit scores;
    #              using 10 HMMs to align a sequence
    python3 ../witch.py -i data/unaligned_all.fasta -d scenarioB_output \
        -o aligned.fasta -w 0 -k 10
elif [[ $scenario == 3 ]]; then
    # 3) Scenario C - backbone alignment available; backbone tree missing;
    #                 query sequences available; also saving weights to local
    python3 ../witch.py --num-cpus -1 -b data/backbone.aln.fasta \
        -q data/unaligned_frag.fasta -d scenarioC_output -o aligned.fasta \
        --save-weight 1
elif [[ $scenario == 4 ]]; then
    # 4) Scenario D - backbone alignment available; backbone tree available;
    #                 query sequences available; saving weights to local;
    #                 also save decomposition results for future use (e.g.,
    #                 faster rerun)
    python3 ../witch.py --num-cpus -1 -b data/backbone.aln.fasta \
        -e data/backbone.tre \
        -q data/unaligned_frag.fasta -d scenarioD_output -o aligned.fasta \
        --save-weight 1 --keep-decomposition 1
elif [[ $scenario == 5 ]]; then
    # 5) Scenario E - same as Scenario D, but with a user-specified config file
    python3 ../witch.py --num-cpus -1 -b data/backbone.aln.fasta \
        -e data/backbone.tre \
        -q data/unaligned_frag.fasta -d scenarioE_output -o aligned.fasta \
        --save-weight 1 --keep-decomposition 1 \
        -c user.config
fi
