#!/bin/bash

all_unaligned_path=data/unaligned_all.fasta
backbone_aln_path=data/backbone.aln.fasta.gz
backbone_tre_path=data/backbone.tre
query_path=data/unaligned_frag.fasta
outname=aligned.fasta

scenario=4
if [[ $1 != "" ]]; then
    scenario=$1
fi

if [[ $scenario == 1 ]]; then
    # Scenario A - unaligned sequences only
    python ../witch.py -i ${all_unaligned_path} -d scenarioA_output \
        -o ${outname}
elif [[ $scenario == 2 ]]; then
    # Scenario B - unaligned sequences only; using bit scores;
    #              using 10 HMMs to align a sequence
    python ../witch.py -i ${all_unaligned_path} -d scenarioB_output \
        -o ${outname} -w 0 -k 10
elif [[ $scenario == 3 ]]; then
    # 3) Scenario C - backbone alignment available; backbone tree missing;
    #                 query sequences available; also saving weights to local
    python ../witch.py --num-cpus -1 -b ${backbone_aln_path} \
        -q ${query_path} -d scenarioC_output -o ${outname} \
        --save-weight 1
elif [[ $scenario == 4 ]]; then
    # 4) Scenario D - backbone alignment available; backbone tree available;
    #                 query sequences available; saving weights to local;
    #                 also save decomposition results for future use (e.g.,
    #                 faster rerun)
    python ../witch.py --num-cpus -1 -b ${backbone_aln_path} \
        -e ${backbone_tre_path} \
        -q ${query_path} -d scenarioD_output -o ${outname} \
        --save-weight 1 --keep-decomposition 1
elif [[ $scenario == 5 ]]; then
    # 5) Scenario E - same as Scenario D, but with a user-specified config file
    python ../witch.py --num-cpus -1 -b ${backbone_aln_path} \
        -e ${backbone_tre_path} \
        -q ${query_path} -d scenarioE_output -o ${outname} \
        --save-weight 1 --keep-decomposition 1 \
        -c user.config
fi
