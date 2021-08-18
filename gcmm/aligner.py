import os
import time
from collections import defaultdict
from configs import Configs
from gcmm.weighting import Weights
from helpers.alignment_tools import Alignment, ExtendedAlignment
from multiprocessing import Lock

'''
Helper function to generate backbone alignments for the constraint sets
'''
def getBackbones(index_to_hmm, unaligned, workdir, backbone_dir):
    # return is a map from HMM index to a list of fragments that have high
    # bit scores at the HMM.
    ret = defaultdict(list)
    if Configs.use_weight:
        weights = Weights.weights 
    else:
        weights = Weights.ranked_bitscores

    # for each taxon and its scores on HMMs, sort it in decreasing order
    for taxon, sorted_weights in weights.items():
        if taxon not in unaligned:
            continue

        top_k_hmms = sorted_weights[:Configs.num_hmms]
        if Configs.use_weight:
            if Configs.weight_adjust == 'normalize':
                cur_total_w = sum([w[1] for w in top_k_hmms])
                top_k_hmms = [(w[0], w[1] * (1. / cur_total_w))
                        for w in top_k_hmms]
            elif Configs.weight_adjust == 'maxto1':
                max_w = top_k_hmms[0][1]
                top_k_hmms = [(w[0], w[1] / max_w) for w in top_k_hmms]
        else:
            top_k_hmms = [(w[0], 1) for w in top_k_hmms]
        Configs.log('weights for {}: {}'.format(taxon, top_k_hmms))
        for item in top_k_hmms:
            # append taxon to corresponding HMM i (defined by the index)
            ret[item[0]].append(taxon)

    # now we can split the fragments into chunks, each chunk corresponds to
    # all fragments added to an HMM i
    # For each HMM, run HMMAlign on its fragment chunk to get a backbone
    # alignment
    weights_path = workdir + '/weights.txt'
    weights_file = open(weights_path, 'w')
    
    # initialize the backbone index at backbone_start_index
    bb_index = 0
    for i, names in ret.items():
        #Configs.log("Generating fragment chunks/alignment for HMM {}".format(i))
        hmm_dir = workdir
        if not os.path.isdir(hmm_dir):
            os.system('mkdir -p {}'.format(hmm_dir))
        if os.path.isdir(hmm_dir + '/fragments'):
            os.system('rm -r {}/fragments'.format(hmm_dir))
        os.system('mkdir -p {}/fragments'.format(hmm_dir))

        # [NEW] save each single sequence to a fasta
        this_hmm = index_to_hmm[i].hmm_model_path
        for taxon in names:
            frag_path = '{}/fragments/{}.fasta'.format(hmm_dir, taxon)
            frag = unaligned.sub_alignment([taxon])
            frag.write(frag_path, 'FASTA')
        
            # hmmalign
            hmmalign_result_path = '{}/hmmalign.results.{}.{}.out'.format(
                    hmm_dir, taxon, i)
            cmd = '{} -o {} {} {}'.format(Configs.hmmalign_path,
                    hmmalign_result_path, this_hmm, frag_path)
            if not (os.path.exists(hmmalign_result_path) and os.path.getsize(
                hmmalign_result_path) > 0):
                os.system(cmd)
        
            # Extended alignment
            ap_aln = ExtendedAlignment(frag.get_sequence_names())
            ap_aln.build_extended_alignment(index_to_hmm[i].alignment,
                    [hmmalign_result_path], True)
            #ap_aln.relabel_original_columns(remaining_cols)
            for key in ap_aln.keys():
                ap_aln[key] = ap_aln[key].upper()
            # save extended alignment to backbone directory
            this_bb_path = backbone_dir + '/{}_{}.fasta'.format(taxon, i)
            ap_aln.write(this_bb_path, 'FASTA')

            # add the weights of the bb to weights.txt
            if Configs.use_weight:
                real_this_bb_path = os.popen('realpath -s {}'.format(
                    this_bb_path)).read().split('\n')[0]
                weights_file.write('{},{}\n'.format(
                        real_this_bb_path, Weights.weights_map[taxon][i]))
    weights_file.close()

'''
The core function of GCM+eHMMs. Merge the HMMAlign results of the set of
queries to the input alignment.

The main differences between GCM+eHMMs and UPP are: 1) GCM+eHMMs uses adjusted
bitscore, which takes the number of queries in an HMM into consideration. 2)
GCM+eHMMs can utilize more than one HMM (while UPP uses the best HMM based on
bitscore) to align the queries; hence, more information is used.
'''
def alignSubQueries(index_to_hmm, lock, index):
    s11 = time.time()
    # add constraints
    constraints_dir = Configs.outdir + '/constraints/{}'.format(index)
    if not os.path.isdir(constraints_dir):
        os.system('mkdir -p {}'.format(constraints_dir))
    # 0-th constraint set comes from the input alignment
    os.system('cp {} {}/c0.fasta'.format(Configs.backbone_path,
        constraints_dir))

    unaligned_dir = Configs.outdir + '/data'
    unaligned_path = unaligned_dir + '/unaligned_frag_{}.txt'.format(index)
    unaligned = Alignment(); unaligned.read_file_object(unaligned_path)
    unaligned_names = unaligned.get_sequence_names()

    c_index = 1
    for name in unaligned_names:
        single_frag = unaligned.sub_alignment([name])
        single_frag.write('{}/c{}.fasta'.format(constraints_dir, c_index),
                'FASTA')
        c_index += 1
    time_obtain_constraints = time.time() - s11

    s12 = time.time()
    # backbone alignments for GCM
    bb_dir = Configs.outdir + '/backbone_alignments/{}'.format(index)
    if not os.path.isdir(bb_dir):
        os.system('mkdir -p {}'.format(bb_dir))
    # hmmsearch directory
    hmmsearch_dir = Configs.outdir + '/search_results/{}'.format(index)
    if not os.path.isdir(hmmsearch_dir):
        os.system('mkdir -p {}'.format(hmmsearch_dir))

    # get backbones with the information we have
    getBackbones(index_to_hmm, unaligned, hmmsearch_dir, bb_dir)
    time_obtain_backbones = time.time() - s12

    s13 = time.time()
    # run GCM (modified MAGUS which takes in weights) on the subset
    est_path = Configs.outdir + '/magus_result_{}.txt'.format(index)
    gcm_outdir = Configs.outdir + '/magus_outputs/{}'.format(index)
    if Configs.use_weight:
        weights_path = hmmsearch_dir + '/weights.txt'
        cmd = 'python3 {} -np {} \
                -d {} -s {} -b {} -o {} \
                -w {} -f {} --graphclustermethod {} \
                --graphtracemethod {} --graphtraceoptimize {}'.format(
                Configs.magus_path, 1,
                gcm_outdir, constraints_dir, bb_dir, est_path,
                weights_path, Configs.inflation_factor, 
                Configs.graphclustermethod,
                Configs.graphtracemethod, Configs.graphtraceoptimize)
    else:
        cmd = 'python3 {} -np {} \
                -d {} -s {} -b {} -o {} -f {} \
                --graphclustermethod {} \
                --graphtracemethod {} --graphtraceoptimize {}'.format(
                Configs.magus_path, 1,
                gcm_outdir, constraints_dir, bb_dir, est_path,
                Configs.inflation_factor,
                Configs.graphclustermethod, Configs.graphtracemethod,
                Configs.graphtraceoptimize)
    os.system(cmd)

    # remove temp folders
    if not Configs.keeptemp:
        os.system('rm -r {}'.format(hmmsearch_dir))
        os.system('rm -r {}'.format(bb_dir))
        os.system('rm -r {}'.format(constraints_dir))
        if not Configs.keepgcmtemp:
            os.system('rm -r {}'.format(gcm_outdir))
    time_gcm = time.time() - s13
    
    lock.acquire()
    try:
        Configs.runtime(' '.join(['(alignSubQueries,',
                'i={}) Time to obtain'.format(index),
                'constraints (s):', str(time_obtain_constraints)]))
        Configs.runtime(' '.join(['(alignSubQueries,',
                'i={}) Time to obtain backbones (s):'.format(index),
                str(time_obtain_backbones)]))
        Configs.runtime(' '.join(['(alignSubQueries,',
                'i={}) Time to run GCM and'.format(index),
                'clean temporary files (s):', str(time_gcm)]))
        Configs.debug("Command used: {}".format(cmd))
        Configs.warning('{} passed to main pipeline!'.format(est_path))
    finally:
        lock.release()
    
    return est_path