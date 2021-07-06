import os, sys, re, math
import numpy as np
import logging
from Bio import AlignIO
sys.path.append("/home/chengze5/tallis")
sys.path.append("/home/chengze5/tallis/softwares/sepp/sepp")
from alignment import *
from alignment_tools import * 
from collections import defaultdict

NUM_THREAD = 8
SUBSET_SIZE = 1 
TOP_FIT_THRESHOLD = 0.75
WEIGHT_THRESHOLD = 0.95
keeptemp = False 

# function to calculate the HMM weighting, given the bitscores and sizes
# of the HMMs (for a given query taxon)
# inputs: ensemble of HMMs H (with their bitscores and sizes)
# outputs: weights for HMMs H
#
def calculateWeights(taxon, indexes, bitscores, sizes):
    #logging.debug('working with: {}'.format(taxon))
    weights = {}
    
    temp = []
    assert len(indexes) == len(bitscores) == len(sizes)
    for i in range(len(bitscores)):
        score_i, size_i = bitscores[i], sizes[i]
        exponents = np.array(bitscores) - score_i \
                + np.log2(np.array(sizes) / size_i)
        denominator = np.sum(np.power(2, exponents))
        temp.append(1. / denominator)
        weights[indexes[i]] = 1. / denominator
    return weights, temp

def getHMMBuildAlignment(dir, index):
    # get hmm build results from target directory
    cmd = "find {} -name hmmbuild.input* -type f".format(dir)
    logging.debug("Command used: %s", cmd)
    hmmbuild_align_path = os.popen(cmd).read().split('\n')[0]
    hmmbuild_align = Alignment(); hmmbuild_align.read_file_object(
            hmmbuild_align_path)
    logging.debug("Read %d sequences from HMMBuild %d",
            hmmbuild_align.get_num_taxa(), index)

    # also get the HMM model path
    cmd = "find {} -name hmmbuild.model.* -type f".format(dir)
    hmmbuild_model_path = os.popen(cmd).read().split('\n')[0]
    return hmmbuild_align, hmmbuild_align_path, hmmbuild_model_path

def getAlignmentSubsets(path):
    cmd = "find {} -name A_0_* -type d".format(path)
    logging.debug("Command used: %s", cmd)
    align_dirs = os.popen(cmd).read().split('\n')[:-1]
    logging.info("Number of alignment subsets: %d", len(align_dirs))

    # create an alignment for each align_dir
    align_subsets, index_to_alignment, index_to_model, index_to_dir \
            = {}, {}, {}, {}
    for align_dir in align_dirs:
        index = int(align_dir.split('/')[-1].split('_')[2])
        alignment, align_path, model_path = getHMMBuildAlignment(align_dir,
                index)
        index_to_alignment[index] = align_path
        index_to_model[index] = model_path
        index_to_dir[index] = align_dir
        align_subsets[index] = alignment

    return align_subsets, index_to_alignment, index_to_model, index_to_dir

def getFragments(paths):
    names = []
    for path in paths:
        a = Alignment(); a.read_file_object(path)
        names += a.get_sequence_names()
    return names

def getBackbones(k, index_to_alignment, index_to_model, index_to_dir, ranks,
        unaligned_path, unaligned_frag, tempdir, outdir, weights, weights_map,
        weights_path, strategy='top_k'):
    # return is a map from HMM index to a list of fragments that have high
    # bit scores at the HMM.
    ret = defaultdict(list)

    # for each taxon and its scores on HMMs, sort it in decreasing order
    for taxon, sorted_weights in weights.items():
        if taxon not in unaligned_frag:
            continue

        # use different strategy to consider which HMMs to use for 
        # a query
        if strategy == 'top_k':
            # a strategy to select the top k HMMs, k is user-defined
            for item in sorted_weights[:k]:
                # append taxon to corresponding HMM i (defined by the index)
                ret[item[0]].append(taxon)
        else:
            print("Strategy {} not supported!".format(strategy))

    # now we can split the fragments into chunks, each chunk corresponds to
    # all fragments added to an HMM i
    # For each HMM, run HMMAlign on its fragment chunk to get a backbone
    # alignment
    weights_file = open(weights_path, 'w')

    for index, names in ret.items():
        logging.warning("Generating fragment chunks/alignment for HMM {}".format(
            index))
        hmm_dir = tempdir + '/A_0_{}'.format(index)
        if not os.path.isdir(hmm_dir):
            os.system('mkdir -p {}'.format(hmm_dir))
        if os.path.isdir(hmm_dir + '/fragments'):
            os.system('rm -r {}/fragments'.format(hmm_dir))
        os.system('mkdir -p {}/fragments'.format(hmm_dir))
        if os.path.isdir(hmm_dir + '/hmmalign'):
            os.system('rm -r {}/hmmalign'.format(hmm_dir))
        os.system('mkdir -p {}/hmmalign'.format(hmm_dir))

        # [NEW] save each single sequence to a fasta
        this_hmm = index_to_model[index]
        for taxon in names:
            frag_path = '{}/fragments/fragment_{}.fasta'.format(hmm_dir,
                    index)
            frag = unaligned_frag.sub_alignment([taxon])
            frag.write(frag_path, 'FASTA')
        
            # hmmalign
            hmmalign_result_path = '{}/hmmalign/hmmalign.results.{}.out'.format(
                    hmm_dir, index)
            cmd = 'hmmalign -o {} {} {}'.format(hmmalign_result_path,
                    this_hmm, frag_path)
            if not (os.path.exists(hmmalign_result_path) and os.path.getsize(
                hmmalign_result_path) > 0):
                os.system(cmd)
        
            # Extended alignment
            ap_aln = ExtendedAlignment(frag.get_sequence_names())
            ap_aln.build_extended_alignment(index_to_alignment[index],
                    [hmmalign_result_path], True)
            #ap_aln.relabel_original_columns(remaining_cols)
            for key in ap_aln.keys():
                ap_aln[key] = ap_aln[key].upper()
            # save extended alignment
            this_bb_path = outdir + '/bb{}.fasta'.format(index)
            ap_aln.write(this_bb_path, 'FASTA')

            # add the weights of the bb to weights.txt
            real_this_bb_path = os.popen('realpath -s {}'.format(this_bb_path))\
                                    .read().split('\n')[0]
            weights_file.write('{},{}\n'.format(
                    real_this_bb_path, weights_map[taxon][index]))
    weights_file.close()
           

# function to split unaligned sequences to subsets of size (at most) 25
# so that GCM can run for each subset
def getSubQueries(unaligned, outdir, num, num_subset):
    if not os.path.isdir(outdir):
        os.system('mkdir -p {}'.format(outdir))

    # get names
    frag_names = unaligned.get_sequence_names()

    # get #[num] names from top of frag_names
    ind = 0
    for i in range(0, num_subset):
        start_ind, end_ind = i*SUBSET_SIZE, min(num, (i+1)*SUBSET_SIZE)
        chosen_frags = frag_names[start_ind:end_ind]
        # get subsets of sequences
        sub_unaligned = unaligned.sub_alignment(chosen_frags)
        # save to outdir
        sub_unaligned.write(outdir + '/unaligned_frag_{}.txt'.format(i), 'FASTA')

# function to align a single subset of queries using GCM and HMMs
def alignSubQueries(align_subsets, index_to_alignment, index_to_model,
        index_to_dir, ranks,
        backbone_path, outdir, unaligned_dir, k, index, num_subset,
        strategy, weights, weights_map):
    # add constraints
    constraints_dir = outdir + '/constraints/{}'.format(index)
    if not os.path.isdir(constraints_dir):
        os.system('mkdir -p {}'.format(constraints_dir))
    #backbone = os.popen('find {} -name pastajob.marker001.* -type f'.format(
    #    pasta_dir)).read().split('\n')[:-1]
    #assert len(backbone) > 0; backbone = backbone[0]
    #os.system('cp {} {}/c0.fasta'.format(backbone, constraints_dir))
    os.system('cp {} {}/c0.fasta'.format(backbone_path, constraints_dir))

    unaligned_path = unaligned_dir + '/unaligned_frag_{}.txt'.format(index)
    unaligned = Alignment(); unaligned.read_file_object(unaligned_path)
    unaligned_names = unaligned.get_sequence_names()
    c_index = 1
    for name in unaligned_names:
        single_frag = unaligned.sub_alignment([name])
        single_frag.write('{}/c{}.fasta'.format(constraints_dir, c_index),
                'FASTA')
        c_index += 1
    
    # backbone alignments for GCM
    bb_dir = outdir + '/backbone_alignments/{}'.format(index)
    if not os.path.isdir(bb_dir):
        os.system('mkdir -p {}'.format(bb_dir))
    # hmmsearch directory
    hmmsearch_dir = outdir + '/search_results/{}'.format(index)
    if not os.path.isdir(hmmsearch_dir):
        os.system('mkdir -p {}'.format(hmmsearch_dir))

    # get backbones with the information we have
    weights_path = hmmsearch_dir + '/weights.txt'
    getBackbones(k, index_to_alignment, index_to_model, index_to_dir, ranks,
            unaligned_path, unaligned, hmmsearch_dir, bb_dir,
            weights, weights_map, weights_path, strategy=strategy)
    
    # run GCM (modified MAGUS which takes in weights) on the subset
    #magusbin = '/home/chengze5/tallis/softwares/MAGUS/magus.py'
    magusbin = '/home/chengze5/tallis/softwares/modified_MAGUS/magus.py'
    est_path = outdir + '/magus_result_{}.txt'.format(index)
    cmd = 'python3 {} -np {} -d {}/magus_outputs/{} -s {} -b {} -o {} -w {}'.format(
            magusbin, NUM_THREAD, outdir, index, constraints_dir, bb_dir,
            est_path, weights_path)
    os.system(cmd)

    # remove temp folders
    if not keeptemp:
        os.system('rm -r {}'.format(hmmsearch_dir))
        os.system('rm -r {}/magus_outputs/{}'.format(outdir, index))
        os.system('rm -r {}/backbone_alignments/{}'.format(outdir, index))
        os.system('rm -r {}/constraints/{}'.format(outdir, index))

def main():
    if len(sys.argv) < 6:
        print("Parameters required: [hmm path] [backbone path] [outdir] [refpath] [k]")
        exit()
    global keeptemp, TOP_FIT_THRESHOLD
    path = sys.argv[1]
    backbone_path = sys.argv[2]
    outdir = sys.argv[3]
    ref_path = sys.argv[4]
    # check if it is similar_score
    if float(sys.argv[5]) < 1:
        strategy = 'similar_score'
        TOP_FIT_THRESHOLD = float(sys.argv[5])
        k_hmms = 1
    else:
        strategy = 'top_k'
        k_hmms = int(sys.argv[5])
    keeptemp = False 
    logging.basicConfig(level=logging.WARNING)
    
    # alignment logger
    aln_logger = logging.getLogger('alignment')
    aln_logger.setLevel(logging.WARNING)

    if len(sys.argv) > 6:
        global NUM_THREAD
        NUM_THREAD = int(sys.argv[6])

    if not os.path.isdir(path):
        print("Path {} does not exist.".format(sys.argv[1]))
        exit()
    if not os.path.isfile(ref_path):
        print("Reference aln does not exist.")
        exit()

    # get ref and unaligned
    #ref = Alignment(); ref.read_file_object(ref_path)
    unaligned_path = '/'.join(ref_path.split('/')[:-1]) + '/unaligned_frag.txt'
    unaligned = Alignment(); unaligned.read_file_object(unaligned_path)

    # 1) get all sub-queries, write to [outdir]/data
    num_unaligned = len(unaligned)
    num_subset = max(1, math.ceil(num_unaligned / SUBSET_SIZE))
    data_dir = outdir + '/data'
    getSubQueries(unaligned, data_dir, len(unaligned), num_subset)
 
    # 1.2) read in all HMMSearch results (from UPP)
    align_subsets, index_to_alignment, index_to_model, index_to_dir \
            = getAlignmentSubsets(path)
    ranks = defaultdict(list)
    total_num_models, counter = len(index_to_dir), 0
    for index, hmmdir in index_to_dir.items():
        counter += 1
        logging.warning("reading HMMSearch result {}/{}".format(
            counter, total_num_models))

        cmd = 'find {} -name hmmsearch.results.* -type f'.format(hmmdir)
        hmmsearch_paths = os.popen(cmd).read().split('\n')[:-1]
        for each_path in hmmsearch_paths:
            with open(each_path, 'r') as f:
                temp_map = eval(f.read())
                # mapping of taxon->scores for this index
                for taxon, scores in temp_map.items():
                    ranks[taxon].append((index, scores[1]))
    # write to local
    with open(outdir + '/ranked_scores.txt', 'w') as f:
        for taxon, scores in ranks.items():
            sorted_scores = sorted(scores, key = lambda x: x[1], reverse=True)
            f.write(taxon + ':' + ';'.join(
                [str(z) for z in sorted_scores]) + '\n')
    logging.warning("Finished writing ranked scores to local")

    # 1.3) obtain the weights of each query
    # - get sizes of each HMM
    cmd = 'find {} -name hmmbuild.input.* -type f'.format(path)
    hmmbuild_inputs = os.popen(cmd).read().split('\n')[:-1]
    all_sizes = {}
    for hmmbuild_input_path in hmmbuild_inputs:
        index = int(hmmbuild_input_path.split('/')[-2].split('_')[-1])
        a = Alignment(); a.read_file_object(hmmbuild_input_path)
        all_sizes[index] = a.get_num_taxa()
        del a

    # iterate through each query taxon
    weights, weights_map = {}, {}
    for taxon, sorted_scores in ranks.items():
        indexes = [x[0] for x in sorted_scores]
        bitscores = [x[1] for x in sorted_scores]
        sizes = [all_sizes[x] for x in indexes]
        this_weights_map,_ = calculateWeights(taxon, indexes, bitscores, sizes)
        weights[taxon] = sorted([(ind, w) for ind, w in this_weights_map.items()],
                key = lambda x: x[1], reverse=True)
        weights_map[taxon] = this_weights_map

    # write weights to local
    with open('{}/all_weights.txt'.format(outdir), 'w') as f:
        for taxon, weight in weights.items():
            w = [str(x) for x in weight]
            f.write(taxon + ':' + ';'.join(w) + '\n')

    # 2) now deal with each subset alignment individually, and merge them
    # back at the end
    for i in range(num_subset):
        alignSubQueries(align_subsets, index_to_alignment,
                index_to_model, index_to_dir, ranks,
                backbone_path, outdir, data_dir, k_hmms, i, num_subset,
                strategy, weights, weights_map)
    
    # 3) merge all magus results together to form a merged fasta
    logging.info("merging %d MAGUS results...", num_subset)
    merged_path = outdir + '/merged.fasta'
    cmd = 'python3 merger.py {} {}'.format(outdir, merged_path)
    os.system(cmd)

    # 3.2) remove intermediate results
    if not keeptemp:
        os.system('rm {}/magus_result_*'.format(outdir))
    if not keeptemp and os.path.isdir('{}/backbone_alignments'.format(outdir)):
        os.system('rm -r {}/backbone_alignments'.format(outdir))
    if not keeptemp and os.path.isdir('{}/constraints'.format(outdir)):
        os.system('rm -r {}/constraints'.format(outdir))
    if not keeptemp and os.path.isdir('{}/magus_outputs'.format(outdir)):
        os.system('rm -r {}/magus_outputs'.format(outdir))

    # 4) run FastSP to get the alignment metrics
    fastspbin = '/home/chengze5/tallis/softwares/FastSP/FastSP.jar'
    cmd = 'java -Xmx4096m -jar {} -ml -mlr -r {} -e {} -o {}/fastsp.out'.format(
            fastspbin, ref_path, merged_path, outdir)
    os.system(cmd)

if __name__ == "__main__":
    main()
