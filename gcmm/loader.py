import os
import time
import math
from configs import Configs
from collections import defaultdict
from helpers.alignment_tools import Alignment, read_fasta
from multiprocessing import Lock, Pool
from functools import partial

'''
A class holding directory/alignment information for each HMM subset from the
UPP decomposition
'''
class HMMSubset(object):
    def __init__(self, path, index):
        self.alignment_dir = path
        self.index = index

        # get hmm build results from target directory
        cmd = "find {} -name hmmbuild.input* -type f".format(path)
        #Configs.debug("Command used: {}".format(cmd))
        self.alignment_path = os.path.abspath(
                os.popen(cmd).read().split('\n')[0])
        self.alignment = Alignment()
        self.alignment.read_file_object(self.alignment_path)

        # also get the hmm model path
        cmd = "find {} -name hmmbuild.model.* -type f".format(path)
        #Configs.debug("Command used: {}".format(cmd))
        self.hmm_model_path = os.popen(cmd).read().split('\n')[0]

'''
Initialize global lock for writing to log file
'''
def init_lock_loader(l):
    global lock
    lock = l

'''
Function to write a single set of queries to local, given its index, etc
'''
def writeOneQuerySet(frag_names, unaligned, outdir, index, start_ind, end_ind):
    chosen_frags = frag_names[start_ind:end_ind]
    sub_unaligned = unaligned.sub_alignment(chosen_frags)
    sub_unaligned.write(outdir + '/unaligned_frag_{}.txt'.format(index),
            'FASTA')

'''
Split and write sub-queries to local
'''
def writeSubQueries(unaligned, outdir, num_seq, num_subset):
    if not os.path.isdir(outdir):
        os.system('mkdir -p {}'.format(outdir))
    frag_names = unaligned.get_sequence_names()

    args = []
    for i in range(0, num_subset):
        start_ind, end_ind = i * Configs.subset_size, \
                min(num_seq, (i + 1) * Configs.subset_size)
        args.append((i, start_ind, end_ind))
    pool = Pool(Configs.num_cpus)
    func = partial(writeOneQuerySet, frag_names, unaligned, outdir)
    pool.starmap(func, args)
    pool.close()
    pool.join()

    ## get #[num] names from top of frag_names
    #ind = 0
    #for i in range(0, num_subset):
    #    start_ind, end_ind = i * Configs.subset_size, \
    #            min(num_seq, (i + 1) * Configs.subset_size)
    #    # get subsets of sequences
    #    chosen_frags = frag_names[start_ind:end_ind]
    #    sub_unaligned = unaligned.sub_alignment(chosen_frags)

    #    # save to outdir
    #    sub_unaligned.write(outdir + '/unaligned_frag_{}.txt'.format(i),
    #            'FASTA')

'''
Load in UPP decomposition output subsets
'''
def getAlignmentSubsets(path):
    cmd = "find {} -name A_0_* -type d".format(path)
    #Configs.debug("Command used: {}".format(cmd))
    align_dirs = os.popen(cmd).read().split('\n')[:-1]
    Configs.log("Number of alignment subsets: {}".format(len(align_dirs)))

    # create an AlignmentSubset object for each align_dir
    index_to_hmm = {}
    for align_dir in align_dirs:
        index = int(align_dir.split('/')[-1].split('_')[2])
        index_to_hmm[index] = HMMSubset(align_dir, index) 
    return index_to_hmm

'''
Read HMMSearch results given an HMMSubset object
'''
def readHMMSearch(total_num_models, subset):
    global lock
    lock.acquire()
    try:
        Configs.log('Reading HMMSearch result {}/{}'.format(subset.index,
            total_num_models))
    finally:
        lock.release()

    cmd = 'find {} -name hmmsearch.results.* -type f'.format(
            subset.alignment_dir)
    hmmsearch_paths = os.popen(cmd).read().split('\n')[:-1]
    subset_ranks = defaultdict(list)
    for each_path in hmmsearch_paths:
        with open(each_path, 'r') as f:
            temp_map = eval(f.read())
            for taxon, scores in temp_map.items():
                subset_ranks[taxon].append((subset.index, scores[1]))
    return subset_ranks

'''
MP version of reading in and ranking bitscores from HMM subsets (from UPP)
'''
def readAndRankBitscoreMP(index_to_hmm):
    ranks = defaultdict(list)
    total_num_models = len(index_to_hmm)
    l = Lock()

    # submit jobs to Pool
    pool = Pool(Configs.num_cpus, initializer=init_lock_loader, initargs=(l,))
    func = partial(readHMMSearch, total_num_models)
    all_subset_ranks = pool.map(func, index_to_hmm.values())
    pool.close()
    pool.join()

    # merge all MP results
    for subset_ranks in all_subset_ranks:
        for taxon, ind_score_pairs in subset_ranks.items():
            for pair in ind_score_pairs:
                ranks[taxon].append(pair)

    # write to local
    ranked = defaultdict(list)
    with open(Configs.outdir + '/ranked_scores.txt', 'w') as f:
        for taxon, scores in ranks.items():
            sorted_scores = sorted(scores, key = lambda x: x[1], reverse=True)
            ranked[taxon] = sorted_scores
            f.write(taxon + ':' + ';'.join(
                [str(z) for z in sorted_scores]) + '\n')
    Configs.warning("Finished writing ranked scores to local!")
    return ranked

'''
Read in and rank bitscores from UPP decomposition
'''
#def readAndRankBitscore(index_to_hmm):
#    ranks = defaultdict(list)
#    total_num_models, counter = len(index_to_hmm), 0
#    for index, subset in index_to_hmm.items():
#        counter += 1
#        Configs.log("Reading HMMSearch result {}/{}".format(counter,
#            total_num_models))
#
#        cmd = 'find {} -name hmmsearch.results.* -type f'.format(
#            subset.alignment_dir)
#        #Configs.debug("Command used: {}".format(cmd))
#        hmmsearch_paths = os.popen(cmd).read().split('\n')[:-1]
#        for each_path in hmmsearch_paths:
#            with open(each_path, 'r') as f:
#                temp_map = eval(f.read())
#                # mapping of taxon->scores for this index
#                for taxon, scores in temp_map.items():
#                    ranks[taxon].append((index, scores[1]))
#    # write to local
#    ranked = defaultdict(list)
#    with open(Configs.outdir + '/ranked_scores.txt', 'w') as f:
#        for taxon, scores in ranks.items():
#            sorted_scores = sorted(scores, key = lambda x: x[1], reverse=True)
#            ranked[taxon] = sorted_scores
#            f.write(taxon + ':' + ';'.join(
#                [str(z) for z in sorted_scores]) + '\n')
#    Configs.warning("Finished writing ranked scores to local!")
#    return ranked

'''
Split query sequences into batches of a defined size
'''
def loadSubQueries(unaligned):
    s1 = time.time()
    #unaligned = Alignment()
    #unaligned.read_file_object(Configs.query_path)
    num_subset = max(1, math.ceil(len(unaligned) / Configs.subset_size))

    # 1) get all sub-queries, write to [outdir]/data
    data_dir = Configs.outdir + '/data'
    writeSubQueries(unaligned, data_dir, len(unaligned), num_subset)
 
    # 1.2) read in all HMMSearch results (from UPP)
    index_to_hmm = getAlignmentSubsets(Configs.hmmdir)

    # 1.3) read and rank bitscores
    #ranked_bitscore = readAndRankBitscore(index_to_hmm)
    ranked_bitscore = readAndRankBitscoreMP(index_to_hmm)
    time_load_files = time.time() - s1
    Configs.runtime('Time to load files and split queries (s): {}'.format(
        time_load_files))
    return index_to_hmm, ranked_bitscore
