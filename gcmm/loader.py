import os
import time
import math
import tempfile
from configs import Configs
from collections import defaultdict
from helpers.alignment_tools import Alignment, read_fasta

#from multiprocessing import Lock, Pool
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
        cmd = "find {} -maxdepth 1 -name hmmbuild.input* -type f".format(path)
        #Configs.debug("Command used: {}".format(cmd))
        self.alignment_path = os.path.realpath(
                os.popen(cmd).read().split('\n')[0])
        #self.alignment = Alignment()
        #self.alignment.read_file_object(self.alignment_path)

        # also get the hmm model path
        cmd = "find {} -maxdepth 1 -name hmmbuild.model.* -type f".format(path)
        #Configs.debug("Command used: {}".format(cmd))
        self.hmm_model_path = os.path.realpath(
                os.popen(cmd).read().split('\n')[0])

        self.num_taxa = 0
        _map = {}
        with open(self.hmm_model_path, 'r') as f:
            lines = f.read().split('\n')[:15]
            for line in lines:
                key, val = line.split()[0], ''.join(line.split()[1:])
                if key == 'NSEQ':
                    self.num_taxa = int(val)
                    break

'''
Initialize global lock for writing to log file
'''
#def init_lock_loader(l):
#    global lock
#    lock = l

'''
Function to write a single set of queries to local, given its index, etc
'''
def writeOneQuerySet(frag_names, unaligned, outdir, args):
    index, start_ind, end_ind = args
    chosen_frags = frag_names[start_ind:end_ind]
    sub_unaligned = unaligned.sub_alignment(chosen_frags)
    sub_unaligned.write(outdir + '/unaligned_frag_{}.txt'.format(index),
            'FASTA')

'''
Split and write sub-queries to local
'''
def writeSubQueries(unaligned, outdir, pool):
    #if not os.path.isdir(outdir):
    #    os.system('mkdir -p {}'.format(outdir))
    frag_names = unaligned.get_sequence_names()

    # rename taxon name with illegal characters
    renamed_taxa = {}
    for taxon in frag_names:
        if '/' in taxon:
            taxon_name = tempfile.mktemp().split('/')[-1]
            renamed_taxa[taxon] = taxon_name
            unaligned[taxon_name] = unaligned[taxon]
            unaligned.pop(taxon)
    frag_names = list(unaligned.keys())
    num_seq = len(frag_names)

    Configs.log('Started splitting queries (N={}) to subsets...'.format(
        num_seq))
    sid_to_query_names = {}; sid_to_query_seqs = {}
    for i in range(0, num_seq):
        taxon = frag_names[i]; seq = unaligned[taxon]
        #subaln = unaligned.sub_alignment([frag_names[i]])
        sid_to_query_names[i] = taxon; sid_to_query_seqs[i] = seq
    Configs.log('Finished splitting queries in memory.')
    #Configs.log('Index to query map: {}'.format(str(sid_to_query_names)))
    #args = []
    #for i in range(0, num_subset):
    #    start_ind, end_ind = i * Configs.subset_size, \
    #            min(num_seq, (i + 1) * Configs.subset_size)
    #    args.append((i, start_ind, end_ind))
    #func = partial(writeOneQuerySet, frag_names, unaligned, outdir)
    #pool.map(func, args)
    #Configs.log('Finished splitting queries to local.')

    if len(renamed_taxa) > 0:
        Configs.log('The following taxa are renamed '
                '(names will be reverted in the output): '
                '{}'.format(renamed_taxa))
    return num_seq, sid_to_query_names, sid_to_query_seqs, renamed_taxa

'''
Helper function to load one alignment subset
'''
def getOneAlignmentSubset(lock, align_dir):
    this_dict = {}
    index = int(align_dir.split('/')[-1].split('_')[2])
    this_dict[index] = HMMSubset(align_dir, index) 
    lock.acquire()
    try:
        Configs.debug('Finished dealing with subset #{} from {}'.format(
            index, align_dir))
    finally:
        lock.release()

    return this_dict

'''
Load in UPP decomposition output subsets
'''
def getAlignmentSubsets(path, lock, pool):
    cmd = "find {} -name A_0_* -type d".format(path)
    #Configs.debug("Command used: {}".format(cmd))
    align_dirs = os.popen(cmd).read().split('\n')[:-1]
    Configs.log("Number of alignment subsets: {}".format(len(align_dirs)))

    # create an AlignmentSubset object for each align_dir
    index_to_hmm = {}
    func = partial(getOneAlignmentSubset, lock)
    all_dicts = list(pool.map(func, align_dirs))
    for d in all_dicts:
        index_to_hmm.update(d)
    #for align_dir in align_dirs:
    #    index = int(align_dir.split('/')[-1].split('_')[2])
    #    index_to_hmm[index] = HMMSubset(align_dir, index) 
    return index_to_hmm

'''
Read HMMSearch results given an HMMSubset object
'''
def readHMMSearch(lock, total_num_models, subset):
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
def readAndRankBitscoreMP(index_to_hmm, renamed_taxa, lock, pool):
    ranks = defaultdict(list)
    total_num_models = len(index_to_hmm)
    Configs.log('Reading and ranking bit-scores from HMMSearch files')

    # submit jobs to Pool
    func = partial(readHMMSearch, lock, total_num_models)
    all_subset_ranks = list(pool.map(func, index_to_hmm.values()))

    # merge all MP results
    Configs.log('Ranking bit-scores')
    for subset_ranks in all_subset_ranks:
        for taxon, ind_score_pairs in subset_ranks.items():
            for pair in ind_score_pairs:
                ranks[taxon].append(pair)

    # sort the bitscores and write to local
    ranked_bitscores = defaultdict(list)
    for taxon, scores in ranks.items():
        sorted_scores = sorted(scores, key = lambda x: x[1], reverse=True)
        taxon_name = taxon
        if taxon in renamed_taxa:
            taxon_name = renamed_taxa[taxon]
        ranked_bitscores[taxon_name] = sorted_scores
    Configs.log('Finished ranking bit-scores')
    return ranked_bitscores

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
def loadSubQueries(lock, pool):
    s1 = time.time()
    unaligned = Alignment()
    unaligned.read_file_object(Configs.query_path)
    # make sure unaligned sequences are read in upper-case letters
    for key in unaligned.keys():
        unaligned[key] = unaligned[key].upper()

    # 1) get all sub-queries, write to [outdir]/data
    data_dir = Configs.outdir + '/data'
    num_seq, sid_to_query_names, sid_to_query_seqs, \
            renamed_taxa = writeSubQueries(unaligned, data_dir, pool)
 
    # 1.2) read in all HMMSearch results (from UPP)
    index_to_hmm = getAlignmentSubsets(Configs.hmmdir, lock, pool)

    # 1.3) read and rank bitscores
    #ranked_bitscore = readAndRankBitscore(index_to_hmm)
    ranked_bitscore = readAndRankBitscoreMP(index_to_hmm, renamed_taxa,
                                            lock, pool)
    time_load_files = time.time() - s1
    Configs.runtime('Time to split queries and rank bit-scores (s): {}'.format(
        time_load_files))
    return num_seq, index_to_hmm, ranked_bitscore, sid_to_query_names, \
            sid_to_query_seqs, renamed_taxa
