import os, time, math, gzip, tempfile
from collections import defaultdict
from witch_msa.configs import Configs, tqdm_styles
from witch_msa.helpers.alignment_tools import ExtendedAlignment, Alignment, read_fasta
from witch_msa.gcmm.task import getTasks, runTasks

#from multiprocessing import Lock, Pool
import concurrent.futures
import subprocess
from functools import partial
from tqdm import tqdm

'''
A class holding directory/alignment information for each HMM subset from the
UPP decomposition
'''
class HMMSubset(object):
    def __init__(self, path, index):
        self.alignment_dir = path
        self.index = index
        self.num_taxa = 0

        # get hmm build results from target directory
        cmd = "find {} -maxdepth 1 -name hmmbuild.input* -type f".format(path)
        #Configs.debug("Command used: {}".format(cmd))
        #self.alignment_path = os.path.realpath(
        #        os.popen(cmd).read().split('\n')[0])
        self.alignment_path = os.path.realpath(subprocess.run(cmd.split(),
                capture_output=True, text=True).stdout.split('\n')[0])

        # also get the hmm model path
        cmd = "find {} -maxdepth 1 -name hmmbuild.model.* -type f".format(path)
        #Configs.debug("Command used: {}".format(cmd))
        #self.hmm_model_path = os.path.realpath(
        #        os.popen(cmd).read().split('\n')[0])
        self.hmm_model_path = os.path.realpath(subprocess.run(cmd.split(),
                capture_output=True, text=True).stdout.split('\n')[0])

        _map = {}
        with open(self.hmm_model_path, 'r') as f:
            line = f.readline().strip()
            # Modified @ 2.15.2024 - Chengze Shen
            # instead of reading in all lines at once, read one by one until
            # _safe_trigger or encounter field "" NSEQ ""
            # maybe it can alleviate memory issue encountered with large
            # dataset
            _idx, _safe_trigger = 0, 20
            while line and _idx < _safe_trigger:
                key, val = line.split()[0], ''.join(line.split()[1:])
                if key == 'NSEQ':
                    self.num_taxa = int(val)
                    break
                line = f.readline().strip()
                _idx += 1
            if self.num_taxa == 0:
                raise ValueError(
                        'Cannot find field: NSEQ, from {}'.format(
                            path))

            #lines = f.read().split('\n')[:15]
            #for line in lines:
            #    key, val = line.split()[0], ''.join(line.split()[1:])
            #    if key == 'NSEQ':
            #        self.num_taxa = int(val)
            #        break

'''
Function to write an extended alignment output (for a query alignment job)
to local as intermediate outputs. These can be used as checkpoints to avoid
rerunning the jobs if the main task is interrupted.
Writing and reading to compressed format (gzip by default)
This needs to be as compact as possible, possibly only saving a few fields:
    1. the query taxon name
    2. the intermediate query alignment
'''
def writeOneCheckpointAlignment(path, query):
    # make sure we are dealing with ExtendedAlignment object
    if not isinstance(query, ExtendedAlignment):
        return
    # make sure only one taxon is present (hence, query alignment)
    if len(query) != 1:
        return
    taxon = list(query.keys())[0]; seq = query[taxon]
    line = "{}\t{}\n".format(taxon, seq)
    encoded = line.encode('utf-8')

    # every line corresponds to a query taxon and its seq, separated by a tab
    with gzip.open(path, 'ab') as f:
        f.write(encoded)

'''
Helper function to process a splice of lines (query taxa) when reading 
checkpoint alignments
'''
def readOneCheckpointAlignment(lines):
    subset_checkpoint_queries = []
    for line in lines:
        taxon = '\t'.join(line.split('\t')[:-1])
        seq = line.split('\t')[-1]
        query = ExtendedAlignment([])
        query[taxon] = seq; query._reset_col_names()
        insertion = -1; regular = 0
        for i in range(len(seq)):
            if seq[i].islower():
                query._col_labels[i] = insertion; insertion -= 1
            else:
                query._col_labels[i] = regular; regular += 1
        subset_checkpoint_queries.append(query)

    return subset_checkpoint_queries

'''
Function to read from existing checkpoint query alignments, if any exists
(e.g., ./checkpoint_alignments.txt)
Use multiple cores to simulatenously process lines
'''
def readCheckpointAlignments(path, pool, lock):
    Configs.log('Reading existing checkpoint query alignments: {}'.format(
        path))
    checkpoint_queries = {}
    s1 = time.time()

    # read from the file and decode using utf-8 (gzipped)
    raw_inputs = []
    with gzip.open(path, 'rb') as f:
        raw_inputs = f.read().decode('utf-8').split('\n')[:-1]
    num_inputs = len(raw_inputs)

    # split tasks evenly based on the number of cpus
    per_task = int(num_inputs / pool._max_workers)

    futures = []
    for i in range(0, num_inputs, per_task):
        end_ind = min(num_inputs, i + per_task)
        futures.append(pool.submit(readOneCheckpointAlignment,
            raw_inputs[i:end_ind]))
    for future in tqdm(
            concurrent.futures.as_completed(futures),
            total=num_inputs, **tqdm_styles):
        subset_checkpoint_queries = future.result()
        for query in subset_checkpoint_queries:
            checkpoint_queries[list(query.keys())[0]] = query

    runtime_read_checkpoint = time.time() - s1
    Configs.log('Finished reading {} existing checkpoint query alignments'.format(
        len(checkpoint_queries)))
    Configs.runtime(' '.join(['(readCheckpointAlignments) Time to read',
            'existing checkpoint query alignments (s):',
            str(runtime_read_checkpoint)]))
    return checkpoint_queries

'''
Function to create a local copy of the backbone alignment in upper-cases
'''
def writeTempBackbone(outdir, backbone_path):
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    tmp_backbone_path = os.path.join(outdir, 'backbone.aln.fasta')
    Configs.log('Creating a local copy of backbone alignment ' + \
            '(all letters to upper-cases) at: ' + \
            tmp_backbone_path)
    alignment = Alignment(); alignment.read_file_object(backbone_path)
    backbone_length = alignment.sequence_length()

    assert backbone_length != None, \
            'The input backbone {} is not aligned!'.format(backbone_path)

    for key in alignment.keys():
        alignment[key] = alignment[key].upper()
    alignment.write(tmp_backbone_path, 'FASTA')
    del alignment

    return tmp_backbone_path, backbone_length

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
    taxa_names = list(unaligned.keys())
    num_seq = len(taxa_names)

    Configs.log('Started splitting queries (N={}) to subsets...'.format(
        num_seq))
    sid_to_query_names = {}; sid_to_query_seqs = {}
    # also rename taxon name with illegal characters
    renamed_taxa = {}
    for i in range(0, num_seq):
        if '/' in taxa_names[i]:
            # replace the taxon name at the same position in taxa_names
            old_name = taxa_names[i]
            new_name = 'renamed_query_{}'.format(i)
            taxa_names[i] = new_name

            renamed_taxa[old_name] = new_name
            unaligned[new_name] = unaligned[old_name]
            unaligned.pop(old_name)

        seq = unaligned[taxa_names[i]]
        sid_to_query_names[i] = taxa_names[i]; sid_to_query_seqs[i] = seq
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
    subset_args = align_dirs

    mytasks = getTasks(subset_args)
    ret, _, _, _ = runTasks(func, pool, mytasks, len(subset_args),
            max_concurrent_jobs=Configs.max_concurrent_jobs)
    for d in ret:
        index_to_hmm.update(d)

    #all_dicts = list(pool.map(func, align_dirs))
    #for d in all_dicts:
    #    index_to_hmm.update(d)

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
    print('\nReading and ranking bit-scores from HMMSearch files...')

    # submit jobs to Pool
    func = partial(readHMMSearch, lock, total_num_models)
    all_subset_ranks, futures = [], []
    for value in index_to_hmm.values():
        futures.append(pool.submit(func, value))
    for future in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures), **tqdm_styles):
        all_subset_ranks.append(future.result())
    #all_subset_ranks = list(pool.map(func, index_to_hmm.values()))

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
Rank bit-scores
'''
def rankBitscores(index_to_hmm, renamed_taxa, lock, pool):
    s1 = time.time()
    ranked_bitscore = readAndRankBitscoreMP(index_to_hmm, renamed_taxa,
            lock, pool)
    time_rank_bitscore = time.time() - s1
    Configs.runtime(' '.join(['(rankBitscores) Time to rank bit-scores (s)',
        str(time_rank_bitscore)]))
    return ranked_bitscore

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
    #ranked_bitscore = readAndRankBitscoreMP(index_to_hmm, renamed_taxa,
    #                                        lock, pool)

    time_load_files = time.time() - s1
    Configs.runtime(' '.join(['(loadSubQueries) Time to split queries',
            '(s):', str(time_load_files)]))
    return num_seq, index_to_hmm, sid_to_query_names, \
            sid_to_query_seqs, renamed_taxa
