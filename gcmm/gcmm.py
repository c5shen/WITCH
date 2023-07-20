'''
Created on 1.19.2022 by Chengze Shen

Main pipeline of WITCH
'''

import os, sys, math, psutil, shutil 
from configs import * 
from gcmm.algorithm import DecompositionAlgorithm, SearchAlgorithm
from gcmm.loader import loadSubQueries, writeTempBackbone
from gcmm.weighting import writeWeights, writeBitscores, writeWeightsToLocal
from gcmm.aligner import alignSubQueries, alignSubQueriesNew
from gcmm.backbone import BackboneJob
from gcmm.merger import mergeAlignmentsCollapsed

# experimental, customized ProcessPoolExecutor
from gcmm import WITCHProcessPoolExecutor

from helpers.alignment_tools import Alignment
from helpers.general_tools import memoryUsage

import multiprocessing as mp
import concurrent.futures
from multiprocessing import Lock, Queue, Manager#, Pool
from concurrent.futures.process import ProcessPoolExecutor
from functools import partial
from tqdm import tqdm

# max system recursion limit hard encoding to a large number
# a temp fix for dendropy tree recursion issues
sys.setrecursionlimit(10000)

'''
Delete all unnecessary intermediate files
'''
def clearTempFiles():
    # make an empty dir for rsync removal
    blank_dir = os.path.join(Configs.outdir, 'blank')
    if not os.path.isdir(blank_dir):
        os.makedirs(blank_dir)
    #if not Configs.keepsubalignment:
    if os.path.isdir('{}/temp'.format(Configs.outdir)):
        #os.system('rsync -a --delete {}/ {}/temp/'.format(blank_dir,
        #    Configs.outdir))
        os.system('rm -rf {}/temp'.format(Configs.outdir))

    dirs_to_remove = ['tree_decomp',
            'backbone_alignments', 'constraints',
            'search_results', 'data', 'weights', 'bitscores']
    if Configs.keep_decomposition:
        dirs_to_remove = dirs_to_remove[1:]

    for _d in dirs_to_remove:
        if os.path.isdir('{}/{}'.format(Configs.outdir, _d)):
            os.system('rsync -a --delete {}/ {}/{}/'.format(blank_dir,
                Configs.outdir, _d))
            os.system('rmdir {}/{}'.format(Configs.outdir, _d))

    if (not Configs.keepgcmtemp) \
            and os.path.isdir('{}/magus_outputs'.format(Configs.outdir)):
        os.system('rsync -a --delete {}/ {}/magus_outputs/'.format(blank_dir,
            Configs.outdir))
        os.system('rmdir {}/magus_outputs'.format(Configs.outdir))

    if os.path.isdir(blank_dir):
        os.system('rmdir {}'.format(blank_dir))

'''
Init function for a queue and get configurations for each worker
'''
def initiate_pool(args):
    buildConfigs(args)

'''
6.21.2023 - Additionally, init pool for query alignment with two dicts
          - used later
'''
def initiate_pool_query_alignment(q, args,
        subset_to_retained_columns,
        subset_to_nongaps_per_column):
    buildConfigs(args)
    alignSubQueries.q = q
    alignSubQueriesNew.subset_to_retained_columns = subset_to_retained_columns
    alignSubQueriesNew.subset_to_nongaps_per_column = \
            subset_to_nongaps_per_column

'''
Main process for WITCH 
'''
def mainAlignmentProcess(args):
    m = Manager()
    lock = m.Lock()
    #l = Lock()
    q = Queue()

    # initialize the main pool at the start so that it's not memory
    # intensive
    Configs.warning('Initializing ProcessorPoolExecutor instance...')
    pool = ProcessPoolExecutor(Configs.num_cpus,
            initializer=initiate_pool, initargs=(args,))
            #mp_context=mp.get_context('spawn'),

    # if not user provided, default to <outdir>/tree_comp/root
    if not Configs.hmmdir:
        Configs.hmmdir = Configs.outdir + '/tree_decomp/root'

    # 0) obtain the backbone alignment/tree and eHMMs
    # If no UPP eHMM directory detected, decompose from the backbone
    if not os.path.exists(Configs.hmmdir):
        # if both backbone alignment/tree present
        if Configs.backbone_path and Configs.backbone_tree_path:
            pass
        else:
            # if missing backbone alignment, first perform an UPP-like
            # sequence split into backbone/query sets (i.e., randomly select
            # up to 1,000 sequences from 25% of the median length to be
            # backbone sequences)
            # Then, align the backbone sequences using MAGUS/PASTA
            # Finally, generate a FastTree2 backbone tree

            # jobs will only run if the corresponding paths are missing
            print('\nPerforming backbone alignment and/or tree estimation...')
            bb_job = BackboneJob()
            bb_job.setup()

            Configs.backbone_path, Configs.query_path = \
                    bb_job.run_alignment()
            Configs.backbone_tree_path = bb_job.run_tree()

        # after obtaining backbone alignment/tree, perform decomposition
        # and HMMSearches
        print('\nDecomposing the backbone tree...')
        decomp = DecompositionAlgorithm(
                Configs.backbone_path, Configs.backbone_tree_path,
                Configs.alignment_size)
        hmmbuild_paths, subset_to_retained_columns, subset_to_nongaps_per_column = \
                decomp.decomposition(lock, pool)
        print('\nPerforming all-against-all HMMSearches ' \
                'between the backbone and queries...')
        search = SearchAlgorithm(hmmbuild_paths)
        hmmsearch_paths = search.search(lock, pool)
        
    else:
        # go over the given hmm directory and obtain all subset alignment
        # get their retained columns with respect to the backbone alignment
        _dummy_search = SearchAlgorithm(None)
        subset_to_retained_columns, subset_to_nongaps_per_column \
                = _dummy_search.readHMMDirectory(lock, pool)

    # sanity check before moving on
    assert Configs.backbone_path != None \
            and os.path.exists(Configs.backbone_path), 'backbone alignment missing'
    assert Configs.backbone_tree_path != None \
            and os.path.exists(Configs.backbone_tree_path), 'backbone tree missing'
    assert Configs.query_path != None \
            and os.path.exists(Configs.query_path), 'query sequences missing'
    assert Configs.hmmdir != None \
            and os.path.isdir(Configs.hmmdir), 'eHMM directory missing'

    #print('--- Current memory usage: {} MB ---'.format(memoryUsage()))

    # 0) create a temporary (working) backbone alignment,
    # enforcing all letters to be upper-cases
    tmp_backbone_path, backbone_length = writeTempBackbone(
            Configs.outdir + '/tree_decomp/backbone', Configs.backbone_path)

    # 1) get all sub-queries, write to [outdir]/data
    num_seq, index_to_hmm, ranked_bitscores, sid_to_query_names, \
            sid_to_query_seqs, renamed_taxa = loadSubQueries(lock, pool)

    # 2) calculate weights, if needed
    if Configs.use_weight:
        print('\nCalculating weights...')
        taxon_to_weights = writeWeights(index_to_hmm, ranked_bitscores, pool)
    else:
        print('\nLoading bit-scores...')
        taxon_to_weights = writeBitscores(ranked_bitscores, pool)

    # 2a) if saving weights to local
    if Configs.save_weight:
        print('\n(user option) Writing weights to local...')
        writeWeightsToLocal(taxon_to_weights, Configs.outdir + '/weights.txt')

    # 3) solve each subset
    sub_alignment_paths = []
    if not os.path.isdir(Configs.outdir + '/temp'):
        os.makedirs(Configs.outdir + '/temp')

    ####### Addition - 6.13.2023 #######
    # Create lists that record retained column indexes and the
    # number of nongap characters per column
    # of each subset alignment for in-memory merging of HMM-query alignments,
    # instead of calling GCM (the WITCH-ng's way)
    # will be passed to the initiation of process pool as initial arguments

    # close pool and re-initiate pool to run the actual query alignment part
    Configs.warning('Closing ProcessPoolExecutor instance (for backbone)...')
    pool.shutdown()
    #pool = ProcessPoolExecutor(Configs.num_cpus,
    pool = WITCHProcessPoolExecutor(Configs.num_cpus,
            initializer=initiate_pool_query_alignment,
            initargs=(q, args, subset_to_retained_columns,
                subset_to_nongaps_per_column))
            #mp_context=mp.get_context('spawn'),
    Configs.warning('ProcessPoolExecutor instance re-opened (for query alignment).')

    # ProcessPoolExecutor version
    print('\nPerforming GCM alignments on query subsets...')
    index_list = [_i for _i in range(num_seq)]
    subset_query_names = [sid_to_query_names[_i] for _i in range(num_seq)]
    subset_query_seqs = [sid_to_query_seqs[_i] for _i in range(num_seq)]
    
    # deal with queries that has no weight
    subset_weights = []
    ignored_queries = []
    for _i in range(num_seq):
        _q = sid_to_query_names[_i]
        if _q in taxon_to_weights:
            subset_weights.append(taxon_to_weights[_q])
        else:
            Configs.warning('{} does not have any weight/bitscore '.format(
                _q) + 'and will be ignored the final alignment.')
            subset_weights.append(tuple())
            ignored_queries.append(_i)

    #func = partial(alignSubQueries, tmp_backbone_path, index_to_hmm, lock)
    func = partial(alignSubQueriesNew, tmp_backbone_path, backbone_length,
            index_to_hmm, lock)

    # try submit jobs one-by-one and collect results
    futures, success = [], []
    #while len(success) < (num_seq - len(ignored_queries)):
    for i in range(len(subset_query_names)):
        futures.append(pool.submit(func, subset_query_names[i],
            subset_query_seqs[i], subset_weights[i], index_list[i]))

    # iterate over jobs as they complete
    for future in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures), **tqdm_styles):
        _query, _index = future.result()
        # failed job, should be ignored in the output
        if len(_query) == 0:
            ignored_queries.append(subset_query_names[_index])
        else:
            success.append(_query)

    #results = list(pool.map(func, subset_query_names, subset_query_seqs,
    #    subset_weights, index_list, chunksize=Configs.chunksize))
    #retry_results, success, failure, ignored = [], [], [], []
    #while len(success) < (num_seq - len(ignored_queries)):
    #    success.extend([r for r, _i in results if r is not None])
    #    success.extend([r for r, _i in retry_results if r is not None])
    #    ignored.extend([':'.join(r.split(':')[1:]) 
    #                for r, _i in results if 'skipped:' in r])
    #    ignored.extend([':'.join(r.split(':')[1:])
    #                for r, _i in retry_results if 'skipped:' in r])

    #    failed_items = []; failed_item_queries = []; failed_item_weights = []
    #    while not q.empty():
    #        failed_items.append(q.get())
    #    if len(failed_items) > 0:
    #        Configs.log('Rerunning failed jobs: {}'.format(failed_items))
    #        failed_item_query_names = [sid_to_query_names[_i]
    #                                for _i in failed_items]
    #        failed_item_query_seqs = [sid_to_query_seqs[_i]
    #                                for _i in failed_items]
    #        failed_item_weights = [taxon_to_weights[_q]
    #                                for _q in failed_item_query_names]
    #        failure.append(failed_items)
    #        retry_results = list(pool.map(func, failed_item_query_names,
    #            failed_item_query_seqs, failed_item_weights, failed_items,
    #            chunksize=Configs.chunksize))
    queries = success

    # 4) merge all results 
    print('\nAll GCM subproblems finished! Doing merging with transitivity...')
    mergeAlignmentsCollapsed(tmp_backbone_path, queries, renamed_taxa, pool)

    # if there are any ignored queries, write them to local
    ignored_path = Configs.outdir + '/ignored_queries.fasta'
    Configs.log('Writing {} ignored sequences to local at {}'.format(
        len(ignored_queries), ignored_path)) 
    newname_to_oldname = {v: k for k, v in renamed_taxa.items()}
    if len(ignored_queries) > 0:
        with open(ignored_path, 'w') as f:
            for _i in ignored_queries:
                name, seq = sid_to_query_names[_i], sid_to_query_seqs[_i]
                if name in newname_to_oldname:
                    name = newname_to_oldname[name]
                f.write('>{}\n{}\n'.format(name, seq))

    Configs.warning('Closing ProcessPoolExecutor instance...')
    pool.shutdown()
    Configs.warning('ProcessPoolExecutor instance closed.')
    if not Configs.keeptemp:
        clearTempFiles()

    print('\nAll done!')
