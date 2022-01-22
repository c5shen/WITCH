'''
Created on 1.19.2022 by Chengze Shen

Main pipeline of WITCH
'''

import os, math, psutil, shutil 
from configs import Configs
from gcmm.algorithm import DecompositionAlgorithm, SearchAlgorithm
from gcmm.loader import loadSubQueries 
from gcmm.weighting import writeWeights, writeBitscores
from gcmm.aligner import alignSubQueries
from gcmm.backbone import BackboneJob
from gcmm.merger import mergeAlignments

from helpers.alignment_tools import Alignment

import multiprocessing as mp
from multiprocessing import Lock, Queue, Manager#, Pool
from concurrent.futures.process import ProcessPoolExecutor
from functools import partial

'''
Delete all unnecessary intermediate files
'''
def clearTempFiles():
    if not Configs.keepsubalignment:
        os.system('find {}/temp/ -type f -delete'.format(Configs.outdir))
        shutil.rmtree('{}/temp'.format(Configs.outdir))
    if os.path.isdir('{}/backbone_alignments'.format(Configs.outdir)):
        shutil.rmtree('{}/backbone_alignments'.format(Configs.outdir))
    if os.path.isdir('{}/constraints'.format(Configs.outdir)):
        shutil.rmtree('{}/constraints'.format(Configs.outdir))
    if os.path.isdir('{}/search_results'.format(Configs.outdir)):
        shutil.rmtree('{}/search_results'.format(Configs.outdir))
    if os.path.isdir('{}/data'.format(Configs.outdir)):
        shutil.rmtree('{}/data'.format(Configs.outdir))
    if os.path.isdir('{}/weights'.format(Configs.outdir)):
        shutil.rmtree('{}/weights'.format(Configs.outdir))
    if os.path.isdir('{}/bitscores'.format(Configs.outdir)):
        shutil.rmtree('{}/bitscores'.format(Configs.outdir))
    if not Configs.keepgcmtemp \
            and os.path.isdir('{}/magus_outputs'.format(Configs.outdir)):
        shutil.rmtree('{}/magus_outputs'.format(Configs.outdir))

'''
Init function for a queue
'''
def init_queue(q):
    alignSubQueries.q = q

'''
Dummy function
'''
def dummy():
    pass

'''
Main process for WITCH 
'''
def mainAlignmentProcess():
    m = Manager()
    lock = m.Lock()
    #l = Lock()
    q = Queue()

    # initialize the main pool at the start so that it's not memory
    # intensive
    Configs.warning('Initializing ProcessorPoolExecutor instance...')
    pool = ProcessPoolExecutor(Configs.num_cpus,
            mp_context=mp.get_context('fork'),
            initializer=init_queue, initargs=(q,))
    _ = pool.submit(dummy)

    # 0) obtain the backbone alignment/tree and eHMMs
    # If no UPP eHMM directory provided, decompose from the backbone
    if not Configs.hmmdir:
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
        decomp = DecompositionAlgorithm(Configs.backbone_path,
                Configs.backbone_tree_path)
        hmmbuild_paths = decomp.decomposition(lock, pool)
        search = SearchAlgorithm(hmmbuild_paths)
        hmmsearch_paths = search.search(lock, pool)
        
        # default to <outdir>/tree_comp/root
        Configs.hmmdir = Configs.outdir + '/tree_decomp/root'

    # sanity check before moving on
    assert Configs.backbone_path != None \
            and os.path.exists(Configs.backbone_path), 'backbone alignment missing'
    assert Configs.backbone_tree_path != None \
            and os.path.exists(Configs.backbone_tree_path), 'backbone tree missing'
    assert Configs.query_path != None \
            and os.path.exists(Configs.query_path), 'query sequences missing'
    assert Configs.hmmdir != None \
            and os.path.isdir(Configs.hmmdir), 'eHMM directory missing'


    # 1) get all sub-queries, write to [outdir]/data
    num_subset, index_to_hmm, ranked_bitscores = loadSubQueries(lock, pool)

    # 2) calculate weights, if needed
    if Configs.use_weight:
        print('\nCalculating weights...')
        writeWeights(index_to_hmm, ranked_bitscores, pool)
    else:
        writeBitscores(ranked_bitscores, pool)

    # 3) solve each subset
    sub_alignment_paths = []
    if not os.path.isdir(Configs.outdir + '/temp'):
        os.makedirs(Configs.outdir + '/temp')

    ############ multiprocessing with Pool ##########
    # manager version
    #pool = Pool(Configs.num_cpus, initializer=init_queue, initargs=(q,))
    #index_list = [i for i in range(num_subset)]
    #func = partial(alignSubQueries, index_to_hmm, lock)
    #sub_alignment_paths = pool.map(func, index_list)
    #pool.close()
    #pool.join()

    # ProcessPoolExecutor version
    print('\nPerforming GCM alignments on query subsets...')
    index_list = [i for i in range(num_subset)]
    func = partial(alignSubQueries, Configs.backbone_path,
            index_to_hmm, lock)
    results = list(pool.map(func, index_list))
    retry_results, success, failure = [], [], []
    while len(success) < num_subset:
        success.extend([r for r in results if not r is None])
        success.extend([r for r in retry_results if not r is None])
        failed_items = []
        while not q.empty():
            failed_items.append(q.get())
        if failed_items:
            Configs.log('Rerunning failed jobs: {}'.format(failed_items))
            failure.append(failed_items)
            retry_results = list(pool.map(func, failed_items))
    sub_alignment_paths = success

    # global lock version
    #pool = Pool(Configs.num_cpus, initializer=init_lock, initargs=(l))
    #index_list = [i for i in range(num_subset)]
    #func = partial(alignSubQueries, index_to_hmm)
    #sub_alignment_paths = pool.map(func, index_list)
    #pool.close()
    #pool.join()

    ############ sequential version #################
    #for i in range(num_subset):
    #    output_path = alignSubQueries(i, index_to_hmm)
    #    sub_alignment_paths.append(os.path.abspath(output_path))

    # 4) merge all results 
    print('\nAll GCM subproblems finished! Doing merging with transitivity...')
    mergeAlignments(sub_alignment_paths, pool)

    Configs.warning('Closing ProcessPoolExecutor instance...')
    pool.shutdown()
    Configs.warning('ProcessPoolExecutor instance closed.')
    if not Configs.keeptemp:
        clearTempFiles()

    print('\nAll done!')
