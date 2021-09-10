import os, math, psutil, shutil 
from configs import Configs
from gcmm.loader import loadSubQueries 
from gcmm.weighting import writeWeights, writeBitscores
from gcmm.aligner import alignSubQueries
from gcmm.merger import mergeAlignments
from helpers.alignment_tools import Alignment

from multiprocessing import Pool, Lock, Queue, Manager
from concurrent.futures.process import ProcessPoolExecutor
from functools import partial

'''
Delete all unnecessary intermediate files
'''
def clearTempFiles():
    if not Configs.keepsubalignment:
        os.system('rm -r {}/magus_result_*'.format(Configs.outdir))
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
Main process for GCM+eHMMs
'''
def mainAlignmentProcess():
    m = Manager()
    lock = m.Lock()
    #l = Lock()
    q = Queue()

    # initialize the main pool at the start so that it's not memory
    # intensive
    Configs.warning('Initializing ProcessorPoolExecutor instance...')
    pool = ProcessPoolExecutor(Configs.num_cpus, initializer=init_queue,
            initargs=(q,))
    _ = pool.submit(dummy)

    # 1) get all sub-queries, write to [outdir]/data
    num_subset, index_to_hmm, ranked_bitscores = loadSubQueries()

    # 2) calculate weights, if needed 
    if Configs.use_weight:
        writeWeights(index_to_hmm, ranked_bitscores, pool)
    else:
        writeBitscores(ranked_bitscores, pool)

    # 3) solve each subset
    sub_alignment_paths = []

    ############ multiprocessing with Pool ##########
    # manager version
    #pool = Pool(Configs.num_cpus, initializer=init_queue, initargs=(q,))
    #index_list = [i for i in range(num_subset)]
    #func = partial(alignSubQueries, index_to_hmm, lock)
    #sub_alignment_paths = pool.map(func, index_list)
    #pool.close()
    #pool.join()

    # ProcessPoolExecutor version
    index_list = [i for i in range(num_subset)]
    func = partial(alignSubQueries, index_to_hmm, lock)
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
    print("\nAll GCM subproblems finished! Doing merging with transitivity...")
    mergeAlignments(sub_alignment_paths, pool)

    Configs.warning('Closing ProcessPoolExecutor instance...')
    pool.shutdown()
    Configs.warning('ProcessPoolExecutor instance closed.')
    if not Configs.keeptemp:
        clearTempFiles()
