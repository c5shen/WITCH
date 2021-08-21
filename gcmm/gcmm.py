import os, math 
from configs import Configs
from gcmm.loader import loadSubQueries 
from gcmm.weighting import loadWeights, Weights
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
        os.system('rm {}/magus_result_*'.format(Configs.outdir))
    if os.path.isdir('{}/backbone_alignments'.format(Configs.outdir)):
        os.system('rm -r {}/backbone_alignments'.format(Configs.outdir))
    if os.path.isdir('{}/constraints'.format(Configs.outdir)):
        os.system('rm -r {}/constraints'.format(Configs.outdir))
    if os.path.isdir('{}/search_results'.format(Configs.outdir)):
        os.system('rm -r {}/search_results'.format(Configs.outdir))
    if os.path.isdir('{}/data'.format(Configs.outdir)):
        os.system('rm -r {}/data'.format(Configs.outdir))
    if not Configs.keepgcmtemp \
            and os.path.isdir('{}/magus_outputs'.format(Configs.outdir)):
        os.system('rm -r {}/magus_outputs'.format(Configs.outdir))

'''
Init function for a queue
'''
def init_queue(q):
    alignSubQueries.q = q

'''
Main process for GCM+eHMMs
'''
def mainAlignmentProcess():
    m = Manager()
    lock = m.Lock()
    #l = Lock()

    # 1) get all sub-queries, write to [outdir]/data
    unaligned = Alignment(); unaligned.read_file_object(Configs.query_path)
    num_subset = max(1, math.ceil(len(unaligned) / Configs.subset_size))
    index_to_hmm, ranked_bitscores = loadSubQueries(unaligned)

    # 2) calculate weights, if needed 
    Weights.ranked_bitscores = ranked_bitscores
    if Configs.use_weight:
        loadWeights(index_to_hmm, ranked_bitscores)

    # 3) solve each subset
    sub_alignment_paths = []
    q = Queue()

    ############ multiprocessing with Pool ##########
    # manager version
    #pool = Pool(Configs.num_cpus, initializer=init_queue, initargs=(q,))
    #index_list = [i for i in range(num_subset)]
    #func = partial(alignSubQueries, index_to_hmm, lock)
    #sub_alignment_paths = pool.map(func, index_list)
    #pool.close()
    #pool.join()

    # ProcessPoolExecutor version
    pool = ProcessPoolExecutor(Configs.num_cpus, initializer=init_queue,
            initargs=(q,))
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
    Configs.warning('Closing ProcessPoolExecutor instance...')
    pool.shutdown()
    Configs.warning('ProcessPoolExecutor instance closed.')

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
    mergeAlignments(sub_alignment_paths)

    if not Configs.keeptemp:
        clearTempFiles()
