import os, sys, time, itertools

from witch_msa.configs import *
from witch_msa.gcmm.aligner import alignSubQueries, alignSubQueriesNew
from witch_msa.gcmm.loader import readCheckpointAlignments
from witch_msa.gcmm.task import *
from witch_msa.gcmm.callback import callback_queryAlignment

from witch_msa.helpers.alignment_tools import ExtendedAlignment

import concurrent.futures
from functools import partial
from tqdm import tqdm

'''
Function to yield the list of arguments for a job/task to perform.
'''
#def getTasks(remaining_indexes,
#        subset_query_names, subset_query_seqs, subset_weights):
#    for i in range(len(remaining_indexes)):
#        _i = remaining_indexes[i]
#        yield subset_query_names[_i], subset_query_seqs[_i], \
#                subset_weights[_i], _i

'''
Function to handle a single future returned from the main query alignment job
    future.result() should have two fields:
        0 - query alignment in ExtendedAlignment
        1 - query index
    Return: runtime in handling the future object
'''
#def handleFuture(future, success, ignored, retry, checkpoint_path,
#        i_retry=0):
#    s1 = time.time()
#    _query, _index = future.result()
#    if not _query and i_retry > 0:
#        retry.append(_index)
#    else:
#        # failed job indicated in the <witch-ng> or <default> pipelines,
#        # should be ignored in the output
#        if (not _query) or len(_query) == 0:
#            ignored.append(subset_query_names[_index])
#        else:
#            # write to checkpoint file
#            writeOneCheckpointAlignment(checkpoint_path, _query)
#            success.append(_query)
#    return time.time() - s1

'''
Function to handle returned future objects from query-alignments
(alignSubQueries* functions).
    b_retry - whether to retry some failed indexes (indicated by their return types)
    i_retry - the remaining number of retries allowed

    return: (1) a list of indexes for retry
            (2) the remaining number of retries
'''
#def handleFutures(futures, subset_query_names, success, ignored, i_retry,
#        checkpoint_path):
#    runtime_handle_result = 0.
#    retry = []
#    for future in tqdm(
#            concurrent.futures.as_completed(futures),
#            total=len(futures), **tqdm_styles):
#        s1 = time.time()
#        _query, _index = future.result()
#
#        # In the case of the <default> pipeline, there is a chance
#        # that the subprocess failed (e.g., due to timeout).
#        # Need to re-queue these jobs, possibly with a longer timeout threshold
#        if not _query and i_retry > 0:
#            retry.append(_index)
#        else:
#            # failed job indicated in the <witch-ng> or <default> pipelines,
#            # should be ignored in the output
#            if (not _query) or len(_query) == 0:
#                ignored.append(subset_query_names[_index])
#            else:
#                # write to checkpoint file
#                writeOneCheckpointAlignment(checkpoint_path, _query)
#                success.append(_query)
#        runtime_handle_result += time.time() - s1
#
#    return retry, max(0, i_retry - 1), runtime_handle_result


'''
Function to handle submitting and collecting results for all query alignments.
Also handles writing to checkpoints (so their runtimes are recorded)
'''
def submitAndCollectFutures(pool, lock, sid_to_query_names, sid_to_query_seqs,
        num_seq, index_to_hmm, taxon_to_weights,
        tmp_backbone_path, backbone_length):
    handle_runtime = 0.

    # allow for 1 retry if a query alignment failed (most likely due to 
    # running old-witch). Re-run with WITCH-ng pipeline
    i_retry = 1

    # items to return 
    queries, ignored = [], []

    subset_query_names = [sid_to_query_names[_i] for _i in range(num_seq)]
    subset_query_seqs = [sid_to_query_seqs[_i] for _i in range(num_seq)]
    remaining_indexes = [_i for _i in range(num_seq)]

    # set up checkpoints for query-alignments, and load in any finished
    # query alignments from a previous run
    checkpoint_path = Configs.outdir + '/checkpoint_alignments.txt.gz'

    # load from existing checkpoint path, if any exists
    checkpoint_queries = {}
    if os.path.exists(checkpoint_path) and os.stat(checkpoint_path).st_size > 0:
        print('\nFound existing checkpoint query alignments: {}'.format(
            checkpoint_path))
        checkpoint_queries = readCheckpointAlignments(checkpoint_path, pool, lock)
        print('\tSuccessfully read in {} checkpoint query alignments...'.format(
            len(checkpoint_queries)))

        # append checkpoint queries to queries
        queries = list(checkpoint_queries.values())
    
    # deal with queries that needs to be skipped (has no weight) or have
    # completed (from checkpoint)
    subset_weights, completed_indexes = [], []
    for _i in range(num_seq):
        _q = sid_to_query_names[_i]
        if _q in checkpoint_queries:
            completed_indexes.append(_i)
            subset_weights.append(tuple())
            continue

        if _q in taxon_to_weights:
            subset_weights.append(taxon_to_weights[_q])
        else:
            Configs.warning('{} does not have any weight/bitscore '.format(
                _q) + 'and will be ignored the final alignment.')
            subset_weights.append(tuple())
            ignored.append(_i)
    remaining_indexes = list(set(remaining_indexes).difference(set(completed_indexes)))
    remaining_indexes = list(set(remaining_indexes).difference(set(ignored)))

    # early return if no remaining indexes
    if len(remaining_indexes) == 0:
        return queries, ignored

    # Set up either WITCH or WITCH-ng's way of aligning the query sequence
    func_map = {'old-witch': alignSubQueries, 'witch-ng': alignSubQueriesNew}
    if func_map[Configs.mode]:
        func = partial(func_map[Configs.mode], tmp_backbone_path,
                backbone_length, index_to_hmm, lock, Configs.timeout)
    else:
        raise NotImplementedError

    mytasks = getTasksWithIndexes(remaining_indexes, subset_query_names,
            subset_query_seqs, subset_weights)
    success, _ignored, retry, handle_runtime = runTasks(
            func, pool, mytasks, len(remaining_indexes),
            i_retry=i_retry,
            max_concurrent_jobs=Configs.max_concurrent_jobs,
            callback_func=callback_queryAlignment,
            callback_args=[checkpoint_path]
            )
    queries.extend(success); ignored.extend(_ignored)

    # queue up re-try if any
    retry_func = partial(func_map['witch-ng'], tmp_backbone_path,
            backbone_length, index_to_hmm, lock, Configs.timeout)
    while len(retry) > 0 and i_retry > 0:
        i_retry -= 1
        retry_tasks = getTasksWithIndexes(retry, subset_query_names,
                subset_query_seqs, subset_weights)
        success, _ignored, retry, handle_runtime = runTasks(
                retry_func, pool, retry_tasks, len(retry),
                i_retry=i_retry,
                max_concurrent_jobs=Configs.max_concurrent_jobs,
                callback_func=callback_queryAlignment,
                callback_args=[checkpoint_path]
                )
        queries.extend(success); ignored.extend(_ignored)

    
    # try submitting jobs one-by-one and collect results
    #futures, success, retry = [], [], []
    # UPDATE 8.4.2023 - try limiting the total number of existing jobs
    #                 - by pausing submission after queue is full
    #with tqdm(total=len(remaining_indexes), **tqdm_styles) as pbar:
    #    tasks_to_do = getTasks(remaining_indexes, subset_query_names,
    #            subset_query_seqs, subset_weights)
    #    futures = {
    #        pool.submit(func, *task): task[-1]
    #        for task in itertools.islice(tasks_to_do, Configs.max_concurrent_jobs)
    #    }
    #    while futures:
    #        # wait for the next future to complete
    #        done, _ = concurrent.futures.wait(
    #                futures, return_when=concurrent.futures.FIRST_COMPLETED)

    #        # process results
    #        for future in done:
    #            # allow re-adding some failed jobs back to queue
    #            handle_runtime += handleFuture(future, success, ignored,
    #                    retry, checkpoint_path, i_retry=1)
    #            original_task = futures.pop(future)
    #            #print(f"Output of {original_task} is {future.result()}")
    #        pbar.update(len(done))

    #        # schedule the next set of tasks, no more than the number
    #        # of done to only have Configs.max_concurrent_jobs in the pool
    #        for task in itertools.islice(tasks_to_do, len(done)):
    #            future = pool.submit(func, *task)
    #            futures[future] = task[-1]
    
    # run retries if any exist
    # retry will use WTICH-ng's way (which presumably should be faster in the
    # case the subprocess reaches <timeout> seconds
    #while len(retry) > 0:
    #    print('\tQueuing up {} failed alignments, {} retries remaining...'.format(
    #        len(retry), i_retry))
    #    retry_futures = []
    #    retry_func = partial(func_map['witch-ng'], tmp_backbone_path,
    #            backbone_length, index_to_hmm, lock, Configs.timeout)
    #    for i in retry:
    #        retry_futures.append(pool.submit(retry_func, subset_query_names[i],
    #            subset_query_seqs[i], subset_weights[i], i))

    #    retry, i_retry, _runtime = handleFutures(retry_futures,
    #            subset_query_names, success, ignored, i_retry,
    #            checkpoint_path)
    #    handle_runtime += _runtime
    #queries.extend(success)

    Configs.runtime(' '.join(['(submitAndCollectFutures) Time to collect',
                'futures and write checkpoints (s):',
                str(handle_runtime)]))
    return queries, ignored
