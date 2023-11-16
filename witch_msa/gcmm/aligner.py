import re, os, subprocess, psutil, shutil
import time
import math
from collections import defaultdict
from witch_msa.configs import Configs
from witch_msa.gcmm.weighting import readWeights, readBitscores
from witch_msa.helpers.alignment_tools import Alignment, ExtendedAlignment, compressInsertions
from multiprocessing import Queue#, Lock

## initialize the lock for asynchronous safe logging
## Lock is passed from the main process
#def init_lock(l):
#    global lock
#    lock = l

'''
helper function to obtain merged query alignment with insertions marked
'''
def getQueryAlignment(query_name, path, dtype='fasta'):
    query = ExtendedAlignment([])
    query_name, insertions, query_aligned_columns = \
            query.read_query_alignment(query_name, path, dtype)

    # only remove temp file if input type is a fasta file
    if dtype == 'fasta':
        os.system('rm {}'.format(path))
        query_aligned_columns = []
    return query, query_aligned_columns

'''
Helper function to generate backbone alignments for the constraint sets
'''
def getBackbones(index_to_hmm, taxon, taxon_ind, seq, query_path, sorted_weights, 
        workdir, backbone_dir, use_gcm=True):
    ret = dict()

    # for each taxon and its scores on HMMs, sort it in decreasing order
    #for taxon, sorted_weights in weights.items():
    #taxon = next(iter(unaligned)); seq = unaligned[taxon]
    weights_map = {ind: w for (ind, w) in sorted_weights}

    ret_str = ''
    if ret_str != '':
        ret_str += '\n'

    if len(sorted_weights) == 0:
        return 'N/A', None

    # during loading bit-scores/weights, we already trimmed weights/scores
    # so here we just use them
    top_k_hmms = sorted_weights
    if Configs.use_weight:
        pass
        #if Configs.weight_adjust == 'normalize':
        #    cur_total_w = sum([w[1] for w in top_k_hmms])
        #    top_k_hmms = [(w[0], w[1] * (1. / cur_total_w))
        #            for w in top_k_hmms]
        #elif Configs.weight_adjust == 'maxto1':
        #    max_w = top_k_hmms[0][1]
        #    top_k_hmms = [(w[0], w[1] / max_w) for w in top_k_hmms]
    else:
        top_k_hmms = [(w[0], 1) for w in top_k_hmms]
    ret_str += '{}\tpassed to main pipeline with top {} weights: {}'.format(
            taxon, len(top_k_hmms), top_k_hmms)
    for item in top_k_hmms:
        # denote that this HMM will be used
        ret[item[0]] = 1

    # now we can split the fragments into chunks, each chunk corresponds to
    # all fragments added to an HMM i
    # For each HMM, run HMMAlign on its fragment chunk to get a backbone
    # alignment
    weights_path = workdir + '/weights.txt'
    #weights_file = open(weights_path, 'w')
    to_write = ''

    #fn_to_weight = dict()
    
    ### 6.21.2023 - alternative getBackbone without writing to local (used in GCM)
    # return a list of aligned positions for each hmmalign extended alignment
    subset_to_aligned_columns = dict()
    
    # initialize the backbone index at backbone_start_index
    bb_index = 0
    for i in ret.keys():
        if not os.path.isdir(workdir):
            os.makedirs(workdir)

        # [NEW] save each single sequence to a fasta
        this_hmm = index_to_hmm[i].hmm_model_path
    
        # hmmalign
        hmmalign_result_path = '{}/hmmalign.results.{}.{}.out'.format(
                workdir, taxon_ind, i)
        cmd = '{} -o {} {} {}'.format(Configs.hmmalignpath,
                hmmalign_result_path, this_hmm, query_path)
        os.system(cmd)
        #subprocess.run(cmd.split(' '))
    
        # Extended alignment
        ap_aln = ExtendedAlignment([taxon])

        # if subsequently using GCM/MAGUS to merge query-HMM alignments
        # need to save the files to local with weights.txt
        if use_gcm:
            ap_aln.build_extended_alignment(index_to_hmm[i].alignment_path,
                    [hmmalign_result_path], True)
            for key in ap_aln.keys():
                ap_aln[key] = ap_aln[key].upper()
            # save extended alignment to backbone directory
            this_bb_path = backbone_dir + '/{}_{}.fasta'.format(taxon_ind, i)
            ap_aln.write(this_bb_path, 'FASTA')

            # add the weights of the bb to weights.txt
            if Configs.use_weight:
                real_this_bb_path = os.path.realpath(this_bb_path)
                to_write += '{},{}\n'.format(real_this_bb_path,
                    weights_map[i])
                #weights_file.write('{},{}\n'.format(
                #        real_this_bb_path, weights_map[i]))
                #fn_to_weight[real_this_bb_path] = weights_map[i]
        else:
            ap_aln.read_extended_alignment(hmmalign_result_path)
            query_aligned_columns = []
            regular = 0
            for _i in range(len(ap_aln[taxon])):
                char = ap_aln[taxon][_i]
                if char == '-':
                    regular += 1
                    continue
                elif char.islower():
                    query_aligned_columns.append(-1)
                else:
                    query_aligned_columns.append(regular)
                    regular += 1

            #query_only_aln, _query_aligned_columns = \
            #        getQueryAlignment(taxon, ap_aln, dtype='mutablealignment')
            subset_to_aligned_columns[i] = query_aligned_columns

    if use_gcm:
        with open(weights_path, 'w') as weights_file:
            weights_file.write(to_write)
    
    return ret_str, weights_map, subset_to_aligned_columns #fn_to_weight

'''
The core function of GCM+eHMMs. Merge the HMMAlign results of the set of
queries to the input alignment.

The main differences between GCM+eHMMs and UPP are: 1) GCM+eHMMs uses adjusted
bitscore, which takes the number of queries in an HMM into consideration. 2)
GCM+eHMMs can utilize more than one HMM (while UPP uses the best HMM based on
bitscore) to align the queries; hence, more information is used.
'''
def alignSubQueries(backbone_path, backbone_length, index_to_hmm, lock, timeout,
        taxon, seq, query_weights, index):
    # index maps to query in unaligned 
    # query_weights in the form ((hmmX, scoreX), (hmmY, scoreY), ...)

    s11 = time.time()
    # add constraints
    constraints_dir = Configs.outdir + '/constraints/{}'.format(index)
    if not os.path.isdir(constraints_dir):
        os.makedirs(constraints_dir)
    # 0-th constraint set comes from the input alignment
    shutil.copyfile(backbone_path, '{}/c0.fasta'.format(
        constraints_dir))

    # only one query sequence in unaligned -> our second constraint
    unaligned = Alignment()
    unaligned[taxon] = seq
    query_path = constraints_dir + '/c1.fasta'
    #taxon = next(iter(unaligned))
    #unaligned[taxon] = unaligned[taxon].upper()
    unaligned.write(query_path, 'FASTA')

    #c_index = 1
    #for name in unaligned.keys():
    #    single_frag = unaligned.sub_alignment([name])
    #    # enforce upper case letters for all bases
    #    single_frag[name] = single_frag[name].upper()
    #    single_frag.write('{}/c{}.fasta'.format(constraints_dir, c_index),
    #            'FASTA')
    #    c_index += 1
    time_obtain_constraints = time.time() - s11

    s12 = time.time()
    # backbone alignments for GCM
    bb_dir = Configs.outdir + '/backbone_alignments/{}'.format(index)
    if not os.path.isdir(bb_dir):
        os.makedirs(bb_dir)
    # hmmsearch directory
    search_dir = Configs.outdir + '/search_results/{}'.format(index)
    if not os.path.isdir(search_dir):
        os.makedirs(search_dir)

    # get backbones with the information we have
    weights_str, weights_map, _ = getBackbones(index_to_hmm,
            taxon, index, seq, query_path,
            query_weights, search_dir, bb_dir)
    time_obtain_backbones = time.time() - s12

    #if weights_map:
    #    list_of_bb = list(weights_map.keys())
    #    list_of_weights = list(weights_map.values())
    #else:
    #    list_of_bb, list_of_weights = [], []

    s13 = time.time()

    # time out the process if exceeding [timeout] seconds
    task_timed_out = False
    gcm_outdir = Configs.outdir + '/magus_outputs/{}'.format(index)
    gcmpath = Configs.gcm_path

    est_path = None
    if weights_str != 'N/A':
        # run GCM (modified MAGUS which takes in weights) on the subset
        est_path = Configs.outdir + '/temp/magus_result_{}.txt'.format(index)
        cmd = [Configs.magus_path, '-np', '1',
                '-d', gcm_outdir, '-s', constraints_dir, '-b', bb_dir,
                '-o', est_path, '-f', str(Configs.inflation_factor),
                '--graphclustermethod', Configs.graphclustermethod,
                '--graphtracemethod', Configs.graphtracemethod,
                '--graphtraceoptimize', Configs.graphtraceoptimize]
        #import glob
        #files = glob.glob(f'{constraints_dir}/*')
        #cmd = [gcmpath, 'merge',
        #        '-i', *files, '-g', *list_of_bb,
        #        '-o', est_path, '-w', *[str(w) for w in list_of_weights]]
        #print(' '.join(cmd))
        # use macOS version mcl (version 21.257) if system is macOS
        if Configs.mclpath:
            cmd += ['--mclpath', Configs.mclpath]
        if Configs.use_weight:
            weights_path = search_dir + '/weights.txt'
            cmd += ['-w', weights_path] 
        #os.system(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        try:
            p.wait(timeout)
        except subprocess.TimeoutExpired:
            task_timed_out = True
            parent = psutil.Process(p.pid)
            for c in parent.children(recursive=True):
                c.kill()
            parent.kill()
        # additional check for called process error
        except subprocess.CalledProcessError:
            ret_code = p.returncode
            res = p.communicate()
            stderr_out = res[1]

            err_to_write = 'Encountered subprocess.CalledProcessError!' \
                + '\nerror code: {}'.format(ret_code) \
                + '\ndetailed call process output:'
            print(err_to_write)
            Configs.error(err_to_write)
            for line in res[0].decode(encoding='utf-8').split('\n'):
                print(line)
                Configs.error(line)
            exit(1)
    else:
        # no matching HMMs to align the given query sequence, directly
        # return
        Configs.warning('{} in task #{}'.format(taxon, index) \
                + ' does not have any matching HMMs, skipping...')
        return ExtendedAlignment([]), index, taxon

    # remove temp folders
    if not Configs.keeptemp:
        dirs_to_delete = [search_dir, bb_dir, constraints_dir]
        for _dir in dirs_to_delete:
            if os.path.isdir(_dir):
                os.system('rm -rf {}'.format(_dir))
    if not Configs.keepgcmtemp and os.path.isdir(gcm_outdir):
        os.system('rm -rf {}'.format(gcm_outdir))
    time_gcm = time.time() - s13
    
    if not task_timed_out:
        # read in the extended alignment of the query taxon
        s14 = time.time()
        if os.path.exists(est_path):
            query, _ = getQueryAlignment(taxon, est_path)
        else:
            # put the task back to queue
            lock.acquire()
            try:
                Configs.warning(
                        'Task #{}->{} finished within time but ' \
                        'has no output, queuing for retry...'.format(
                            index, taxon))
                alignSubQueries.q.put(index)
            finally:
                lock.release()
            return None, index, taxon

        time_merged_query = time.time() - s14

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
            Configs.runtime(' '.join(['(alignSubQueries,',
                    'i={}) Time to obtain the merged alignment'.format(index),
                    ' (s):', str(time_merged_query)]))
            Configs.debug("[MAGUS] Command used: {}".format(' '.join(cmd)))
            Configs.log(weights_str)
            #Configs.log('{} passed to main pipeline...'.format(taxon))
        finally:
            lock.release()
        return query, index, taxon
    else:
        lock.acquire()
        try:
            Configs.debug("[MAGUS] Command used: {}".format(' '.join(cmd)))
            Configs.warning('Task #{}->{} failed, queuing for retry...'.format(
                index, taxon))
            alignSubQueries.q.put(index)
        finally:
            lock.release()
        return None, index, taxon

'''
New query-HMM alignment merging function that does not rely on GCM/MAGUS.
This one will use the knowledge of the index mapping between each sub-alignment
to the backbone alignment to create an alignment graph, and add weights to edges
between nodes (of backbone and a query).

Then, the merged alignment for a query is obtained by running Smith-Waterman
on the alignment graph (|C|x|q|, where |C| is the # sites of the backbone and 
|q| the length of the query sequence.

All of the above can be achieved without calls to subprocesses such as MAGUS
or the experimental GCM code (by Baqiao), and most operations can be done
in-memory (for speed).
'''
def alignSubQueriesNew(backbone_path, backbone_length, index_to_hmm, lock, timeout,
        taxon, seq, query_weights, index):
    ############## STEP 1 ###############
    # obtain backbones (edges without weights yet)
    s11 = time.time()
    # query constraint to local (for hmmsearch)
    constraints_dir = Configs.outdir + '/constraints/{}'.format(index)
    if not os.path.isdir(constraints_dir):
        os.makedirs(constraints_dir)

    unaligned = Alignment(); unaligned[taxon] = seq
    query_path = constraints_dir + '/c1.fasta'
    unaligned.write(query_path, 'FASTA')

    # hmmsearch directory
    search_dir = Configs.outdir + '/search_results/{}'.format(index)
    if not os.path.isdir(search_dir):
        os.makedirs(search_dir)
    weights_str, subset_to_weight, subset_to_aligned_columns = \
            getBackbones(index_to_hmm, taxon, index, seq, query_path,
                    query_weights, search_dir, None, use_gcm=False)
    #print(taxon, subset_to_weight[next(iter(subset_to_weight))],
    #        subset_to_aligned_columns[next(iter(subset_to_aligned_columns))])
    time_obtain_backbones = time.time() - s11

    s12 = time.time()
    query = ExtendedAlignment([])
    if weights_str != 'N/A':
        ############## STEP 2 ###############
        # added weighted edges to the alignment graph from gluing alignment
        # subset_to_retained_columns -> globally shared/used
        # subset_to_nongaps_per_column -> globally shared/used
        # subset_to_weight -> k=10 then 10 subsets and their weights
        # subset_to_aligned_columns -> aligned columns in subsets for this query

        # technically we are just creating a matrix for DP later, and adding in
        # corresponding weighted edges as they are read in
        alignment_graph = [[0. for _j in range(backbone_length + 1)]
                              for _i in range(len(seq) + 1)]
        # key -> (index of query char, index of backbone col)
        # val -> accumulated weights of the key
        combined_weights = defaultdict(float)

        subset_to_retained_columns = alignSubQueriesNew.subset_to_retained_columns
        subset_to_nongaps_per_column = alignSubQueriesNew.subset_to_nongaps_per_column
        # also record the minimum/maximum of indexes of backbone used
        # this information can be used to limit the DP search space (for speed-up)
        min_col_ind, max_col_ind = backbone_length + 1, -1

        for subset, aligned_columns in subset_to_aligned_columns.items():
            # subset -> a bit complicated to explain here how it is used 
            # indexes of aligned_columns -> row index
            for i in range(len(aligned_columns)):
                subset_col = aligned_columns[i]
                if subset_col == -1:    # skip insertion site of query
                    continue
                
                # compute the column index to add edge weight
                #print(subset, subset_col, len(subset_to_retained_columns[subset]))
                j = subset_to_retained_columns[subset][subset_col]
                
                # compute the weight to accumulated to combined_weights
                weight_to_add = subset_to_nongaps_per_column[subset][subset_col] * \
                        subset_to_weight[subset]
                combined_weights[(i, j)] += weight_to_add

                # update min/max index used for the column (backbone)
                min_col_ind = min(min_col_ind, j)
                max_col_ind = max(max_col_ind, j)
        #print(taxon, seq, combined_weights)
        #exit(0)

        ############## STEP 2.5 ###############
        # run Smith-Waterman on the alignment graph to obtain the highest weighted
        # trace (thus, alignment of the query to the backbone)
        # use the min/max col indexes to bound backtrace region
        backtrace = [[0 for _j in range(backbone_length + 1)]
                              for _i in range(len(seq) + 1)]

        for i in range(0, len(seq) + 1):
            for j in range(min_col_ind, max_col_ind + 2):
                if i == 0 or j == min_col_ind:
                    alignment_graph[i][j] = 0.
                    continue

                cur_max = 0.; cur_bt = 0
                cw = combined_weights.get((i-1, j-1), 0.)
                values = [alignment_graph[i-1][j-1] + cw,
                          alignment_graph[i-1][j],
                          alignment_graph[i][j-1]]
                for _ind, val in enumerate(values):
                    if _ind == 0 and cw <= 0:
                        cur_bt = 1
                        continue
                    if val > cur_max:
                        cur_max = val
                        cur_bt = _ind
                alignment_graph[i][j] = cur_max
                backtrace[i][j] = cur_bt

        # obtain alignment from backtrace
        # the alignment is represented by the set of characters of the query
        result = []; i, j = len(seq), max_col_ind + 1
        while i > 0 and j > min_col_ind:
            # diagonal - match
            bt = backtrace[i][j]
            if bt == 0:    # a match
                result.append(seq[i-1])
                i -= 1; j -= 1
            elif bt == 1:  # a query "insertion", use lowercase
                result.append(seq[i-1].lower())
                i -= 1
            elif bt == 2:  # a query "deletion", use a gap '-'
                result.append('-')
                j -= 1
            else:
                raise ValueError

        # fill up remaining positions (i.e., i and j both to 0)
        # only one "while" will happen
        while i > 0:
            result.append(seq[i-1].lower()); i -= 1
        while j > min_col_ind:
            result.append('-'); j -= 1

        # reverse the result to get correct order
        result = result[::-1]
        #print(taxon, result)

        # append '-' to start, end to fill up the alignment (since we use min/max
        # of backbone columns to reduce search space)
        result = ['-'] * min_col_ind + result + \
                ['-'] * (backbone_length - max_col_ind - 1)

        # compress insertions at front/end of the query alignment
        combined = compressInsertions(''.join(result))

        # return an extended alignment object with just the query
        # sequence and its updated indexes
        query[taxon] = combined; query._reset_col_names()
        insertion = -1; regular = 0
        for i in range(len(combined)):
            if combined[i].islower():
                query._col_labels[i] = insertion; insertion -= 1
            else:
                query._col_labels[i] = regular; regular += 1

    # remove temp folders
    if not Configs.keeptemp:
        dirs_to_delete = [search_dir, constraints_dir]
        for _dir in dirs_to_delete:
            if os.path.isdir(_dir):
                os.system('rm -rf {}'.format(_dir))
    time_alignment_trace = time.time() - s12

    # if query is empty, meaning there is no weight and the whole algorithm
    # did not run. Simply return the empty alignment obj
    if len(query) == 0:
        lock.acquire()
        try:
            Configs.warning('{} in task #{}'.format(taxon, index) \
                    + ' does not have any matching HMMs, ignored in final output...')
        finally:
            lock.release()
        return query, index, taxon
    
    # sanity check for finishing the query-HMM alignment merge
    if query.get_length() >= backbone_length:
        lock.acquire()
        try:
            Configs.runtime(' '.join(['(alignSubQueriesNew,',
                    'i={}) Time to obtain constraint/backbones'.format(index),
                    ' (s):', str(time_obtain_backbones)]))
            Configs.runtime(' '.join(['(alignSubQueriesNew,',
                    'i={}) Time to run alignment and clean-up'.format(index),
                    ' (s):', str(time_alignment_trace)]))
            Configs.log(weights_str)
        finally:
            lock.release()
        return query, index, taxon
    else:
        lock.acquire()
        try:
            Configs.warning('Task #{}->{} failed, ignored in final output...'.format(
                index, taxon))
        finally:
            lock.release()
        # return an empty alignment as indication of failure
        return ExtendedAlignment([]), index, taxon
