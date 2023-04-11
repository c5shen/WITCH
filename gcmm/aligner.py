import os, subprocess, psutil, shutil
import time
from collections import defaultdict
from configs import Configs
from gcmm.weighting import readWeights, readBitscores
from helpers.alignment_tools import Alignment, ExtendedAlignment
from multiprocessing import Queue#, Lock

## initialize the lock for asynchronous safe logging
## Lock is passed from the main process
#def init_lock(l):
#    global lock
#    lock = l

'''
Helper function to generate backbone alignments for the constraint sets
'''
def getBackbones(index_to_hmm, taxon, seq, query_path, sorted_weights, 
        workdir, backbone_dir):
    ret = {} 

    # for each taxon and its scores on HMMs, sort it in decreasing order
    #for taxon, sorted_weights in weights.items():
    #taxon = next(iter(unaligned)); seq = unaligned[taxon]
    weights_map = {ind: w for (ind, w) in sorted_weights}

    ret_str = ''
    if ret_str != '':
        ret_str += '\n'

    #weights_map[taxon] = sorted_weights
    ## load weights from local
    #if Configs.use_weight:
    #    sorted_weights, this_weights_map = readWeights(taxon)
    #    weights_map[taxon] = this_weights_map
    #else:
    #    sorted_weights = readBitscores(taxon)

    if len(sorted_weights) == 0:
        return 'N/A', None

    # during loading bit-scores/weights, we already trimmed weights/scores
    # so here we just use them
    top_k_hmms = sorted_weights
    if Configs.use_weight:
        if Configs.weight_adjust == 'normalize':
            cur_total_w = sum([w[1] for w in top_k_hmms])
            top_k_hmms = [(w[0], w[1] * (1. / cur_total_w))
                    for w in top_k_hmms]
        elif Configs.weight_adjust == 'maxto1':
            max_w = top_k_hmms[0][1]
            top_k_hmms = [(w[0], w[1] / max_w) for w in top_k_hmms]
    else:
        top_k_hmms = [(w[0], 1) for w in top_k_hmms]
    ret_str += 'weights for {}: {}'.format(taxon, top_k_hmms)
    for item in top_k_hmms:
        # denote that this HMM will be used
        ret[item[0]] = 1

    # now we can split the fragments into chunks, each chunk corresponds to
    # all fragments added to an HMM i
    # For each HMM, run HMMAlign on its fragment chunk to get a backbone
    # alignment
    weights_path = workdir + '/weights.txt'
    weights_file = open(weights_path, 'w')

    fn_to_weight = {}
    
    # initialize the backbone index at backbone_start_index
    bb_index = 0
    for i in ret.keys():
        if not os.path.isdir(workdir):
            os.makedirs(workdir)
        #if not os.path.isdir(workdir + '/fragments'):
        #    os.makedirs('{}/fragments'.format(workdir))

        # [NEW] save each single sequence to a fasta
        this_hmm = index_to_hmm[i].hmm_model_path

        #query_path = '{}/{}.fasta'.format(workdir, taxon)
        #unaligned.write(query_path, 'FASTA')
    
        # hmmalign
        hmmalign_result_path = '{}/hmmalign.results.{}.{}.out'.format(
                workdir, taxon, i)
        cmd = '{} -o {} {} {}'.format(Configs.hmmalignpath,
                hmmalign_result_path, this_hmm, query_path)
        #if not (os.path.exists(hmmalign_result_path) and os.path.getsize(
        #    hmmalign_result_path) > 0):
        #    os.system(cmd)
        subprocess.run(cmd.split(' '))
    
        # Extended alignment
        ap_aln = ExtendedAlignment([taxon])
        ap_aln.build_extended_alignment(index_to_hmm[i].alignment_path,
                [hmmalign_result_path], True)
        for key in ap_aln.keys():
            ap_aln[key] = ap_aln[key].upper()
        # save extended alignment to backbone directory
        this_bb_path = backbone_dir + '/{}_{}.fasta'.format(taxon, i)
        ap_aln.write(this_bb_path, 'FASTA')

        # add the weights of the bb to weights.txt
        if Configs.use_weight:
            #real_this_bb_path = os.popen('realpath -s {}'.format(
            #    this_bb_path)).read().split('\n')[0]
            real_this_bb_path = os.path.realpath(this_bb_path)
            weights_file.write('{},{}\n'.format(
                    real_this_bb_path, weights_map[i]))
            fn_to_weight[real_this_bb_path] = weights_map[i]
    weights_file.close()
    
    return ret_str, fn_to_weight

'''
helper function to obtain merged query alignment with insertions marked
'''
def getQueryAlignment(query_name, path):
    query = ExtendedAlignment([])
    query_name, insertions = query.read_query_alignment(query_name, path)

    only_query = ExtendedAlignment([])
    only_query._col_labels = query._get_col_labels()
    only_query[query_name] = query[query_name]
    del query
    os.system('rm {}'.format(path))
    return only_query

'''
The core function of GCM+eHMMs. Merge the HMMAlign results of the set of
queries to the input alignment.

The main differences between GCM+eHMMs and UPP are: 1) GCM+eHMMs uses adjusted
bitscore, which takes the number of queries in an HMM into consideration. 2)
GCM+eHMMs can utilize more than one HMM (while UPP uses the best HMM based on
bitscore) to align the queries; hence, more information is used.
'''
def alignSubQueries(backbone_path, index_to_hmm, lock,
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
    weights_str, weights_map = getBackbones(index_to_hmm, taxon, seq, query_path,
            query_weights, search_dir, bb_dir)
    time_obtain_backbones = time.time() - s12

    if weights_map:
        list_of_bb = list(weights_map.keys())
        list_of_weights = list(weights_map.values())
    else:
        list_of_bb, list_of_weights = [], []

    s13 = time.time()

    # time out the process if exceeding [timeout] seconds
    task_timed_out = False
    gcm_outdir = Configs.outdir + '/magus_outputs/{}'.format(index)
    gcmpath = Configs.gcm_path

    est_path = None
    if weights_str != 'N/A':
        # run GCM (modified MAGUS which takes in weights) on the subset
        est_path = Configs.outdir + '/temp/magus_result_{}.txt'.format(index)
        cmd = ['python3', Configs.magus_path, '-np', '1',
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
        if Configs.mclpath is not None:
            cmd += ['--mclpath', Configs.mclpath]
        if Configs.use_weight:
            weights_path = search_dir + '/weights.txt'
            cmd += ['-w', weights_path] 
        #os.system(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        try:
            p.wait(Configs.timeout)
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

            print('Encountered subprocess.CalledProcessError!' +
                '\nerror code: {}'.format(ret_code) +
                '\ndetailed call process output:')
            for line in res[0].decode(encoding='utf-8').split('\n'):
                print(line)
            exit(1)
    else:
        # no matching HMMs to align the given query sequence, directly
        # return
        Configs.warning('{} in task #{}'.format(taxon, index) \
                + ' does not have any matching HMMs, skipping...')
        return 'skipped'

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
            query = getQueryAlignment(taxon, est_path)
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
            return None

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
            Configs.log('{} passed to main pipeline...'.format(taxon))
        finally:
            lock.release()
        return query
    else:
        lock.acquire()
        try:
            Configs.warning('Task #{}->{} failed, queuing for retry...'.format(
                index, taxon))
            alignSubQueries.q.put(index)
        finally:
            lock.release()
        return None
