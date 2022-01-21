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
def getBackbones(index_to_hmm, unaligned, workdir, backbone_dir):
    # return is a map from HMM index to a list of fragments that have high
    # bit scores at the HMM.
    ret = defaultdict(list)
    #if Configs.use_weight:
    #    weights = Weights.weights 
    #else:
    #    weights = Weights.ranked_bitscores

    # for each taxon and its scores on HMMs, sort it in decreasing order
    #for taxon, sorted_weights in weights.items():
    ret_str = ''
    weights_map = {}
    for taxon, seq in unaligned.items():
        if ret_str != '':
            ret_str += '\n'

        # load weights from local
        if Configs.use_weight:
            sorted_weights, this_weights_map = readWeights(taxon)
            weights_map[taxon] = this_weights_map
        else:
            sorted_weights = readBitscores(taxon)

        if sorted_weights == None:
            return 'N/A'
        #if not taxon in weights:
        #    return 'N/A'
        #sorted_weights = weights[taxon]

        top_k_hmms = sorted_weights[:Configs.num_hmms]
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
        #Configs.log('weights for {}: {}'.format(taxon, top_k_hmms))
        ret_str += 'weights for {}: {}'.format(taxon, top_k_hmms)
        for item in top_k_hmms:
            # append taxon to corresponding HMM i (defined by the index)
            ret[item[0]].append(taxon)

    # now we can split the fragments into chunks, each chunk corresponds to
    # all fragments added to an HMM i
    # For each HMM, run HMMAlign on its fragment chunk to get a backbone
    # alignment
    weights_path = workdir + '/weights.txt'
    weights_file = open(weights_path, 'w')
    
    # initialize the backbone index at backbone_start_index
    bb_index = 0
    for i, names in ret.items():
        #Configs.log("Generating fragment chunks/alignment for HMM {}".format(i))
        hmm_dir = workdir
        if not os.path.isdir(hmm_dir):
            os.makedirs(hmm_dir)
        if not os.path.isdir(hmm_dir + '/fragments'):
            os.makedirs('{}/fragments'.format(hmm_dir))

        # [NEW] save each single sequence to a fasta
        this_hmm = index_to_hmm[i].hmm_model_path
        for taxon in names:
            frag_path = '{}/fragments/{}.fasta'.format(hmm_dir, taxon)
            frag = unaligned.sub_alignment([taxon])
            frag.write(frag_path, 'FASTA')
        
            # hmmalign
            hmmalign_result_path = '{}/hmmalign.results.{}.{}.out'.format(
                    hmm_dir, taxon, i)
            cmd = '{} -o {} {} {}'.format(Configs.hmmalign_path,
                    hmmalign_result_path, this_hmm, frag_path)
            #if not (os.path.exists(hmmalign_result_path) and os.path.getsize(
            #    hmmalign_result_path) > 0):
            #    os.system(cmd)
            subprocess.run(cmd.split(' '))
        
            # Extended alignment
            ap_aln = ExtendedAlignment(frag.get_sequence_names())
            ap_aln.build_extended_alignment(index_to_hmm[i].alignment,
                    [hmmalign_result_path], True)
            #ap_aln.relabel_original_columns(remaining_cols)
            for key in ap_aln.keys():
                ap_aln[key] = ap_aln[key].upper()
            # save extended alignment to backbone directory
            this_bb_path = backbone_dir + '/{}_{}.fasta'.format(taxon, i)
            ap_aln.write(this_bb_path, 'FASTA')

            # add the weights of the bb to weights.txt
            if Configs.use_weight:
                real_this_bb_path = os.popen('realpath -s {}'.format(
                    this_bb_path)).read().split('\n')[0]
                #weights_file.write('{},{}\n'.format(
                #        real_this_bb_path, Weights.weights_map[taxon][i]))
                weights_file.write('{},{}\n'.format(
                        real_this_bb_path, weights_map[taxon][i]))
    weights_file.close()
    
    return ret_str

'''
The core function of GCM+eHMMs. Merge the HMMAlign results of the set of
queries to the input alignment.

The main differences between GCM+eHMMs and UPP are: 1) GCM+eHMMs uses adjusted
bitscore, which takes the number of queries in an HMM into consideration. 2)
GCM+eHMMs can utilize more than one HMM (while UPP uses the best HMM based on
bitscore) to align the queries; hence, more information is used.
'''
def alignSubQueries(backbone_path, index_to_hmm, lock, index):
    #global lock

    s11 = time.time()
    # add constraints
    constraints_dir = Configs.outdir + '/constraints/{}'.format(index)
    if not os.path.isdir(constraints_dir):
        os.makedirs(constraints_dir)
    # 0-th constraint set comes from the input alignment
    shutil.copyfile(backbone_path, '{}/c0.fasta'.format(
        constraints_dir))

    unaligned_dir = Configs.outdir + '/data'
    unaligned_path = unaligned_dir + '/unaligned_frag_{}.txt'.format(index)
    unaligned = Alignment(); unaligned.read_file_object(unaligned_path)
    unaligned_names = unaligned.get_sequence_names()

    c_index = 1
    for name in unaligned_names:
        single_frag = unaligned.sub_alignment([name])
        single_frag.write('{}/c{}.fasta'.format(constraints_dir, c_index),
                'FASTA')
        c_index += 1
    time_obtain_constraints = time.time() - s11

    s12 = time.time()
    # backbone alignments for GCM
    bb_dir = Configs.outdir + '/backbone_alignments/{}'.format(index)
    if not os.path.isdir(bb_dir):
        os.makedirs(bb_dir)
    # hmmsearch directory
    hmmsearch_dir = Configs.outdir + '/search_results/{}'.format(index)
    if not os.path.isdir(hmmsearch_dir):
        os.makedirs(hmmsearch_dir)

    # get backbones with the information we have
    weights_str = getBackbones(index_to_hmm, unaligned, hmmsearch_dir, bb_dir)
    time_obtain_backbones = time.time() - s12

    s13 = time.time()

    # time out the process if exceeding [timeout] seconds
    task_timed_out = False
    gcm_outdir = Configs.outdir + '/magus_outputs/{}'.format(index)

    if weights_str != 'N/A':
        # run GCM (modified MAGUS which takes in weights) on the subset
        est_path = Configs.outdir + '/temp/magus_result_{}.txt'.format(index)
        cmd = ['python3', Configs.magus_path, '-np', '1',
                '-d', gcm_outdir, '-s', constraints_dir, '-b', bb_dir,
                '-o', est_path, '-f', str(Configs.inflation_factor),
                '--graphclustermethod', Configs.graphclustermethod,
                '--graphtracemethod', Configs.graphtracemethod,
                '--graphtraceoptimize', Configs.graphtraceoptimize]
        if Configs.use_weight:
            weights_path = hmmsearch_dir + '/weights.txt'
            cmd += ['-w', weights_path] 
        #os.system(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL)

        try:
            p.wait(Configs.timeout)
        except subprocess.TimeoutExpired:
            task_timed_out = True
            parent = psutil.Process(p.pid)
            for c in parent.children(recursive=True):
                c.kill()
            parent.kill()
    else:
        est_path = 'skipped'

    # remove temp folders
    if not Configs.keeptemp:
        # trying shutil.rmtree to remove the directory
        shutil.rmtree(hmmsearch_dir)
        shutil.rmtree(bb_dir)
        shutil.rmtree(constraints_dir)
        if not Configs.keepgcmtemp and os.path.isdir(gcm_outdir):
            shutil.rmtree(gcm_outdir)
    time_gcm = time.time() - s13
    
    if not task_timed_out:
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
            if weights_str != 'N/A':
                Configs.debug("[MAGUS] Command used: {}".format(' '.join(cmd)))
                Configs.log(weights_str)
                Configs.log('{} passed to main pipeline...'.format(est_path))
            else:
                Configs.warning('Some taxa in task #{}'.format(index) \
                        + ' do not have any matching HMMs, skipping...')
        finally:
            lock.release()
        return est_path
    else:
        lock.acquire()
        try:
            Configs.warning('Task #{} failed, queuing for retry...'.format(
                index))
            alignSubQueries.q.put(index)
        finally:
            lock.release()
        return None
