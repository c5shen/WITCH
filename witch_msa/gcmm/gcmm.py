'''
Created on 1.19.2022 by Chengze Shen

Main pipeline of WITCH
'''

import os, sys, math, psutil, shutil 
from witch_msa.configs import * 
from witch_msa.gcmm.algorithm import DecompositionAlgorithm, SearchAlgorithm
from witch_msa.gcmm.loader import loadSubQueries, rankBitscores, writeTempBackbone
from witch_msa.gcmm.weighting import writeWeights, writeBitscores, \
        readWeightsFromLocal, writeWeightsToLocal
from witch_msa.gcmm.aligner import alignSubQueries, alignSubQueriesNew
from witch_msa.gcmm.backbone import BackboneJob
from witch_msa.gcmm.merger import mergeAlignmentsCollapsed
from witch_msa.gcmm.results_handler import submitAndCollectFutures

# experimental, customized ProcessPoolExecutor
from witch_msa.gcmm import *

from witch_msa.helpers.alignment_tools import Alignment
from witch_msa.helpers.general_tools import memoryUsage

import multiprocessing as mp
import concurrent.futures
from multiprocessing import Lock, Queue, Manager#, Pool
from concurrent.futures.process import ProcessPoolExecutor
from functools import partial
from tqdm import tqdm

# max system recursion limit hard encoding to a large number
# a temp fix for dendropy tree recursion issues
sys.setrecursionlimit(30000)

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
def initiate_pool(parser, cmdline_args):
    buildConfigs(parser, cmdline_args)

'''
6.21.2023 - Additionally, init pool for query alignment with two dicts
          - used later
'''
def initiate_pool_query_alignment(q, parser, cmdline_args,
        subset_to_retained_columns,
        subset_to_nongaps_per_column):
    buildConfigs(parser, cmdline_args)
    alignSubQueries.q = q
    alignSubQueriesNew.subset_to_retained_columns = subset_to_retained_columns
    alignSubQueriesNew.subset_to_nongaps_per_column = \
            subset_to_nongaps_per_column

'''
Main process for WITCH 
'''
def mainAlignmentProcess(parser, cmdline_args):
    m = Manager()
    lock = m.Lock()
    #l = Lock()
    q = Queue()

    # initialize the main pool at the start so that it's not memory
    # intensive
    Configs.warning('Initializing ProcessorPoolExecutor instance...')
    pool = ProcessPoolExecutor(Configs.num_cpus,
            initializer=initiate_pool, initargs=(parser, cmdline_args,))
            #mp_context=mp.get_context('spawn'),

    # if not user provided, default to <outdir>/tree_decomp/root
    #if not Configs.hmmdir:
    #    Configs.hmmdir = Configs.outdir + '/tree_decomp/root'

    # 0) obtain the backbone alignment/tree and eHMMs
    # If no UPP eHMM directory detected, decompose from the backbone
    if not Configs.hmmdir:
        # default to <outdir>/tree_decomp/root
        Configs.hmmdir = Configs.outdir + '/tree_decomp/root'
    else:
        # if user provides a directory but it does not exist, exit program
        assert os.path.isdir(Configs.hmmdir), \
            'Provided HMM directory does not exist'

    if not os.path.isdir(Configs.hmmdir):
        # if both backbone alignment/tree are provided by the user
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
            bb_job = BackboneJob(Configs.backbone_path, Configs.query_path,
                    Configs.backbone_tree_path)
            bb_job.setup()

            Configs.backbone_path, Configs.query_path = \
                    bb_job.run_alignment()
            Configs.backbone_tree_path = bb_job.run_tree()
    
        # after obtaining backbone alignment/tree, perform decomposition
        # and HMMSearches
        print('\nDecomposing the backbone tree...')
        decomp = DecompositionAlgorithm(
                Configs.backbone_path, Configs.backbone_tree_path,
                Configs.alignment_size, Configs.alignment_upper_bound)
        hmmbuild_paths, subset_to_retained_columns, subset_to_nongaps_per_column = \
                decomp.decomposition(lock, pool)
        print('\nPerforming all-against-all HMMSearches ' \
                'between the backbone and queries...')
        search = SearchAlgorithm(hmmbuild_paths)
        hmmsearch_paths = search.search(lock, pool)
        del decomp; del search
    else:
        # go over the given hmm directory and obtain all subset alignment
        # get their retained columns with respect to the backbone alignment
        print('\nFound existing HMM directory: {}'.format(Configs.hmmdir))        
        _dummy_search = SearchAlgorithm(None)
        backbone_path, subset_to_retained_columns, subset_to_nongaps_per_column \
                = _dummy_search.readHMMDirectory(lock, pool)
        if not Configs.backbone_path:
            Configs.backbone_path = backbone_path

    # sanity check before moving on
    # minimum required at this point:
    #   1. hmmdir   2. backbone_path    3. query_path
    assert Configs.backbone_path != None \
            and os.path.exists(Configs.backbone_path), 'backbone alignment missing'
    #assert Configs.backbone_tree_path != None \
    #        and os.path.exists(Configs.backbone_tree_path), 'backbone tree missing'
    assert Configs.query_path != None \
            and os.path.exists(Configs.query_path), 'query sequences missing'
    assert Configs.hmmdir != None \
            and os.path.isdir(Configs.hmmdir), 'eHMM directory missing'

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
            initargs=(q, parser, cmdline_args, subset_to_retained_columns,
                subset_to_nongaps_per_column))
            #mp_context=mp.get_context('spawn'),
    Configs.warning('ProcessPoolExecutor instance re-opened (for query alignment).')

    # 0) create a temporary (working) backbone alignment,
    # enforcing all letters to be upper-cases
    tmp_backbone_path, backbone_length = writeTempBackbone(
            Configs.outdir + '/tree_decomp/backbone', Configs.backbone_path)

    # 1) get all sub-queries
    num_seq, index_to_hmm, sid_to_query_names, \
            sid_to_query_seqs, renamed_taxa = loadSubQueries(lock, pool)

    # 2) calculate weights, if needed
    # read in existing weights.txt file if it exists
    weight_path = Configs.outdir + '/weights.txt'
    if os.path.exists(weight_path): 
        print('\nFound existing weights: {}'.format(weight_path))
        taxon_to_weights = readWeightsFromLocal(weight_path)
    else:
        ranked_bitscores = rankBitscores(index_to_hmm, renamed_taxa, lock, pool)
        if Configs.use_weight:
            print('\nCalculating weights (adjusted bit-scores)...')
            taxon_to_weights = writeWeights(index_to_hmm, ranked_bitscores, pool)
        else:
            print('\nLoading bit-scores...')
            taxon_to_weights = writeBitscores(ranked_bitscores, pool)
        del ranked_bitscores

        # 2a) if saving weights to local
        if Configs.save_weight:
            print('\n(user option) Writing weights to local...')
            writeWeightsToLocal(taxon_to_weights, Configs.outdir + '/weights.txt')

    # 3) solve each subset
    sub_alignment_paths = []
    if not os.path.isdir(Configs.outdir + '/temp'):
        os.makedirs(Configs.outdir + '/temp')

    # ProcessPoolExecutor version
    print('\nPerforming GCM alignments on query subsets...')

    queries, ignored_indexes = submitAndCollectFutures(pool, lock,
            sid_to_query_names, sid_to_query_seqs,
            num_seq, index_to_hmm, taxon_to_weights,
            tmp_backbone_path, backbone_length)

    # 4) merge all results 
    print('\nAll GCM subproblems finished! Doing merging with transitivity...')
    mergeAlignmentsCollapsed(tmp_backbone_path, queries, renamed_taxa, pool)

    # if there are any ignored indexes (queries), write them to local
    ignored_path = Configs.outdir + '/ignored_queries.fasta'
    if len(ignored_indexes) > 0:
        Configs.log('Writing {} ignored sequences to local at {}'.format(
            len(ignored_indexes), ignored_path)) 
        newname_to_oldname = {v: k for k, v in renamed_taxa.items()}
        with open(ignored_path, 'w') as f:
            for _i in ignored_indexes:
                name, seq = sid_to_query_names[_i], sid_to_query_seqs[_i]
                if name in newname_to_oldname:
                    name = newname_to_oldname[name]
                f.write('>{}\n{}\n'.format(name, seq))

    Configs.warning('Closing ProcessPoolExecutor instance...')
    pool.shutdown()
    Configs.warning('ProcessPoolExecutor instance closed.')
    if not Configs.keeptemp:
        clearTempFiles()
