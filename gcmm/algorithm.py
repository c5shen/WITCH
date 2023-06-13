'''
Created on 1.22.2022 by Chengze Shen

Algorithms for tree decomposition and hmmsearch.
'''

import os, subprocess, math, time, psutil, shutil 
from configs import Configs
from helpers.alignment_tools import Alignment, MutableAlignment 
from gcmm.tree import PhylogeneticTree
import dendropy
from dendropy.datamodel.treemodel import Tree
import tempfile
from functools import partial
import re

from helpers.math_utils import lcm
#from multiprocessing import Queue, Lock

'''
***** ADOPTED from sepp/exhaustive.py - ExhaustiveAlgorithm class *****
Class to perform backbone tree decomposition as does in UPP
'''
class DecompositionAlgorithm(object):
    def __init__(self, backbone_path, backbone_tree_path, alignment_size=10):
        self.symfrac = 0.0
        self.ere = 0.59
        self.informat = 'afa'
        self.molecule = Configs.molecule
        self.path = Configs.hmmbuildpath

        self.strategy = 'centroid'              # default in SEPP/UPP
        self.decomp_strategy = 'hierarchical'   # ensemble of HMMs
        self.alignment_size = alignment_size    # default: 10 as UPP
        self.minsubsetsize = 2
        self.pdistance = 1                      # default in SEPP/UPP
        self.distances = {}
        self.maxDiam = None                     # default in SEPP/UPP

        self.backbone_path = backbone_path
        self.backbone_tree_path = backbone_tree_path

        self.outdir = '{}/tree_decomp'.format(Configs.outdir)
        #self.timeout = 30

    '''
    Read in the backbone alignment and tree
    '''
    def read_alignment_and_tree(self):
        Configs.log('Reading backbone alignment: {}'.format(
            self.backbone_path))
        assert os.path.exists(self.backbone_path), \
                'Backbone alignment [{}] does not exist!'.format(
                        self.backbone_path)
        alignment = Alignment()
        alignment.read_file_object(self.backbone_path)

        Configs.log('Reading backbone tree: {}'.format(
            self.backbone_tree_path))
        assert os.path.exists(self.backbone_tree_path), \
                'Backbone tree [{}] does not exist'.format(
                        self.backbone_tree_path)
        tree = PhylogeneticTree(dendropy.Tree.get_from_stream(
            open(self.backbone_tree_path, 'r'),
            schema='newick',
            preserve_underscores=True))
        
        return alignment, tree

    '''
    Do the decomposition on the backbone alignment/tree
    Take in a ProcessPoolExecutor for parallelism
    '''
    def decomposition(self, lock, pool):
        start = time.time()
        Configs.log('Started decomposing the backbone to eHMM')
        alignment, tree = self.read_alignment_and_tree()

        assert isinstance(alignment, Alignment)
        assert isinstance(tree, PhylogeneticTree)

        tree.get_tree().resolve_polytomies()
        
        # label edges with numbers so that we can assemble them back
        # at the end
        tree.label_edges()

        # decompose the tree into alignment subsets
        alignment_tree_map = PhylogeneticTree(
                Tree(tree.den_tree)).decompose_tree(
                        self.alignment_size,
                        strategy=self.strategy,
                        minSize=self.minsubsetsize,
                        tree_map={},
                        decomp_strategy=self.decomp_strategy,
                        pdistance=self.pdistance,
                        distances=self.distances,
                        maxDiam=self.maxDiam)
        assert len(alignment_tree_map) > 0, (
                'Tree could not be decomposed '
                'given the following settings: '
                'strategy: {}\nminsubsetsize: {}\nalignment_size: {}'.format(
                    self.strategy, self.minsubsetsize, self.alignment_size))
        
        Configs.debug('Alignment subsets: {}'.format(len(alignment_tree_map)))

        subset_args = []
        for a_key, a_tree in alignment_tree_map.items():
            assert isinstance(a_tree, PhylogeneticTree)
            label = 'A_0_{}'.format(a_key)
            subset_taxa = a_tree.leaf_node_names()
            subset_args.append((label, subset_taxa))
        
        Configs.log('Creating an ensemble of HMMs: {} subsets'.format(
            len(subset_args)))

        # create all subset alignments and HMMBuild them
        outdirprefix = self.outdir + '/root'
        func = partial(subset_alignment_and_hmmbuild, lock, 
                self.path, outdirprefix,
                self.molecule, self.ere, self.symfrac,
                self.informat, self.backbone_path)
        ret = list(pool.map(func, subset_args))
        hmmbuild_paths = []
        
        # record the retained columns in each HMM subset
        # (used in mapping hmmalign results to alignment graph)
        subset_to_retained_columns = {}
        for item in ret:
            hmmbuild_paths.append(item[0])
            # item[1] -- A_0_*
            subset_to_retained_columns[int(item[1].split('_')[-1])] = item[2]
        
        assert len(hmmbuild_paths) == len(subset_args), \
                'Number of HMMs created does not match ' \
                'the original number of subsets'
        Configs.log('Finished creating {} HMMs to {}'.format(
            len(hmmbuild_paths), outdirprefix))
        
        dur = time.time() - start
        Configs.runtime('Time to decompose the backbone (s): {}'.format(
            dur))
        return hmmbuild_paths, subset_to_retained_columns

'''
Class to perform HMMSearch on all hmmbuild subsets and fragment sequences
(Also will break fragments to fragment chunks for better parallelism)
'''
class SearchAlgorithm(object):
    def __init__(self, hmmbuild_paths):
        self.hmmbuild_paths = hmmbuild_paths
        self.unaligned = None
        self.path = Configs.hmmsearchpath

        self.max_chunk_size = 20000             # default in SEPP/UPP

        self.filters = False                    # No filters
        self.elim = 99999999                    # elim value for hmmsearch
        self.piped = False                      # default in SEPP/UPP
        
        self.outdir = Configs.outdir + '/tree_decomp'

    # dummy constructor for pure usage of obtaining retained columns from 
    # user-provided hmm directory
    def __init__(self):
        pass

    ####### ONLY USED WHEN USER PROVIDES AN HMM DIRECTORY #######
    def readHMMDirectory(self, lock, pool):
        subset_to_retained_columns = dict()
        
        # use "find" command to find all subset directories
        cmd = 'find {} -maxdepth 2 -name A_0_* -type d'.format(Configs.hmmdir)
        subset_dirs = [os.path.realpath(x) 
                for x in os.popen(cmd).read().split('\n')[:-1]]
        Configs.log('User-provided HMM directory: {}'.format(Configs.hmmdir))
        Configs.log('Reading {} subsets...'.format(len(subset_dirs)))

        # use pool to process all subset dirs
        subset_args = []
        tmp_map = dict()
        for sdir in subset_dirs:
            ind = int(sdir.split('/')[-1].split('_')[-1])
            subset_args.append((ind, sdir))
            tmp_map[ind] = sdir

        # obtain the backbone alignment, which should be of A_0_0
        cmd = 'find {} -maxdepth 1 -name hmmbuild.input.*fasta -type f'.format(
                tmp_map[0])
        bb_path = os.path.realpath(os.popen(cmd).read().split('\n')[0])
        bb_aln = Alignment(); bb_aln.read_file_object(bb_path)
        
        func = partial(subset_obtain_retained_columns, bb_aln)
        ret = list(pool.map(func, subset_args))

        for item in ret:
            subset_to_retained_columns[int(item[0])] = item[1]

        del bb_aln
        Configs.log('Finished reading decomposition subsets...')
        return subset_to_retained_columns
    
    # main function to perform all-against-all query-HMM searches in MP 
    def search(self, lock, pool):
        Configs.log('Running all-against-all HMM searches between queries ' \
                'and HMMs...')
        start = time.time()
        # create fragment chunks and run HMMSearch on all subsets
        # **** USING SEPP exhaustive implementation to find the best number
        # **** of chunks (Least Common Multiple of #subsets & #cpus // #subsets)
        num_chunks = lcm(
                len(self.hmmbuild_paths),
                Configs.num_cpus) // len(self.hmmbuild_paths)
        #num_chunks = Configs.num_cpus
        frag_chunk_paths = self.read_and_divide_unaligned(num_chunks)

        # create subproblems of HMMSearch
        subset_args = []
        for hmmbuild_path in self.hmmbuild_paths:
            hmmsearch_outdir = '/'.join(hmmbuild_path.split('/')[:-1])
            hmm_label = hmmbuild_path.split('/')[-1].split('.')[-1]
            for i in range(0, len(frag_chunk_paths)):
                frag_chunk_path = frag_chunk_paths[i]
                subset_args.append((hmmsearch_outdir, hmm_label,
                    hmmbuild_path, frag_chunk_path, i))

        func = partial(subset_frag_chunk_hmmsearch, lock, self.path,
                self.piped, self.elim, self.filters)
        hmmsearch_paths = list(pool.map(func, subset_args))
        assert len(hmmsearch_paths) == len(subset_args), \
                'It seems that some HMMSearch jobs failed'
        Configs.log('Finished {} HMMSearch jobs.'.format(len(hmmsearch_paths)))

        dur = time.time() - start
        Configs.runtime('Time to run all-against-all HMMSearches (s): {}'.format(
            dur))
        return frag_chunk_paths 

    def read_and_divide_unaligned(self, num_chunks, extra_frags={}):
        max_chunk_size = self.max_chunk_size
        Configs.log('Reading in fragment files and breaking ' + 
                'them into at most {} chunks '.format(num_chunks) +
                '(each having at most {} sequences)'.format(max_chunk_size))

        # read in query paths (either provided or produced by separating
        # original input sequences to backbone and query sequences)
        unaligned = MutableAlignment()
        unaligned.read_file_object(Configs.query_path)
        ids_unaligned = unaligned.keys()

        # test if input fragment names contain whitespaces or tabs
        ids_unaligned_violations = [id_ for id_ in ids_unaligned
                if (' ' in id_) or ('\t' in id_)]
        if len(ids_unaligned_violations) > 0:
            raise ValueError(
                "Your input fragment file contains {} sequences, ".format(
                    len(ids_unaligned_violations)) + 
                "which names contain either whitespaces or tabs '\\t'. " +
                "Their names are:\n {}".format(
                    "'\n'  ".join(ids_unaligned_violations)))

        # add in extra fragments if provided
        for k, v in extra_frags.items():
            unaligned[k] = v.replace("-", "")

        alg_chunks = unaligned.divide_to_equal_chunks(num_chunks,
                max_chunk_size)
        
        # write to local (<WITCH_OUTPUT_DIR>/tree_decomp/fragment_chunks)
        # and save all frag chunk paths
        ret = []
        fc_outdir = self.outdir + '/fragment_chunks'
        if not os.path.isdir(fc_outdir):
            os.makedirs(fc_outdir)

        for i in range(0, len(alg_chunks)):
            temp_file = None
            if alg_chunks[i]:
                temp_file = '{}/fragment_chunk_{}.fasta'.format(fc_outdir, i)
                #temp_file = tempfile.mktemp(
                #        prefix='fragment_chunk_{}'.format(i),
                #        suffix='.fasta', dir=fc_outdir)
                alg_chunks[i].write(temp_file, 'FASTA')
                Configs.debug('Writing alignment chunk #{} to {}'.format(
                    i, temp_file))
                ret.append(temp_file)
        Configs.log("Finished breaking fragments into {} chunks to {}".format(
            len(ret), fc_outdir))
        return ret


################ HELPER FUNCTIONS FOR MULTI-PROCESSING ################

'''
Obtain subset alignment given taxa, and run hmmbuild on the subset
alignment.
'''
def subset_alignment_and_hmmbuild(lock, binary, outdirprefix, molecule, 
        ere, symfrac, informat, backbone_path, args):
    label, taxa = args
    alignment = Alignment(); alignment.read_file_object(backbone_path)
    subalignment = alignment.sub_alignment(taxa)
    retained_columns = subalignment.delete_all_gaps()
    del alignment

    outdir = os.path.join(outdirprefix, label)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    # write subalignment to outdir
    subalignment_path = '{}/hmmbuild.input.{}.fasta'.format(outdir, label)
    #subalignment_path = tempfile.mktemp(prefix='hmmbuild.input.',
    #        suffix='.fasta', dir=outdir)
    subalignment.write(subalignment_path, 'FASTA')

    # run HMMBuild with 1 cpu given the subalignment
    hmmbuild_path = '{}/hmmbuild.model.{}'.format(outdir, label)
    #hmmbuild_path = tempfile.mktemp(prefix='hmmbuild.model.',
    #        dir=outdir)
    cmd = [binary, '--cpu', '1',
            '--{}'.format(molecule),
            '--ere', str(ere),
            '--symfrac', str(symfrac),
            '--informat', informat,
            '-o', '/dev/null']
    cmd.extend([hmmbuild_path, subalignment_path])
    os.system(' '.join(cmd))
    lock.acquire()
    try:
        Configs.debug('[HMMBuild] Command used: {}'.format(' '.join(cmd)))
    finally:
        lock.release()
        return (hmmbuild_path, label, tuple(retained_columns))

'''
a single HMMSearch job between a frag chunk and an hmm
'''
def subset_frag_chunk_hmmsearch(lock, binary, piped, elim, filters, args):
    outdir, hmm_label, hmm, unaligned, frag_index = args
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    hmmsearch_path = '{}/hmmsearch.results.{}.fragment_chunk_{}'.format(
            outdir, hmm_label, frag_index)
    #hmmsearch_path = tempfile.mktemp(
    #        prefix='hmmsearch.results.fragment_chunk_{}.'.format(frag_index),
    #        dir=outdir)
    cmd = [binary, '--cpu', '1', '--noali', '-E', str(elim)]
    if not piped:
        cmd.extend(['-o', hmmsearch_path])
    if not filters:
        cmd.extend(['--max'])
    cmd.extend([hmm, unaligned])
    os.system(' '.join(cmd))

    # modify the output file and only retain taxon name, E-value, bit-score
    res = evalHMMSearchOutput(hmmsearch_path)
    with open(hmmsearch_path, 'w') as f:
        f.write(str(res))

    lock.acquire()
    try:
        Configs.debug('[HMMSearch] Command used: {}'.format(' '.join(cmd)))
    finally:
        lock.release()
        return hmmsearch_path

'''
helper function to read alignment from a subset alignment and obtain
its retained columns with respect to the backbone alignment
'''
def subset_obtain_retained_columns(bbaln, args):
    ind, sdir = args
    # find the subset alignment path
    cmd = 'find {} -maxdepth 1 -name hmmbuild.input.*fasta -type f'.format(
            sdir)
    try:
        path = os.path.realpath(os.popen(cmd).read().split('\n')[0])
        aln = Alignment(); aln.read_file_object(path)
        keys = aln.keys()
        subaln = bbaln.sub_alignment(keys)
        retained_columns = subaln.delete_all_gaps()
        return ind, retained_columns
    except:
        # cannot find the hmm input path in the corresponding directory
        raise Exception(
            'Cannot load subset {} at {}! Decomposition may be incomplete.'.format(
                ind, sdir))

'''
helper function for modifying HMMSearch output file
'''
def evalHMMSearchOutput(path):
    outfile = open(path, 'r')
    results = {}

    pattern = re.compile(
        r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+"
        r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)")
    start_reading = False
    for line in outfile:
        line = line.strip()
        if not start_reading and line.startswith("E-value") is True:
            start_reading = True
        elif start_reading and line == "":
            start_reading = False
            break
        elif start_reading:
            matches = pattern.search(line)
            if matches is not None and matches.group(0).find("--") == -1:
                results[matches.group(9).strip()] = (
                    float(matches.group(1).strip()),
                    float(matches.group(2).strip()))
                # _LOG.debug("Fragment scores;"
                #           "fragment:%s E-Value:%s BitScore:%s" %(matches
                # .group(9).strip(),matches.group(1).strip(), matches.
                # group(2).strip()))
    outfile.close()
    return results
