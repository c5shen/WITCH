'''
Created on 1.22.2022 by Chengze Shen

Algorithms for tree decomposition and hmmsearch.
'''

import re, os, subprocess, math, time, psutil, shutil 
from operator import add
from functools import partial
from tqdm import tqdm

import dendropy
from dendropy.datamodel.treemodel import Tree

from witch_msa.configs import Configs, tqdm_styles
from witch_msa.gcmm import *
from witch_msa.gcmm.tree import PhylogeneticTree
from witch_msa.helpers.alignment_tools import Alignment, MutableAlignment 
from witch_msa.helpers.math_utils import lcm
from witch_msa.gcmm.task import getTasks, runTasks

import concurrent.futures

#import pyhmmer
#from helpers.pyhmmer_tools import *


'''
***** ADOPTED from sepp/exhaustive.py - ExhaustiveAlgorithm class *****
Class to perform backbone tree decomposition as does in UPP
'''
class DecompositionAlgorithm(object):
    def __init__(self, backbone_path, backbone_tree_path,
            alignment_size=10, alignment_upper_bound=None):
        self.symfrac = 0.0
        self.ere = 0.59
        self.informat = 'afa'
        self.molecule = Configs.molecule
        self.path = Configs.hmmbuildpath

        self.strategy = 'centroid'              # default in SEPP/UPP
        self.decomp_strategy = 'hierarchical'   # ensemble of HMMs
        self.alignment_size = alignment_size    # default: 10 as UPP
        self.upper_bound = alignment_upper_bound# default: None as UPP
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
        Configs.log(' '.join(['The ensemble of HMM',
                'subsets have the following range (# sequences): [{}, {}]'.format(
                    self.alignment_size, self.upper_bound)]))
        for a_key, a_tree in alignment_tree_map.items():
            #assert isinstance(a_tree, PhylogeneticTree)
            label = 'A_0_{}'.format(a_key)
            subset_taxa = a_tree.leaf_node_names()
            # Added @ 1.14.2024 - Chengze Shen
            # if Configs.alignment_upper_bound is set, then use it to limit
            # the subset alignments based on their sizes (# of sequences)
            if self.upper_bound is not None:
                # ignore the subset if it has too many sequences
                if len(subset_taxa) > self.upper_bound:
                    continue
            subset_args.append((label, set(subset_taxa)))
        
        Configs.log('Creating an ensemble of HMMs: {} subsets'.format(
            len(subset_args)))

        # create all subset alignments and HMMBuild them
        if self.molecule is None:
            self.molecule = Configs.inferDataType(self.backbone_path)
        outdirprefix = self.outdir + '/root'
        func = partial(subset_alignment_and_hmmbuild, lock, 
                self.path, outdirprefix,
                self.molecule, self.ere, self.symfrac,
                self.informat, self.backbone_path)

        # Updated on 1.11.2024 - Chengze Shen
        #   - clearly this is an oversight when dealing with large backbone
        #   - and many subsets (creating too many HMMBuild jobs simultaneously
        #   - and encounters memory issue (since each job needs to read in
        #   - the backbone alignment
        #ret = list(pool.map(func, subset_args))
        mytasks = getTasks(subset_args)
        ret, _, _, _ = runTasks(func, pool, mytasks, len(subset_args),
                max_concurrent_jobs=Configs.max_concurrent_jobs)

        hmmbuild_paths = []
        
        # record the retained columns in each HMM subset
        # (used in mapping hmmalign results to alignment graph)
        subset_to_retained_columns = dict()
        subset_to_nongaps_per_column = dict()
        for item in ret:
            hmmbuild_paths.append(item[0])
            # item[1] -- A_0_*
            ind = int(item[1].split('_')[-1])
            subset_to_retained_columns[ind] = item[2]
            subset_to_nongaps_per_column[ind] = item[3]

        # convert to list
        # Modified @ 1.14.2024 - Chengze Shen
        # instead of using list, use the original dict with the label as key
        # this is to avoid idx out of bound issue if we set an upper bound
        # for subsets
        #subset_to_retained_columns = [subset_to_retained_columns[k]
        #            for k in sorted(subset_to_retained_columns.keys())]
        #subset_to_nongaps_per_column = [subset_to_nongaps_per_column[k]
        #            for k in sorted(subset_to_nongaps_per_column.keys())]
        assert len(hmmbuild_paths) == len(subset_args), \
                'Number of HMMs created does not match ' \
                'the original number of subsets'

        # Added @ 5.11.2024 - Chengze Shen
        # sanity check for whether all HMMs are created correctly
        probs = sanityCheckFileCreation(hmmbuild_paths)
        if len(probs) > 0:
            Configs.error('Some HMM files are not created correctly: {}'.format(
                ', '.join(probs)))
            notifyError(getLineInfo())
        else:
            Configs.log('Successfully creating {} HMMs at {}'.format(
                len(hmmbuild_paths), outdirprefix))
        
        dur = time.time() - start
        Configs.runtime(' '.join(['(DecompositionAlgorithm.decomposition)',
                'Time to decompose the backbone (s):', str(dur)]))
        return hmmbuild_paths, subset_to_retained_columns, \
                subset_to_nongaps_per_column

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
        self.molecule = Configs.molecule

        self.filters = False                    # No filters
        self.elim = 99999999                    # elim value for hmmsearch
        self.piped = False                      # default in SEPP/UPP
        
        self.outdir = Configs.outdir + '/tree_decomp'

    ####### ONLY USED WHEN PRESENTS WITH AN HMM DIRECTORY #######
    def readHMMDirectory(self, lock, pool):
        subset_to_retained_columns = dict()
        subset_to_nongaps_per_column = dict()
        
        # use "find" command to find all subset directories
        cmd = 'find {} -maxdepth 2 -name A_0_* -type d'.format(Configs.hmmdir)
        subset_dirs = [os.path.realpath(x) 
                for x in os.popen(cmd).read().split('\n')[:-1]]
        Configs.log('Found existing HMM directory: {}'.format(Configs.hmmdir))
        Configs.log('Reading {} subsets...'.format(len(subset_dirs)))
        
        # terminate if not finding any HMMs in the current directory
        if len(subset_dirs) == 0:
            msg = 'Cannot find any pre-existing HMMs in {}!'.format(
                    Configs.hmmdir) + ' Please remove the directory and rerun.'
            Configs.error(msg)
            raise FileNotFoundError(msg)
                    
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
            subset_to_nongaps_per_column[int(item[0])] = item[2]

        # convert to list
        # Modified @ 1.14.2024 - Chengze Shen
        # instead of using list, use the original dict with the label as key
        # this is to avoid idx out of bound issue if we set an upper bound
        # for subsets
        #subset_to_retained_columns = [subset_to_retained_columns[k]
        #            for k in sorted(subset_to_retained_columns.keys())]
        #subset_to_nongaps_per_column = [subset_to_nongaps_per_column[k]
        #            for k in sorted(subset_to_nongaps_per_column.keys())]

        #del bb_aln
        Configs.log('Finished reading decomposition subsets...')
        return bb_path, subset_to_retained_columns, subset_to_nongaps_per_column
    
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

        # run each subproblem
        if self.molecule is None:
            self.molecule = Configs.inferDataType(self.backbone_path)
        func = partial(subset_frag_chunk_hmmsearch, lock, self.path,
                self.piped, self.elim, self.filters, self.molecule)
        # Updated on 1.24.2024 - Chengze Shen
        #   - updated this part of code as well to make sure memory won't
        #   - be an issue (as the update around line 150)
        hmmsearch_paths, futures = [], []
        mytasks = getTasks(subset_args)
        ret, _, _, _ = runTasks(func, pool, mytasks, len(subset_args),
                max_concurrent_jobs=Configs.max_concurrent_jobs)
        for item in ret:
            if item:
                hmmsearch_paths.append(item)
        #for subset_arg in subset_args:
        #    futures.append(pool.submit(func, subset_arg))
        #for future in tqdm(
        #        concurrent.futures.as_completed(futures),
        #        total=len(subset_args), **tqdm_styles):
        #    res = future.result()
        #    if res:
        #        hmmsearch_paths.append(res)
        assert len(hmmsearch_paths) == len(subset_args), \
                'It seems that some HMMSearch jobs failed'
    
        # Added @ 5.11.2024 - Chengze Shen
        # sanity check for whether all HMMSearch jobs completed correctly
        probs = sanityCheckFileCreation(hmmsearch_paths)
        if len(probs) > 0:
            Configs.error('Some HMMSearch jobs are not completed correctly: {}'.format(
                ', '.join(probs)))
            notifyError(getLineInfo())
        else:
            Configs.log('Successfully completed {} HMMSearch jobs.'.format(
                len(hmmsearch_paths)))

        dur = time.time() - start
        Configs.runtime(' '.join(['(SearchAlgorithm.search) Time to run',
                'all-against-all HMMSearches (s):', str(dur)]))
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
    #alignment = Alignment(); alignment.read_file_object(backbone_path)
    #subalignment = alignment.sub_alignment(taxa)
    #retained_columns = subalignment.delete_all_gaps()

    # Updated @ 1.22.2024 - Chengze Shen
    # Changed the behavior for reading in backbone alignment from
    # "at once" to "line by line" (to avoid memory issues with large
    # backbone alignment)
    subalignment = Alignment()
    with open(backbone_path, 'r') as f:
        taxon, seq = f.readline(), ''
        while taxon:
            taxon = taxon.strip()
            if seq == '':
                seq = f.readline().strip()  # first entry
            # only retain taxon defined in taxa
            if taxon[1:] in taxa:
                subalignment[taxon[1:]] = seq
            # update to next entry
            taxon = f.readline()
            if taxon:
                seq = f.readline().strip()
    retained_columns = subalignment.delete_all_gaps()
    
    # also accumulate the number of nongaps for each column
    nongaps_per_column = [0] * subalignment.sequence_length() 
    for seq in subalignment.values():
        nongaps_per_column = map(add, nongaps_per_column,
                [int(x != '-') for x in seq])
    #del alignment

    outdir = os.path.join(outdirprefix, label)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    # write subalignment to outdir
    subalignment_path = '{}/hmmbuild.input.{}.fasta'.format(outdir, label)
    subalignment.write(subalignment_path, 'FASTA')

    ############# 9.19.2023 ############
    ## added version to use pyhmmer for building
    #hmmbuild_path = '{}/hmmbuild.model.{}'.format(outdir, label)
    #alphabet = moleculeToAlphabet(molecule)
    #msa_d = alignmentToTextMSA(
    #        subalignment, name=f'hmmbuild.model.{label}').digitize(alphabet)

    #builder = pyhmmer.plan7.Builder(alphabet, ere=ere, symfrac=symfrac)
    #background = pyhmmer.plan7.Background(alphabet)
    #hmm, _, _ = builder.build_msa(msa_d, background)
    #with open(hmmbuild_path, 'wb') as fh:
    #    hmm.write(fh, binary=False)
    #lock.acquire()
    #try:
    #    Configs.debug('[pyhmmer.plan7.Builder] ere={}, symfrac={}'.format(
    #        ere, symfrac))
    #finally:
    #    lock.release()
    #    return (hmmbuild_path, label, tuple(retained_columns),
    #            tuple(nongaps_per_column))

    # run HMMBuild with 1 cpu given the subalignment
    hmmbuild_path = '{}/hmmbuild.model.{}'.format(outdir, label)
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
        return (hmmbuild_path, label, tuple(retained_columns),
                tuple(nongaps_per_column))

'''
a single HMMSearch job between a frag chunk and an hmm
'''
def subset_frag_chunk_hmmsearch(lock, binary, piped, elim, filters,
        molecule, args):
    outdir, hmm_label, hmmbuild_path, unaligned, frag_index = args
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    ########### 9.19.2023 ##############
    ## use pyhmmer.plan7.Pipeline to do the HMMSearch
    #hmmsearch_path = '{}/hmmsearch.results.{}.fragment_chunk_{}'.format(
    #        outdir, hmm_label, frag_index)
    #alphabet = moleculeToAlphabet(molecule)
    #background = pyhmmer.plan7.Background(alphabet)

    ## load HMM
    #hmm = None
    #with pyhmmer.plan7.HMMFile(hmmbuild_path) as hmm_file:
    #    hmm = hmm_file.read()
    #
    ## set up the pipeline
    #kwargs = {'E': elim}
    #if not filters:
    #    kwargs.update({'bias_filter': False, 'F1': 1.0, 'F2': 1.0, 'F3': 1.0})
    #pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background,
    #        **kwargs)

    ## run the pipeline
    #tophits = None
    #with pyhmmer.easel.SequenceFile(unaligned,
    #        digital=True, alphabet=alphabet) as seq_file:
    #    tophits = pipeline.search_hmm(hmm, seq_file)
    #res = evalHMMSearchOutputPyhmmer(tophits)
    #with open(hmmsearch_path, 'w') as f:
    #    f.write(str(res))
    #lock.acquire()
    #try:
    #    Configs.debug('[pyhmmer.plan7.Pipeline.search_hmm] query={}, model={}, outpath={}, kwargs={}'.format(
    #        unaligned, hmmbuild_path, hmmsearch_path, kwargs))
    #finally:
    #    lock.release()
    #    return hmmsearch_path

    hmmsearch_path = '{}/hmmsearch.results.{}.fragment_chunk_{}'.format(
            outdir, hmm_label, frag_index)
    cmd = [binary, '--cpu', '1', '--noali', '-E', str(elim)]
    if not piped:
        cmd.extend(['-o', hmmsearch_path])
    if not filters:
        cmd.extend(['--max'])
    cmd.extend([hmmbuild_path, unaligned])
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
Additionally, return the number of non-gap characters for each retained column
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
        
        # accumulate the number of nongaps for each column
        nongaps_per_column = [0] * subaln.sequence_length() 
        for seq in subaln.values():
            nongaps_per_column = map(add, nongaps_per_column,
                    [int(x != '-') for x in seq])

        return ind, tuple(retained_columns), tuple(nongaps_per_column)
    except FileNotFoundError as e:
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

'''
helper function to keep just taxon name, e-value, and bit-score from tophits
'''
def evalHMMSearchOutputPyhmmer(tophits):
    results = {}
    for hit in tophits:
        name, evalue, score = hit.name.decode(), hit.evalue, hit.score
        results[name] = (evalue, f'{score:.1f}')
    return results
