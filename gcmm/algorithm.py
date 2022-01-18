import os, subprocess, math, psutil, shutil 
from configs import Configs
from helpers.alignment_tools import Alignment
from gcmm.tree import PhylogeneticTree
import dendropy
from dendropy.datamodel.treemodel import Tree
import tempfile
from functools import partial

from multiprocessing import Queue, Lock

'''
***** ADOPTED from sepp/exhaustive.py - ExhaustiveAlgorithm class *****
Class to perform backbone tree decomposition as does in UPP
'''
class DecompositionAlgorithm(object):
    def __init__(self):
        self.elim = 99999999    # elim value for hmmsearch
        self.filters = False
        self.symfrac = True
        self.informat = 'fasta'
        self.molecule = Configs.molecule
        self.path = Configs.hmmbuild_path

        self.strategy = 'centroid'              # default in SEPP/UPP
        self.decomp_strategy = 'hierarchical'   # ensemble of HMMs
        self.alignment_size = 10                # default in UPP
        self.minsubsetsize = 1
        self.pdistance = 1                      # default in SEPP/UPP
        self.distances = {}
        self.maxDiam = None                     # default in SEPP/UPP

        self.backbone_path = Configs.backbone_path
        self.backbone_tree_path = Configs.backbone_tree_path

        self.outdir = '{}/tree_decomp/root'.format(Configs.outdir)
        self.timeout = 30

    '''
    Read in the backbone alignment and tree
    '''
    def read_alignment_and_tree(self):
        Configs.log('Reading backbone alignment: {}'.format(
            self.backbone_path))
        alignment = Alignment()
        alignment.read_file_object(Configs.backbone_path)

        Configs.log('Reading backbone tree: {}'.format(
            Configs.backbone_tree_path))
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
        func = partial(subset_alignment_and_hmmbuild, lock, 
                self.path, self.outdir,
                self.molecule, alignment)
        results = list(pool.map(func, subset_args))
        for r in results:
            pass
        
'''
Obtain subset alignment given taxa, and run hmmbuild on the subset
alignment.
'''
def subset_alignment_and_hmmbuild(lock, binary, outdirprefix, molecule, 
        alignment, args):
    label, taxa = args
    subalignment = alignment.sub_alignment(taxa)
    subalignment.delete_all_gaps()

    outdir = os.path.join(outdirprefix, label)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    # write subalignment to outdir
    subalignment_path = tempfile.mktemp(prefix='hmmbuild.input.',
            suffix='.fasta', dir=outdir)
    subalignment.write(subalignment_path, 'FASTA')

    # run HMMBuild with 1 cpu given the subalignment
    hmmbuild_path = tempfile.mktemp(prefix='hmmbuild.model.',
            dir=outdir)
    cmd = [binary, '--cpu', '1', '--{}'.format(molecule),
            '--ere', '0.59', '--symfrac', '0.0', '--informat', 'afa',
            '-o', '/dev/null']
    cmd.extend([hmmbuild_path, subalignment_path])
    os.system(' '.join(cmd))
    lock.acquire()
    try:
        Configs.debug('Command used: {}'.format(' '.join(cmd)))
    finally:
        lock.release()
