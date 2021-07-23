import os, math 
from configs import Configs
from gcmm.loader import loadSubQueries 
from gcmm.weighting import loadWeights, Weights
from gcmm.aligner import alignSubQueries
from gcmm.merger import mergeAlignments
from helpers.alignment_tools import Alignment

'''
Delete all unnecessary intermediate files
'''
def clearTempFiles():
    os.system('rm {}/magus_result_*'.format(Configs.outdir))
    if os.path.isdir('{}/backbone_alignments'.format(Configs.outdir)):
        os.system('rm -r {}/backbone_alignments'.format(Configs.outdir))
    if os.path.isdir('{}/constraints'.format(Configs.outdir)):
        os.system('rm -r {}/constraints'.format(Configs.outdir))
    if os.path.isdir('{}/search_results'.format(Configs.outdir)):
        os.system('rm -r {}/search_results'.format(Configs.outdir))
    if not Configs.keepgcmtemp \
            and os.path.isdir('{}/magus_outputs'.format(Configs.outdir)):
        os.system('rm -r {}/magus_outputs'.format(Configs.outdir))

'''
Main process for GCM+eHMMs
'''
def mainAlignmentProcess():
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
    for i in range(num_subset):
        output_path = alignSubQueries(i, index_to_hmm)
        sub_alignment_paths.append(os.path.abspath(output_path))

    # 4) merge all results 
    print("\nAll GCM subproblems finished! Doing merging with transitivity...")
    mergeAlignments(sub_alignment_paths)

    if not Configs.keeptemp:
        clearTempFiles()
