'''
Created on 10.28.2021 by Chengze Shen

Merger of all query alignments to form the final alignment. The merging step
is exactly the same one from PASTA and UPP (by transitivity).
'''

import os, sys, re
import time
from configs import Configs
from helpers.alignment_tools import Alignment, read_fasta, \
        CompactAlignment, compact, ExtendedAlignment 
from functools import partial
#from concurrent.futures.process import ProcessPoolExecutor
from math import ceil

'''
function to merge a set of input paths (to alignments) sequentially
'''
def sequential_merger(inpaths):
    init_index = 0
    while init_index < len(inpaths) and inpaths[init_index] == 'skipped':
        init_index += 1
    init_aln = Alignment(); init_aln.read_file_object(inpaths[init_index])
    new_aln = compact(init_aln)
    for i in range(init_index + 1, len(inpaths)):
        inpath = inpaths[i]
        
        # skip these ones
        if inpath == 'skipped':
            continue
        frag_aln = Alignment(); frag_aln.read_file_object(inpath)
        new_aln.merge_in(compact(frag_aln))
    return new_aln

'''
function to take in a set of result paths for merging, and write
the merged alignment to an output path
'''
def mergeAlignments(inpaths, renamed_taxa, pool):
    Configs.log('(Naive merging) Merging all GCM subproblems with transitivity...')
    start = time.time()
    outpath = Configs.output_path
    assert len(inpaths) > 0

    # split paths into NUM_CPUS chunks
    chunks = []
    chunk_size = ceil(len(inpaths) / Configs.num_cpus)
    for i in range(0, len(inpaths), chunk_size):
        chunks.append(inpaths[i:min(i+chunk_size, len(inpaths))])

    # initialize Pool for multiprocessing
    #pool = Pool(Configs.num_cpus)
    merged_alns = list(pool.map(sequential_merger, chunks))
    #pool.close()
    #pool.join()

    # for the merged chunks, merge them into one alignment 
    final_aln = merged_alns[0]
    for i in range(1, len(merged_alns)):
        final_aln.merge_in(merged_alns[i])
    
    # revert back names of renamed taxa
    for original_name, rename in renamed_taxa.items():
        if rename in final_aln:
            final_aln[original_name] = final_aln[rename]
            final_aln.pop(rename)
    final_aln.write(outpath, 'FASTA')
    end = time.time()

    Configs.log('Finished merging all GCM subproblems, output file: {}'.format(
        outpath))
    Configs.runtime('Time to merge all outputs (s): {}'.format(end - start))

'''
helper function for collapsed merging (multiprocessing)
'''
def getQueryAlignment(backbone_keys, path):
    query = ExtendedAlignment([])
    query_name, insertions = query.read_query_alignment(backbone_keys, path)
    
    only_query = ExtendedAlignment([])
    only_query._col_labels = query._get_col_labels()
    only_query[query_name] = query[query_name]
    return only_query

'''
function to merge all subalignments to one alignment and with all singletons
in queries collapsed (in lower cases). This is the same behavior as UPP.
An additional "masked" version of the final alignment with all lower cases
removed will also be written to disk.
'''
def mergeAlignmentsCollapsed(backbone_alignment_path, inpaths,
        renamed_taxa, pool):
    Configs.log('(UPP-style merging) Merging all GCM subproblems ' \
            'with transitivity and singletons from queries collapsed...')
    start = time.time()
    outpath = Configs.output_path
    masked_outpath = Configs.output_path + '.masked'
    assert len(inpaths) > 0

    # read in all backbone sequences/alignment
    full_aln = ExtendedAlignment([])
    full_aln.read_file_object(backbone_alignment_path)
    full_aln.from_string_to_bytearray()
    backbone_keys = {x: 1 for x in full_aln.keys()}
    
    # read in queries so that insertions are marked
    func = partial(getQueryAlignment, backbone_keys)
    queries = list(pool.map(func, inpaths))

    # merge all queries to the backbone
    for query in queries:
        full_aln.merge_in(query, False)
        del query
    full_aln.from_bytearray_to_string()
    full_aln.write(outpath, 'FASTA')
    Configs.log('Finished merging all GCM subproblems, output file: {}'.format(
        outpath))
    
    # write a masked version of full alignment
    full_aln.remove_insertion_columns()
    full_aln.write(masked_outpath, 'FASTA')
    Configs.log('Masked final alignment written to : {}'.format(masked_outpath))

    end = time.time()
    Configs.runtime('Time to merge all outputs (s): {}'.format(end - start))
