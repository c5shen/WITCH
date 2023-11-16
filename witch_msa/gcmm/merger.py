'''
Created on 10.28.2021 by Chengze Shen

Merger of all query alignments to form the final alignment. The merging step
is exactly the same one from PASTA and UPP (by transitivity).
'''

import os, sys, re
import time
from witch_msa.configs import Configs
from witch_msa.helpers.alignment_tools import Alignment, read_fasta, \
        CompactAlignment, compact, ExtendedAlignment 
from functools import partial
#from concurrent.futures.process import ProcessPoolExecutor
from math import ceil

'''
function to merge a set of input paths (to alignments) sequentially
'''
def sequential_merger(queries, inpaths):
    init_index = 0
    while init_index < len(inpaths) and inpaths[init_index] == 'skipped':
        init_index += 1
    init_aln = Alignment(); init_aln.read_file_object(queries[init_index])
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
function to merge all subalignments to one alignment and with all singletons
in queries collapsed (in lower cases). This is the same behavior as UPP.
An additional "masked" version of the final alignment with all lower cases
removed will also be written to disk.
'''
def mergeAlignmentsCollapsed(backbone_alignment_path, queries,
        renamed_taxa, pool):
    Configs.log('(UPP-style merging) Merging all GCM subproblems ' \
            'with transitivity and singletons from queries collapsed...')
    start = time.time()
    outpath = Configs.output_path
    masked_outpath = Configs.output_path + '.masked'
    if not (len(queries) > 0):
        print('No query alignment provided to merger!')
        exit(1)

    # read in all backbone sequences/alignment
    full_aln = ExtendedAlignment([])
    full_aln.read_file_object(backbone_alignment_path)
    full_aln.from_string_to_bytearray()
    
    # read in queries so that insertions are marked
    #backbone_keys = {x: 1 for x in full_aln.keys()}
    #func = partial(getQueryAlignment, backbone_keys)
    #queries = list(pool.map(func, inpaths))

    # merge all queries to the backbone
    for query in queries:
        if query != 'skipped':
            full_aln.merge_in(query, False)
            #del query
    full_aln.from_bytearray_to_string()
    
    # rename back taxa
    name_map = {v: k for k, v in renamed_taxa.items()}
    count = 0
    for name in list(name_map.keys()):
        ori_name = name_map[name]
        if name in full_aln:
            full_aln[ori_name] = full_aln[name]
            full_aln.pop(name)
            count += 1
    if count > 0:
        Configs.log('Converted {} names back to their originals'.format(count))
    Configs.log('Finished merging all GCM subproblems, output file: {}'.format(
        outpath))
    full_aln.write(outpath, 'FASTA')
    
    # write a masked version of full alignment
    full_aln.remove_insertion_columns()
    full_aln.write(masked_outpath, 'FASTA')
    Configs.log('Masked final alignment written to: {}'.format(masked_outpath))

    end = time.time()
    Configs.runtime('Time to merge all outputs (s): {}'.format(end - start))
