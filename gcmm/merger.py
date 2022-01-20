'''
Created on 10.28.2021 by Chengze Shen

Merger of all query alignments to form the final alignment. The merging step
is exactly the same one from PASTA and UPP (by transitivity).
'''

import os, sys, re
import time
from configs import Configs
from helpers.alignment_tools import Alignment, read_fasta, \
        CompactAlignment, compact 
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
def mergeAlignments(inpaths, pool):
    Configs.log('Merging all GCM subproblems with transitivity...')
    start = time.time()
    outpath = Configs.outdir + '/merged.fasta'
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
    final_aln.write(outpath, 'FASTA')
    end = time.time()

    Configs.log('Finished merging all GCM subproblems, output file: {}'.format(
        outpath))
    Configs.runtime('Time to merge all outputs (s): {}'.format(end - start))
