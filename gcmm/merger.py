import os, sys, re
import time
from configs import Configs
from helpers.alignment_tools import Alignment, read_fasta, \
        CompactAlignment, compact 

# function to take in a set of result paths for merging, and write
# the merged alignment to an output path
def mergeAlignments(inpaths):
    start = time.time()
    outpath = Configs.outdir + '/merged.fasta'
    assert len(inpaths) > 0
    init_aln = Alignment(); init_aln.read_file_object(inpaths[0])
    new_aln = compact(init_aln)

    for i in range(1, len(inpaths)):
        inpath = inpaths[i]
        Configs.debug('Merging {}'.format(inpath))
        frag_aln = Alignment(); frag_aln.read_file_object(inpath)
        new_aln.merge_in(compact(frag_aln))
    new_aln.write(outpath, 'FASTA')
    end = time.time()

    Configs.runtime('Time to merge all outputs (s): {}'.format(end - start))

#assert len(sys.argv) == 3, "need 2 fields: [result dir] [output path]"
#result_dir = sys.argv[1]
#out_path = sys.argv[2]
#
## read in all output files with prefix "magus_result_"
#cmd = 'find {} -name magus_result_* -type f'.format(result_dir)
#paths = os.popen(cmd).read().split('\n')[:-1]
#
## initial alignment 
#init_aln = Alignment(); init_aln.read_file_object(paths[0])
#
#new_aln = compact(init_aln)
#for i in range(1, len(paths)):
#    path = paths[i]
#    frag_aln = Alignment(); frag_aln.read_file_object(path)
#    new_aln.merge_in(compact(frag_aln))
#new_aln.write(out_path, 'FASTA')
