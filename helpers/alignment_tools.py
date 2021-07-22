import time
import numpy as np
import re, os, sys
import itertools
from copy import deepcopy
from random import random
from functools import reduce
from dendropy.datamodel.taxonmodel import Taxon
import time
import io
from collections import defaultdict, Mapping

try:
    filetypes = (io.IOBase, file)
except NameError:
    filetypes = io.IOBase

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
    from io import BytesIO

gappat = re.compile(r'[-?]+')
nogappat = re.compile(r'[^-?]+')
_INDEL = re.compile(r"[-]")
_COACH_INS = re.compile(r"[.]")
_DANGEROUS_NAME_CHARS = re.compile(r"[^a-zA-Z0-9]")
ILLEGAL_CHARS = re.compile(r"[^a-zA-Z?-]")

# helper function, prints out all functionalities in this python file
def printHelper():
    print("1.\tdef readFastSP(path) --\tRead and parse FastSP output." +
            "\n\treturns a map of SPFN, SPFP, and Expansion.")
    print("2.\tdef readHMMSearch(path) --\tRead and parse HMMSearch output." +
            "\n\treturns a map of taxa to tuples: (e-value, bitscore)")
    print("3.\tdef read_fasta --\tRead in a fasta file, used in Alignment()")
    print("4.\tdef write_fasta --\tWrite an Alignment object to local.")
    print("5.\tclass Alignment() --\tAlignment class, for jobs with alignment")

# Function to read in a log file, which has runtime information attached
# at the bottom (using /usr/bin/time -v to generate)
# return: the runtime in minutes
def readRuntime(file_path):
    if os.path.isfile(file_path):
        runtime = 0
        with open(file_path, 'r') as f:
            lines = f.read().split('\n')
            if len(lines) < 23:
                return np.nan
            raw = lines[-24:-1][4].split('):')[1].strip().split(':')
            if len(raw) == 2:
                runtime = int(raw[0]) + float(raw[1])/60
            elif len(raw) == 3:
                runtime = int(raw[0]) * 60 + int(raw[1]) + float(raw[2])/60
            else:
                runtime = np.nan
        return runtime
    else:
        return np.nan

# Function to read in a fastsp file and parse the results
# return: a map of SPFN, SPFP, Expansion score
def readFastSP(file_path):
    ret = {'SPFN': np.nan, 'SPFP': np.nan, 'error': np.nan,
            'Expansion': np.nan}
    if os.path.isfile(file_path):
        f = open(file_path, 'r')
        lines = f.read().split('\n')[:-1]
    else:
        return ret
    if len(lines) > 0:
        ret['SP-Score'] = float(lines[0].split(' ')[1])
        ret['Modeler-Score'] = float(lines[1].split(' ')[1])
        ret['SPFN'] = float(lines[2].split(' ')[1])
        ret['SPFP'] = float(lines[3].split(' ')[1])
        ret['Expansion'] = float(lines[4].split(' ')[1])
        ret['TC'] = float(lines[5].split(' ')[1])
        ret['error'] = (ret['SPFN'] + ret['SPFP']) / 2
    return ret


################ Copied from PASTA-Alignment        ##################
################ CompactAlignment class definition  ##################
def write_compact_to_fasta(alignment, dest):
    """Writes the `alignment` in FASTA format to either a file object or file"""
    file_obj = None
    if isinstance(dest, str):
        file_obj = open(dest, "w")
    else:
        file_obj = dest
    for name in list(alignment.keys()):
        s = alignment.as_string_sequence(name)
        file_obj.write('>{}\n{}\n'.format(name, s))
    if isinstance(dest, str):
        file_obj.close()

def read_fasta(src):
    """generator that returns (name, sequence) tuples from either a FASTA
    formatted file or file object.
    """
    file_obj = None
    if isinstance(src, str):
        try:
            file_obj = open(src, "r")
        except IOError:
            print(("The file `%s` does not exist, exiting gracefully" % src))
    elif isinstance(src, filetypes):
            file_obj = src
    else:
        raise TypeError('FASTA reader cannot recognize the source of %s, %s, %s' % (src,type(src),isinstance(src, filetypes)))
    name = None
    seq_list = list()
    for line_number, i in enumerate(file_obj):
        if i.startswith('>'):
            if name:
                yield name, ''.join(seq_list)
                seq_list = list()
            name = i[1:].strip()
        else:
            #seq = ''.join(i.strip().upper().split())
            seq = ''.join(i.strip().split())
            #if not is_sequence_legal(seq):
            #    raise Exception("Error: illegal characeters in sequence at line %d" % line_number)
            seq_list.append(seq)
    if name:
        yield name, ''.join(seq_list)
    if isinstance(src, str):
        file_obj.close()

def write_fasta(alignment, dest):
    """Writes the `alignment` in FASTA format to either a file object or file"""
    file_obj = None
    if isinstance(dest, str):
        file_obj = open(dest, "w")
    else:
        file_obj = dest
    for name, seq in list(alignment.items()):
        file_obj.write('>%s\n%s\n' % (name, seq) )
    if isinstance(dest, str):
        file_obj.close()

class Alignment(dict, object):
    """A simple class that maps taxa names to sequences.
    TODO: switch to dendropy character_matrix
    """
    def __init__(self):
        "creates an empty matrix"
        dict.__init__(self)
        self.datatype = None

    def get_datatype(self):
        return self._datatype

    def set_datatype(self, d):
        if d is None:
            self._datatype = None
        else:
            self._datatype = d.upper()

    datatype = property(get_datatype, set_datatype)

    def get_sequence_names(self):
        "returns a list of sequence names"
        return list(self.keys())

    def get_num_taxa(self):
        "returns the number sequences"
        return len(self.get_sequence_names())

    def read_file_object(self, file_obj, file_format='FASTA'):
        """Augments the matrix by reading the file object.
        If duplicate sequence names are encountered then the old name will be replaced.
        """
        if ( file_format.upper() == 'FASTA' ):
            read_func = read_fasta
        #elif ( file_format.upper() == 'NEXUS' ):
        #    read_func = read_nexus
        #elif ( file_format.upper() == 'PHYLIP' ):
        #    read_func = read_phylip
        #elif ( file_format.upper() == 'COMPACT3' ):
        #    read_func = read_compact3
        else:
            raise NotImplementedError("Unknown file format (%s) is not supported" % file_format)
        for name, seq in read_func(file_obj):
            self[name] = seq

    def write_filepath(self, filename, file_format='FASTA', zipout=False):
        """Writes the sequence data in the specified `file_format` to `filename`"""
        
        file_obj = open_with_intermediates(filename,'w')
        if zipout:
            file_obj.close()
            file_obj = StringIO()
        self.write(file_obj, file_format=file_format)
        if zipout:
            import gzip
            file_obj_gz = gzip.open(filename, "wb", 6)
            file_obj_gz.write(str.encode(file_obj.getvalue()))
            file_obj_gz.close()
        file_obj.close()

    def write(self, file_obj, file_format):
        """Writes the sequence data in the specified `file_format` to `file_obj`"""
        if ( file_format.upper() == 'FASTA' ):
            write_func = write_fasta
        #elif ( file_format.upper() == 'NEXUS' ):
        #    write_func = write_nexus
        #elif ( file_format.upper() == 'PHYLIP' ):
        #    write_func = write_phylip
        #elif ( file_format.upper() == 'COMPACT' ):
        #    write_func = write_compact            
        #elif ( file_format.upper() == 'COMPACT2' ):
        #    write_func = write_compact2       
        #elif ( file_format.upper() == 'COMPACT3' ):
        #    write_func = write_compact3
        else:
            write_func = write_fasta
        write_func(self, file_obj)

    def unaligned(self):
        """
        Returns a new alignment with all gaps and missing sequences removed.
        """
        new_alignment = Alignment()
        new_alignment.datatype = self.datatype
        for name, seq in self.items():
            new_seq = re.sub(_INDEL, '', str(seq))
            if new_seq != '':
                new_alignment[name] = new_seq
        return new_alignment

    def sub_alignment(self, sub_keys):
        "Creates an new alignment with a subset of the taxa."
        new_alignment = Alignment()
        new_alignment.datatype = self.datatype
        for key in sub_keys:
            if key in self:
                new_alignment[key] = self[key]
        return new_alignment

    def is_empty(self):
        return self.__len__() < 1

    def is_aligned(self):
        if self.is_empty():
            raise ValueError("The alignment is empty.\n")
        else:
            v = list(self.values())
            first_seq_len = len(v[0])
            return all([len(i) == first_seq_len for i in v])

    def partition_info(self, base=0):
        return (self.datatype, 1+base, self.sequence_length() + base)

    def sequence_length(self):
        if self.is_aligned():
            return len(list(self.values())[0])

    def max_sequence_length(self):
        return max(len(re.sub(_INDEL, '', v)) for v in list(self.values()))

    def get_all_gap_cols(self):
        all_gaps = list(range(0, self.sequence_length()))
        for seq in self.values():
            all_gaps[:] = [i for i in all_gaps if seq[i] == '-']
        return all_gaps

    def remove_columns(self, indexes):
        for name, seq in self.items():
            self[name] = ''.join((char for idx, char in enumerate(seq)
                if idx not in indexes))

    def delete_all_gaps(self):
        rem = set(self.get_all_gap_cols())
        subset = [x for x in range(0, self.sequence_length()) if x not in rem]
        self.remove_columns(set(rem))
        return subset
    
    def from_bytearray_to_string(self):
        for k,v in self.items():
            self[k] = str(v)

    def from_string_to_bytearray(self):
        for k,v in self.items():
            self[k] = bytearray(v)   

    def mask_gapy_sites(self,minimum_seq_requirement):        
        n = len(list(self.values())[0])
        # The following implements column-based masking. Seems to be more efficient than row based
        masked = []
        allseqs = list(self.values())
        allseqs.sort(key=lambda x: x.count("-"))
        for c in range(0,n):
            r = minimum_seq_requirement
            for seq in allseqs:
                if seq[c] != "-":
                    r -= 1
                if r <= 0:
                    break
            if r > 0:
                masked.append(c)
                 
        #_LOG.debug("%d Columns identified for masking" %len(masked))
        if not masked:
            return
        included = [z for z in reduce(lambda x,y: x+[(x[-1][1]+1,y)],masked,[(-1,-1)]) if z[0]!=z[1]]
        if not included:
            included.append((masked[-1]+1,n))
        if included[-1][1] < n and masked[-1]+1 != n:
            included.append((masked[-1]+1,n))
        for k,seq in self.items():
            tmp = []
            for (i,j) in included:
                tmp.append(seq[i:j])
            self[k] = "".join(tmp)
        nn = len(list(self.values())[0])
        assert (len(masked) == n - nn), "Masking results is not making sense: %d %d %d" %(len(masked), n , nn)
        #_LOG.debug("Masking done. Before masking: %d; After masking: %d; minimum requirement: %d;" %(n,nn,minimum_seq_requirement))

    def merge_in(self, she):
        merge_in(self,she)


class AlignmentSequence:
    def __init__(self):
        self.seq = None
        self.pos = []
    
    def replace(self, match_char, replace_char):
        s = AlignmentSequence()
        s.pos.extend(self.pos)
        s.seq = self.seq.replace(match_char, replace_char)
        return s
    
    def as_string(self, pad_to):
        s=[]
        nxt = 0
        for i,p in enumerate(self.pos):
            if nxt != p:
                s.append("-" * (p-nxt))
            s.append(self.seq[i]) 
            nxt = p+1
        if (pad_to > nxt):
            s.append("-" * (pad_to-nxt))
        return "".join(s)
    
    def __str__(self):
        return self.as_string(0)
    
    def __repr__(self):
        return self.__str__()

class CompactAlignment(dict,object):
    def __init__(self):
        self.colcount = 0
        self.datatype = None
    
    def sub_alignment(self, sub_keys):
        "Creates an new alignment with a subset of the taxa."
        new_alignment = CompactAlignment()
        new_alignment.datatype = self.datatype
        for key in sub_keys:
            if key in self:
                new_alignment[key] = self[key]
        return new_alignment
   
    def is_aligned(self):
        return True
    def sequence_length(self):
        return self.colcount
    def get_num_taxa(self):
        return len(self)
    
    #def unaligned(self):
    #    n = Alignment()
    #    n.datatype = self.datatype
    #    for k,seq in self.items():
    #        n[k] = seq.seq
    #    return n
    
    def iter_column_character_count(self, seqsubset = None):
        if seqsubset is None:
            seqsubset = list(self.keys())
        
        counts = [0] * self.colcount
        for k in seqsubset:
            for pos in self[k].pos:
                counts[pos] += 1
        
        for c in counts:
            yield c
            
    def iter_columns_with_minimum_char_count(self, x, seqsubset = None):
        for i,c in enumerate(self.iter_column_character_count(seqsubset)):
            if c >= x:
                yield i

    def iter_columns_with_maximum_char_count(self, x, seqsubset = None):
        for i,c in enumerate(self.iter_column_character_count(seqsubset)):
            if c <= x:
                yield i
            
                
    def get_insertion_columns(self,shared):
        return set(i for i in self.iter_columns_with_maximum_char_count(0,shared)) 

    def merge_in(self, she):
        '''
        Merges she inside self, assuming we share some common taxa, and the 
        alignment of common taxa is identical across both alignments.
        
        When assumptions are not met, behavior is largely undefined. 
        '''      
        #global _T_ID
        #_T_ID += 1
        #ID = _T_ID    
        #TIMING_LOG.info("transitivitymerge (%d) started" %ID )    
        mykeys = set(self.keys())
        herkeys = set(she.keys())
        #_LOG.debug("Transitive Merge Started. ID:%d - Rows: %d,%d" %(ID,len(mykeys),len(herkeys)))    
        shared = mykeys.intersection(herkeys)
        #_LOG.debug("Shared seq: %d" %(len(shared)))        
        onlyhers = herkeys - shared
        me_ins = self.get_insertion_columns(shared)
        she_ins = she.get_insertion_columns(shared)
        #_LOG.debug("Insertion Columns: %d,%d" %(len(me_ins),len(she_ins)))
        
        memap=[]
        shemap=[]
        
        melen = self.colcount
        shelen =  she.colcount
    
        ime = 0
        ishe = 0
        inew = 0
        while ime < melen or ishe < shelen:
            if ime in me_ins:              
                memap.append(inew)
                ime += 1
            elif ishe in she_ins:
                shemap.append(inew)
                ishe += 1
            else:       
                memap.append(inew)
                shemap.append(inew)
                ishe += 1
                ime += 1
            inew += 1
            
        self.colcount = inew
        
        for seq in self.values():
            seq.pos = [memap[p] for p in seq.pos]
            
        for k in onlyhers:
            she[k].pos = [shemap[p] for p in she[k].pos]
            self[k] = she[k]        
            
        #TIMING_LOG.info("transitivitymerge (%d) finished" %ID )
        #_LOG.debug("Transitive Merge Finished. ID:%d; cols after: %d" %(ID,self.colcount))

    def mask_gapy_sites(self,minimum_seq_requirement):                
        #_LOG.debug("Masking alignment sites with fewer than %d characters from alignment with %d columns" 
        #           %(minimum_seq_requirement,self.colcount))

        masked = set(self.iter_columns_with_maximum_char_count(minimum_seq_requirement-1))
        
        #_LOG.debug("%d Columns identified for masking" %len(masked))
        
        self.mask_sites(masked)
        
    def mask_unaligned_sites(self):
        #_LOG.debug("Masking alignment sites with lower case letters from an alignment with %d sites" 
        #           %(self.colcount))

        masked = set()
        for seq in self.values():
            for c,i in zip(seq.seq,seq.pos):
                if c > 'a' and c < 'z':
                    masked.add(i)
        
        #_LOG.debug("%d Columns identified for masking" %len(masked))
        
        self.mask_sites(masked)     
        
        return self   
        
    def mask_sites(self, masked):
        if not masked:
            return
        
        off = 0
        colmap = []
        for c in range(0,self.colcount):
            if c in masked:
                off += 1
                colmap.append(-1)
            else:
                colmap.append(c - off)
                
        #_LOG.debug("Column index mapping calculated.")
        for seq in self.values():
            # Find ID of char positions that would be masked
            maskind = [i for i,c in enumerate(seq.pos) if c in masked]
            n = len(seq.pos)
            included = [z for z in reduce(lambda x,y: x+[(x[-1][1]+1,y)],maskind,[(-1,-1)]) if z[0]!=z[1]]
            if not maskind:
                included.append((0,n))
            elif (not included or included[-1][1] < n) and maskind[-1]+1 != n:
                included.append((maskind[-1]+1,n))
            tmp = []
            for (i,j) in included:
                tmp.extend(seq.seq[i:j])
            seq.seq = "".join(tmp)
            seq.pos = [colmap[x] for x in (p for p in seq.pos if p not in masked)]
                
        self.colcount -= off
        #_LOG.debug("Masking done. Before masking: %d; After masking: %d;" 
        #           %(self.colcount+off,self.colcount))
        
    def read_filepath(self, filename, file_format='FASTA'):
        """Augments the matrix by reading the filepath.
        If duplicate sequence names are encountered then the old name will be replaced.
        """
        file_obj = open(filename, 'r')
        ret = self.read_file_object(file_obj, file_format=file_format)
        file_obj.close()
        return ret

    def read_file_object(self, file_obj, file_format='FASTA'):
        """Augments the matrix by reading the file object.
        If duplicate sequence names are encountered then the old name will be replaced.
        """
        if ( file_format.upper() == 'FASTA' ):
            read_func = read_fasta        
        #elif ( file_format.upper() == 'COMPACT' ):
        #    read_func = read_compact
        #elif ( file_format.upper() == 'COMPACT3' ):
        #    read_func = read_compact3
        else:
            raise NotImplementedError("Unknown file format (%s) is not supported" % file_format)
        self.colcount = 0
        for name, seq in read_func(file_obj):
            cseq, l = self.get_alignment_seq_object(seq)
            self[name] = cseq
            self.colcount = max(l, self.colcount)
        
    def as_string_sequence(self,name):
        seq = self[name]
        return seq.as_string(self.colcount)
            
    def get_alignment_seq_object(self, seq):
        cseq = AlignmentSequence()
        l = 0
        if isinstance(seq, str):            
            for m in re.finditer(nogappat, seq):
                cseq.pos.extend(range(m.start(),m.end()))
            cseq.seq = re.sub(gappat,"",seq)
            l = len(seq)
        else:
            cseq.seq = seq[0]
            cseq.pos = seq[1]
            l = seq[1][-1] + 1
        return (cseq,l)

    def update_dict_from(self, alignment):
        for k in self.keys():
            alignment[k] = self.as_string_sequence(k)
        alignment.datatype = self.datatype
            
    def update_from_alignment(self, alignment):
        for k,v in alignment.items():
            self[k],l = self.get_alignment_seq_object(v)
        self.colcount = l
        self.datatype = alignment.datatype
            
    def write_filepath(self, filename, file_format='FASTA', zipout=False):
        """Writes the sequence data in the specified `file_format` to `filename`"""
        
        file_obj = open_with_intermediates(filename,'w')
        if zipout:
            file_obj.close()
            file_obj = StringIO()
        self.write(file_obj, file_format=file_format)
        if zipout:
            import gzip
            file_obj_gz = gzip.open(filename, "wb", 6)
            file_obj_gz.write(str.encode(file_obj.getvalue()))
            file_obj_gz.close()
        file_obj.close()

    def write(self, file_obj, file_format):
        """Writes the sequence data in the specified `file_format` to `file_obj`"""
        if ( file_format.upper() == 'FASTA' ):
            write_func = write_compact_to_fasta        
        #elif ( file_format.upper() == 'COMPACT' ):
        #    write_func = write_compact_to_compact
        #elif ( file_format.upper() == 'COMPACT3' ):
        #    write_func = write_compact_to_compact3 
        #elif ( file_format.upper() == 'PHYLIP' ):
        #    write_func = write_compact_to_phylip 
        else:
            write_func = write_compact_to_fasta
        write_func(self, file_obj)

def compact(alg):
    comp = CompactAlignment()
    comp.update_from_alignment(alg)
    return comp

################# Adapted from SEPP/alignment.py ################
class ReadOnlyAlignment(Mapping, object):
    def get_num_taxa(self):
        return len(self.get_sequence_names())
    
    def get_length(self):
        if self.is_empty():
            raise ValueError("The alignment is empty.\n")
        else:
            return len(next(iter(self.values())))

    def get_sequence_names(self):
        return list(self.keys())

    def 



################## Adapted from SEPP/jobs.py ####################
def readHMMSearch(path):
    assert isinstance(path, str), "Need to give me a path to HMMSearch result"
    assert os.path.isfile(path), "File path does not exist"
    f = open(path, 'r')
    lines = f.read().split('\n')

    results = {}
    pattern = re.compile(
    r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+"
    r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)")
    start_reading = False

    for line in lines:
        line = line.strip()
        # start reading when we see "E-value"
        if not start_reading and line.startswith("E-value") is True:
            start_reading = True
        # end reading if encounter empty line
        elif start_reading and line == "":
            start_reading = False
            break
        # extend the reading
        elif start_reading:
            matches = pattern.search(line)
            if matches is not None and matches.group(0).find("--") == -1:
                results[matches.group(9).strip()] = (
                        float(matches.group(1).strip()),
                        float(matches.group(2).strip()))
    f.close()
    return results
