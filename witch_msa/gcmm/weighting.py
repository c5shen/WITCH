'''
Created on 10.28.2021 by Chengze Shen

Bitscore to weight calculation.
'''

import os
import time
import numpy as np
import concurrent.futures
from witch_msa.configs import Configs, tqdm_styles
from tqdm import tqdm

class Weights(object):
    weights = dict()
    weights_map = dict()
    ranked_bitscores = dict()
    def __init__(self):
        pass

'''
Function to read in weights from local given the taxon name
'''
def readWeights(taxon):
    infile = Configs.outdir + '/weights/w_{}.txt'.format(taxon)
    if not os.path.isfile(infile):
        return None, None
    else:
        weights, weights_map = [], [] 
        with open(infile, 'r') as f:
            line = f.read()
            taxon, raw = line.split(':')
            weights = [eval(x) for x in raw.split(';')]
            weights_map = {ind: w for (ind, w) in weights}
        return weights, weights_map

'''
Function to read in bitscores from local given the taxon name
'''
def readBitscores(taxon):
    infile = Configs.outdir + '/bitscores/b_{}.txt'.format(taxon)
    if not os.path.isfile(infile):
        return None, None
    else:
        bitscores = [] 
        with open(infile, 'r') as f:
            line = f.read()
            taxon, raw = line.split(':')
            bitscores = [eval(x) for x in raw.split(';')]
        return bitscores

'''
Function to calculate the HMM weighting, given the bitscores and sizes
of the HMMs (for a given query taxon)
inputs: ensemble of HMMs H (with their bitscores and sizes)
outputs: weights for HMMs H
'''
def calculateWeights(packed_data):
    taxon, indexes, bitscores, sizes = packed_data
    #logging.debug('working with: {}'.format(taxon))
    weights = {}
    
    assert len(indexes) == len(bitscores) == len(sizes)
    for i in range(len(bitscores)):
        score_i, size_i = bitscores[i], sizes[i]
        exponents = np.array(bitscores) - score_i \
                + np.log2(np.array(sizes) / size_i)
        denominator = np.sum(np.power(2, exponents))
        weights[indexes[i]] = 1. / denominator
    
    num_to_retain = min(Configs.num_hmms, len(weights))
    sorted_weights = sorted([(ind, w) for ind, w in weights.items()],
            key = lambda x: x[1], reverse=True)[:num_to_retain]
    return {taxon: tuple(sorted_weights)}

    ## write weights to local (only top k ones)
    #sorted_weights = [str(x) for x in sorted_weights]
    #with open(Configs.outdir + '/weights/w_{}.txt'.format(taxon), 'w') as f:
    #    f.write(taxon + ':' + ';'.join(sorted_weights) + '\n')
    #return None

'''
Function to write a single taxon with its ranked bitscore to local
'''
def writeQueryBitscores(packed_data):
    taxon, sorted_scores = packed_data
    str_sorted_scores = [str(x) for x in sorted_scores]

    with open(Configs.outdir + '/bitscores/b_{}.txt'.format(taxon), 'w') as f:
        f.write(taxon + ':' + ';'.join(str_sorted_scores) + '\n')
    return None

'''
Write bitscores to local (the same way as we write weights)
'''
def writeBitscores(ranked_bitscores, pool):
    s2 = time.time()
    Configs.warning('Starting to load bitscores...')
    #if not os.path.isdir(Configs.outdir + '/bitscores'):
    #    os.makedirs(Configs.outdir + '/bitscores')

    taxon_to_bitscores = {}
    for taxon, sorted_scores in ranked_bitscores.items():
        num_to_retain = min(Configs.num_hmms, len(sorted_scores))
        taxon_to_bitscores[taxon] = tuple(sorted_scores[:num_to_retain])

    #args = []
    #for taxon, sorted_scores in ranked_bitscores.items():
    #    args.append((taxon, sorted_scores))
    #all_score_temps = list(pool.map(writeQueryBitscores, args))

    time_write_scores = time.time() - s2
    Configs.warning('Finished loading bitscores in memory.')
    Configs.runtime(' '.join(['(writeBitscores) Time to write ranked bitscores',
            'to local (s):', str(time_write_scores)]))
    return taxon_to_bitscores

'''
Obtain and write weights to local based on bitscores
'''
def writeWeights(index_to_hmm, ranked_bitscores, pool):
    s2 = time.time()
    Configs.warning('Starting to calculate weights...')
    #pool = Pool(Configs.num_cpus)

    # - get sizes of each HMM
    all_sizes = {}
    for index, subset in index_to_hmm.items():
        all_sizes[index] = subset.num_taxa
        #all_sizes[index] = subset.alignment.get_num_taxa()

    # iterate through each query taxon
    # write to local for each taxon and its weights
    #if not os.path.isdir(Configs.outdir + '/weights'):
    #    os.makedirs(Configs.outdir + '/weights')
    weights, weights_map = {}, {}
    args = []
    for taxon, sorted_scores in ranked_bitscores.items():
        indexes = [x[0] for x in sorted_scores]
        bitscores = [x[1] for x in sorted_scores]
        sizes = [all_sizes[x] for x in indexes]
        args.append((taxon, indexes, bitscores, sizes))

        ## sequential version to calculate weights
        #this_weights_map,_ = calculateWeights(taxon, indexes, bitscores, sizes)
        #weights[taxon] = sorted([(ind, w) for ind, w in this_weights_map.items()],
        #        key = lambda x: x[1], reverse=True)
        #weights_map[taxon] = this_weights_map
    #all_taxon_to_weights = list(pool.map(calculateWeights, args,
    #    chunksize=Configs.chunksize))
    all_taxon_to_weights, futures = [], []
    for arg in args:
        futures.append(pool.submit(calculateWeights, arg))
    for future in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(args), **tqdm_styles):
        res = future.result()
        if res:
            all_taxon_to_weights.append(res)

    taxon_to_weights = {}
    for item in all_taxon_to_weights:
        taxon_to_weights.update(item)

    time_obtain_weights = time.time() - s2
    Configs.warning('Finished calculating weights!')
    Configs.runtime(' '.join(['(writeWeights) Time to obtain weights',
            'given bitscores (s):', str(time_obtain_weights)]))
    return taxon_to_weights

'''
Write weights to local as [outdir]/weights.txt
'''
def writeWeightsToLocal(taxon_to_weights, path):
    Configs.log('Writing weights to {}'.format(path))
    with open(path, 'w') as f:
        for taxon, weights in taxon_to_weights.items():
            f.write('{}:{}\n'.format(taxon, weights))

'''
Function to read weights from a given weights path (e.g., ./weights.txt)
Return a dictionary of taxon to weights
'''
def readWeightsFromLocal(path):
    Configs.log('Reading weights from {}'.format(path))
    taxon_to_weights = {}
    with open(path, 'r') as f:
        line = f.readline()
        while line:
            # split by ':'
            taxon, taxon_weight = line.split(':')
            taxon_to_weights[taxon] = eval(taxon_weight)
            line = f.readline()
    return taxon_to_weights
