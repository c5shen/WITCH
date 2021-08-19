import os
import time
import numpy as np
from configs import Configs
from multiprocessing import Pool

class Weights(object):
    weights = dict()
    weights_map = dict()
    ranked_bitscores = dict()
    def __init__(self):
        pass

'''
Function to calculate the HMM weighting, given the bitscores and sizes
of the HMMs (for a given query taxon)
inputs: ensemble of HMMs H (with their bitscores and sizes)
outputs: weights for HMMs H
'''
def calculateWeights(taxon, indexes, bitscores, sizes):
    #logging.debug('working with: {}'.format(taxon))
    weights = {}
    
    #temp = []
    assert len(indexes) == len(bitscores) == len(sizes)
    for i in range(len(bitscores)):
        score_i, size_i = bitscores[i], sizes[i]
        exponents = np.array(bitscores) - score_i \
                + np.log2(np.array(sizes) / size_i)
        denominator = np.sum(np.power(2, exponents))
        #temp.append(1. / denominator)
        weights[indexes[i]] = 1. / denominator
    return taxon, weights

'''
Obtain and write weights to local based on bitscores
'''
def loadWeights(index_to_hmm, ranks):
    s2 = time.time()
    pool = Pool(Configs.num_cpus)

    # - get sizes of each HMM
    all_sizes = {}
    for index, subset in index_to_hmm.items():
        all_sizes[index] = subset.alignment.get_num_taxa()

    # iterate through each query taxon
    weights, weights_map = {}, {}
    args = []
    for taxon, sorted_scores in ranks.items():
        indexes = [x[0] for x in sorted_scores]
        bitscores = [x[1] for x in sorted_scores]
        sizes = [all_sizes[x] for x in indexes]
        args.append([taxon, indexes, bitscores, sizes])

        #this_weights_map,_ = calculateWeights(taxon, indexes, bitscores, sizes)
        #weights[taxon] = sorted([(ind, w) for ind, w in this_weights_map.items()],
        #        key = lambda x: x[1], reverse=True)
        #weights_map[taxon] = this_weights_map
    all_weights_and_temps = pool.starmap(calculateWeights, args)
    pool.close()
    pool.join()

    for item in all_weights_and_temps:
        weights[item[0]] = sorted([(ind, w) for ind, w in item[1].items()],
                key = lambda x: x[1], reverse=True)
        weights_map[item[0]] = item[1]
    Weights.weights = weights
    Weights.weights_map = weights_map

    # write weights to local
    with open('{}/all_weights.txt'.format(Configs.outdir), 'w') as f:
        for taxon, weight in weights.items():
            w = [str(x) for x in weight]
            f.write(taxon + ':' + ';'.join(w) + '\n')
    time_obtain_weights = time.time() - s2
    Configs.warning('Finished calculating weights!')
    Configs.runtime('Time to obtain weights given bitscores (s): {}'.format(
        time_obtain_weights))
