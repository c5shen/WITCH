import os, sys, time
import argparse
import logging

from configs import *
from gcmm.gcmm import mainAlignmentProcess

def main():
    args = parseArgs()
    buildConfigs(args)
    Configs.log('GCM+eHMMs is running with: {}'.format(' '.join(sys.argv)))
    getConfigs()

    # run the codes
    s1 = time.time()
    mainAlignmentProcess()
    s2 = time.time()

    Configs.log('GCM+eHMMs finished in {} seconds...'.format(s2 - s1))

def parseArgs():
    parser = argparse.ArgumentParser()

    # requirements
    required_group = parser.add_argument_group(
            "Required Parameters".upper(),
            ' '.join(["These are required fields, assuming that an UPP",
                "output and its temporary files are available."]))
    parser.groups = dict()
    parser.groups['required_group'] = required_group
    required_group.add_argument('-p', '--hmmdir', type=str,
            help='Path to the HMMs directory generated from UPP', required=True)
    required_group.add_argument('-b', '--backbone-path', type=str,
            help='Path to the backbone alignment from/for UPP', required=True)
    required_group.add_argument('-q', '--query-path', type=str,
            help='Path to the queries file that we want to align', required=True)
    required_group.add_argument('-o', '--outdir', type=str,
            help='Output directory, default: gcm+eHMMs_output/', required=False,
            default='gcm+eHMMs_output')

    # GCM+eHMMs options
    gcm_eHMMs_group = parser.add_argument_group(
            "GCM+eHMMs options".upper(),
            ' '.join(["These options are used to customize the GCM+eHMMs",
                "pipeline."]))
    parser.groups['gcm_eHMMs_group'] = gcm_eHMMs_group
    gcm_eHMMs_group.add_argument('--keeptemp', action='store_const', const=True,
            help='Keep temporary files in the process (constraints' \
                    + ', HMMSearch results, etc.)', default=False)
    gcm_eHMMs_group.add_argument('--keepsubalignment',
            action='store_const', const=True,
            help='Keep all subalignments by MAGUS/GCM', default=False)
    gcm_eHMMs_group.add_argument('-k', '--num-hmms', type=int,
            help='The number of top HMMs used for aligning a query',
            required=False, default=4)
    gcm_eHMMs_group.add_argument('-w', '--use-weight',
            action='store_const', const=True, required=False,
            help='Whether to use adjusted bitscore (weights), default: 0',
            default=False)
    gcm_eHMMs_group.add_argument('-s', '--subset-size', type=int,
            help='Number of queries in a single GCM run, default: 1',
            required=False, default=1)
    gcm_eHMMs_group.add_argument('--weight-adjust', type=str, required=False,
            default='none', choices=['none', 'normalize', 'maxto1'],
            help='Optional adjustment of weights, default: none')
    gcm_eHMMs_group.add_argument('-t', '--num-cpus', type=int,
            help='Number of cpus for multi-processing, default: -1 (all)',
            required=False, default=-1)

    # GCM option
    gcm_group = parser.add_argument_group(
            'MAGUS/GCM options'.upper(),
            ' '.join(['These options are used to customize MAGUS/GCM,',
                'for example the graph trace method.']))
    parser.groups['gcm_group'] = gcm_group
    gcm_group.add_argument('--keepgcmtemp', action='store_const',
            const=True, default=False, required=False,
            help='Keep temporary files generated from MAGUS/GCM')
    gcm_group.add_argument('-f', '--inflation-factor', type=float,
            default=4, required=False,
            help="Inflation factor for MCL, default: 4") 
    gcm_group.add_argument('--graphclustermethod', type=str,
            choices=['mcl', 'mlrmcl', 'rg', 'none'],
            default='mcl', required=False,
            help="Method for initial clustering of the alignment graph, default: mcl")
    gcm_group.add_argument('--graphtracemethod', type=str,
            choices=['minclusters', 'mwtgreedy', 'mwtsearch', 'fm', 'rg', 'rgfast'],
            default='minclusters', required=False,
            help="Method for finding a trace from the alignment graph, default: minclusters")
    gcm_group.add_argument('--graphtraceoptimize', type=str,
            choices=['true', 'false'], required=False, default='false',
            help="Run an optimization step on the graph trace, default: false")
    return parser.parse_args()


if __name__ == '__main__':
    main()
