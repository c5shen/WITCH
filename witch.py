import os, sys, time
from argparse import ArgumentParser, Namespace
import logging

from configs import _read_config_file 
from configs import *
from gcmm.gcmm import mainAlignmentProcess

version = "0.1.1"

def main():
    parser = _init_parser()
    cmdline_args = sys.argv[1:]
    
    global main_config_path 
    opts = Namespace()
    main_cmd_defaults = []

    if os.path.exists(main_config_path):
        with open(main_config_path, 'r') as cfile:
            main_cmd_defaults = _read_config_file(cfile, opts)
    input_args = main_cmd_defaults + cmdline_args
    args = parser.parse_args(input_args, namespace=opts)

    buildConfigs(args)
    getConfigs()

    Configs.log('WITCH is running with: {}'.format(' '.join(cmdline_args)))
    if os.path.exists(main_config_path):
        Configs.log('Main configuration loaded from {}'.format(
            main_config_path))

    # run the codes
    s1 = time.time()
    mainAlignmentProcess(args)
    s2 = time.time()

    Configs.log('WITCH finished in {} seconds...'.format(s2 - s1))

def _init_parser():
    parser = ArgumentParser(
            description=(
                "This program runs WITCH, an alignment method "
                "extended from UPP and MAGUS."),
            conflict_handler='resolve')
    parser.add_argument('-v', '--version', action='version',
            version="%(prog)s " + version)

    # basic settings
    basic_group = parser.add_argument_group(
            "Basic parameters".upper(),
            ' '.join(["These are basic fields.",
                "Users can choose to provide the backbone alignment/tree,",
                "or even the path to the eHMM from a previous UPP run.",
                "Otherwise, WITCH will generate them."]))
    parser.groups = dict()
    parser.groups['basic_group'] = basic_group
    basic_group.add_argument('-i', '--input-path', type=str,
            help='Path to the input unaligned file (all sequences).', required=False)
    basic_group.add_argument('-p', '--hmmdir', type=str,
            help='Path to the HMMs directory generated from UPP', required=False)
    basic_group.add_argument('-b', '--backbone-path', type=str,
            help='Path to the backbone alignment', required=False)
    basic_group.add_argument('-e', '--backbone-tree-path', type=str,
            help='Path to the backbone tree', required=False)
    basic_group.add_argument('-q', '--query-path', type=str,
            help='Path to the queries file that we want to align', required=False)
    basic_group.add_argument('-d', '--outdir', type=str,
            help='Output directory, default: witch_output/', required=False,
            default='witch_output')
    basic_group.add_argument('-o', '--output-path', type=str,
            help='Output file name, default: merged.fasta', required=False,
            default='merged.fasta')
    #basic_group.add_argument('--collapse-singletons', type=int,
    #        default=1, required=False,
    #        help='Collapse singleton columns and make them lower cases. default: 1')

    ## backbone alignment/tree options
    #backbone_group = parser.add_argument_group(
    #        "backbone options".upper(),
    #        ' '.join(["These options control how backbone sequences",
    #            "are selected and aligned."]))
    #parser.groups['backbone_group'] = backbone_group
    #backbone_group.add_argument('--backbone-size', type=int,
    #        help='Number of sequences to include in the backbone, ' \
    #                + 'default: min(1000, len(all_taxa))', default=None)
    #backbone_group.add_argument('--selection-strategy', type=str,
    #        choices=['median_length', 'random'],
    #        help='Backbone sequence selection strategy, default: ' \
    #                + 'sequences with lengths 25% around the median',
    #                default=None)
    #backbone_group.add_argument('--backbone-method', type=str,
    #        choices=['magus', 'pasta', 'mafft'],
    #        help='Alignment method on the backbone sequences, default: ' \
    #                + 'magus', default=None)

    # WITCH options
    witch_group = parser.add_argument_group(
            "WITCH options".upper(),
            ' '.join(["These options are used to customize WITCH",
                "pipeline."]))
    parser.groups['witch_group'] = witch_group
    witch_group.add_argument('--keeptemp', action='store_const', const=True,
            help='Keep ALL temporary files in the process (constraints' \
                    + ', backbones, HMMSearch results, GCM results, etc.)',
            default=False)
    #witch_group.add_argument('--keepsubalignment',
    #        action='store_const', const=True,
    #        help='Keep all subalignments by MAGUS/GCM', default=False)
    witch_group.add_argument('-k', '--num-hmms', type=int,
            help='The number of top-scored HMMs used for aligning a query, default: 10',
            required=False, default=10)
    witch_group.add_argument('-w', '--use-weight',
            type=int, required=False,
            help='Whether to use adjusted bitscore (weights), default: 1',
            default=1)
    #witch_group.add_argument('-s', '--subset-size', type=int,
    #        help='Number of queries in a single GCM run, default: 1',
    #        required=False, default=1)
    witch_group.add_argument('--weight-adjust', type=str, required=False,
            default='none', choices=['none', 'normalize', 'maxto1'],
            help='Optional adjustment of weights, default: none')
    witch_group.add_argument('-t', '--num-cpus', type=int,
            help='Number of cpus for multi-processing, default: -1 (all)',
            required=False, default=-1)
    witch_group.add_argument('--molecule', type=str,
            help='Whether input is amino/dna/rna, default: dna',
            required=False, default='dna', choices=['amino', 'dna', 'rna'])

    # GCM option
    gcm_group = parser.add_argument_group(
            'MAGUS/GCM options'.upper(),
            ' '.join(['These options are used to customize MAGUS/GCM,',
                'for example the graph trace method.']))
    parser.groups['gcm_group'] = gcm_group
    gcm_group.add_argument('--keepgcmtemp', action='store_const',
            const=True, default=False, required=False,
            help='Keep temporary files generated from MAGUS/GCM')
    gcm_group.add_argument('--timeout', type=int,
            default=60, required=False,
            help='Retry a MAGUS/GCM subtask after [timeout] seconds, default: 60')
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
    return parser


if __name__ == '__main__':
    main()
