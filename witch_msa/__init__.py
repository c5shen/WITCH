#!/usr/bin/env python3
import os, sys, time
from argparse import ArgumentParser, Namespace

from witch_msa.configs import _root_dir, main_config_path, _read_config_file 
from witch_msa.configs import *
from witch_msa.gcmm.gcmm import mainAlignmentProcess
from witch_msa.helpers.general_tools import SmartHelpFormatter

__version__ = "1.0.6"

def witch_runner():
    global _root_dir, main_config_path 
    parser = _init_parser()
    cmdline_args = sys.argv[1:]

    buildConfigs(parser, cmdline_args)
    getConfigs()

    Configs.log('WITCH is running with: {}'.format(' '.join(cmdline_args)))
    if os.path.exists(main_config_path):
        Configs.log('Main configuration loaded from {}'.format(
            main_config_path))

    # run the codes
    s1 = time.time()
    mainAlignmentProcess(parser, cmdline_args)
    s2 = time.time()

    Configs.log('WITCH finished in {} seconds...'.format(s2 - s1))
    print('\nAll done! WITCH finished in {} seconds...'.format(s2 - s1))

def _init_parser():
    parser = ArgumentParser(
            description=(
                "This program runs WITCH, an alignment method "
                "extended from UPP and MAGUS."),
            conflict_handler='resolve',
            formatter_class=SmartHelpFormatter)
    parser.add_argument('-v', '--version', action='version',
            version="%(prog)s " + __version__)

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
    basic_group.add_argument('-t', '--num-cpus', type=int,
            help='Number of cpus for multi-processing. Default: -1 (all)',
            required=False, default=-1)
    basic_group.add_argument('--max-concurrent-jobs', type=int, metavar='INT',
            help=' '.join(['Maximum number of concurrently running jobs.'
                    'Default: min(50, 10 * num_cpus)']),
            required=False, default=None)
    basic_group.add_argument('--timeout', type=int, metavar='INT',
            help=' '.join(['Retry a query alignment after [timeout] seconds.',
                    'The retry will always use \'witch-ng mode\' (see WTICH',
                    'option below). Default: 120']),
            default=120, required=False)
    # Added @ 7.25.2024 - Chengze Shen
    # new arg parameter -c to read in user-customized config files
    # usage --> override any existing parameter settings (except cmdline)
    # priority: cmdline > user.config > main.config
    basic_group.add_argument('-c', '--config-file', type=str, required=False,
            help=' '.join(['User-specified config file. Override main.config',
                'but commandline parameters still have the highest priority.']),
            default=None)

    #basic_group.add_argument('--chunksize', type=int,
    #        help='Chunksize for multiprocessing', required=False,
    #        default=1)
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
    witch_group.add_argument('-m', '--mode', type=str,
            help=' '.join(['Mode for aligning query sequences.',
                    '\nIn the case of adding query sequences into a large',
                    'backbone alignment, the original script used by WITCH',
                    'is inefficient and can lead to extremely long runtime.',
                    '\nThe re-implementation of WITCH-ng\'s way of adding',
                    'sequences (i.e., with DP algorithm) now enables WITCH',
                    'to be much faster and memory-efficient (with some',
                    'additional tweaks). It is advised to use',
                    '\'witch-ng mode\' (enabled as default) within WITCH.',
                    '\nDefault: witch-ng']),
            required=False, default='witch-ng',
            choices=['old-witch', 'witch-ng'])
    witch_group.add_argument('--keeptemp', action='store_const', const=True,
            help='Keep ALL temporary files in the process (constraints' \
                    + ', backbones, HMMSearch results, GCM results, etc.)',
            default=False)
    witch_group.add_argument('--keep-decomposition', type=int, metavar='{0,1}',
            help=' '.join(['Whether to keep the tree decomposition',
                    '(including temp backbone/frag files and',
                    'HMMBuild/HMMSearch results). Default: 1']),
            default=1, required=False)
    #witch_group.add_argument('--keepsubalignment',
    #        action='store_const', const=True,
    #        help='Keep all subalignments by MAGUS/GCM', default=False)
    witch_group.add_argument('-k', '--num-hmms', type=int,
            help=' '.join(['The number of top-scored HMMs used for aligning',
                    'a query. Default: 10']),
            required=False, default=10)
    witch_group.add_argument('-w', '--use-weight',
            metavar='{0,1}', type=int, required=False,
            help='Whether to use adjusted bitscore (weights). Default: 1',
            default=1)
    witch_group.add_argument('--save-weight',
            metavar='{0,1}', type=int, required=False,
            help='Whether to save weights to [outdir]/weights.txt. Default: 0',
            default=0)
    witch_group.add_argument('-A', '--alignment-size', type=int,
            help=' '.join(['Minimum decomposition subset size for tree',
                    'decomposition. Default: 10 (as UPP)']),
            default=10, required=False)
    # Added @ 1.14.2024 - Chengze Shen
    # option to set an upper bound to decomposition subset size (complementary
    # to "-A" but as an upper bound).
    witch_group.add_argument('-Z', '--alignment-upper-bound', type=int,
            help=' '.join(['Maximum decomposition subset size for tree',
                    'decomposition. Default: None (as UPP)']),
            default=None, required=False)
    witch_group.add_argument('--molecule', type=str,
            help='Whether input is amino/dna/rna. Default: infer from data',
            required=False, default=None, choices=['amino', 'dna', 'rna'])
    #witch_group.add_argument('-s', '--subset-size', type=int,
    #        help='Number of queries in a single GCM run, default: 1',
    #        required=False, default=1)
    #witch_group.add_argument('--weight-adjust', type=str, required=False,
    #        default='none', choices=['none', 'normalize', 'maxto1'],
    #        help='(DEPRECATED) Optional adjustment of weights, default: none')

    # miscellaneous group
    misc_group = parser.add_argument_group(
            "Miscellaneous options".upper(),
            ' '.join(["Optional parameters for WITCH setup/config etc."]))
    parser.groups['misc_group'] = misc_group
    misc_group.add_argument('-y', '--bypass-setup', action='store_const',
            const=True, default=False,
            help=' '.join(['(Optional) include to bypass the initial',
                    'step when running WITCH to set up the configuration',
                    'directory (will use ~/.witch_msa).',
                    'Note: you only need to use this option once.']))

    # GCM option
    #gcm_group = parser.add_argument_group(
    #        'MAGUS/GCM options'.upper(),
    #        ' '.join(['These options are used to customize MAGUS/GCM,',
    #            'for example the graph trace method.']))
    #parser.groups['gcm_group'] = gcm_group
    #gcm_group.add_argument('--keepgcmtemp', action='store_const',
    #        const=True, default=False, required=False,
    #        help='Keep temporary files generated from MAGUS/GCM')
    #gcm_group.add_argument('-f', '--inflation-factor', type=float,
    #        default=4, required=False,
    #        help="Inflation factor for MCL.\n[default: 4]") 
    #gcm_group.add_argument('--graphclustermethod', type=str,
    #        choices=['mcl', 'mlrmcl', 'rg', 'none'],
    #        default='mcl', required=False,
    #        help="Method for initial clustering of the alignment graph.\n[default: mcl]")
    #gcm_group.add_argument('--graphtracemethod', type=str,
    #        choices=['minclusters', 'mwtgreedy', 'mwtsearch', 'fm', 'rg', 'rgfast'],
    #        default='minclusters', required=False,
    #        help="Method for finding a trace from the alignment graph.\n[default: minclusters]")
    #gcm_group.add_argument('--graphtraceoptimize', type=str,
    #        choices=['true', 'false'], required=False, default='false',
    #        help="Run an optimization step on the graph trace.\n[default: false]")
    return parser
