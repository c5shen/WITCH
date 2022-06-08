'''
Created on 1.22.2022 by Chengze Shen

Global configuration class.
'''

import os
import time
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
from argparse import ArgumentParser, Namespace
from platform import platform

_root_dir = os.path.dirname(os.path.abspath(__file__))
main_config_path = os.path.join(_root_dir, 'main.config')

'''
Configurations defined by users
'''
class Configs:
    global _root_dir

    hmmdir = None
    input_path = None
    backbone_path = None
    backbone_tree_path = None
    query_path = None
    outdir = None
    output_path = None

    keeptemp = False
    keepsubalignment = False
    
    num_hmms = 4
    use_weight = True 
    weight_adjust = 'none'
    subset_size = 1
    num_cpus = -1
    molecule = 'dna'

    # hmmalign/hmmsearch/magus paths
    magus_path = os.path.join(_root_dir, 'tools/magus/magus.py')
    if 'macOS' in platform():
        hmmer_dir = os.path.join(_root_dir, 'tools/hmmer_macOS')
        fasttree_path = os.path.join(_root_dir, 'tools/fasttree_macOS/FastTreeMP')
    else:
        hmmer_dir = os.path.join(_root_dir, 'tools/hmmer')
        fasttree_path = os.path.join(_root_dir, 'tools/fasttree/FastTreeMP')
    hmmalign_path = os.path.join(hmmer_dir, 'hmmalign')
    hmmsearch_path = os.path.join(hmmer_dir, 'hmmsearch')
    hmmbuild_path = os.path.join(hmmer_dir, 'hmmbuild')

    log_path = None
    error_path = None
    debug_path = None
    runtime_path = None

    # MAGUS/GCM configurations
    keepgcmtemp = False
    inflation_factor = 4
    graphclustermethod = 'mcl'
    graphtracemethod = 'minclusters'
    graphtraceoptimize = 'false'
    timeout = 60

    @staticmethod
    def warning(msg, path=None):
        path = Configs.log_path if path is None else path
        Configs.write(msg, 'WARNING', path)

    @staticmethod
    def log(msg, path=None):
        path = Configs.log_path if path is None else path
        Configs.write(msg, 'LOG', path)

    @staticmethod
    def debug(msg, path=None):
        path = Configs.debug_path if path is None else path
        Configs.write(msg, 'DEBUG', path)

    @staticmethod
    def error(msg, path=None):
        path = Configs.error_path if path is None else path
        Configs.write(msg, 'ERROR', path)
    
    @staticmethod
    def runtime(msg, path=None):
        path = Configs.runtime_path if path is None else path
        with open(path, 'a') as f:
            f.write('{}\n'.format(msg))

    @staticmethod
    def write(msg, level, path):
        if path is not None:
            with open(path, 'a') as f:
                f.write('{}\t[{}] {}\n'.format(time.strftime('%Y-%m-%d %H:%M:%S'),
                    level, msg))

# valid configuration check
def valid_configuration(name, conf):
    assert isinstance(conf, Namespace), \
            'Looking for Namespace object but find {}'.format(type(conf))

    # backbone alignment settings
    if name.lower() == 'backbone':
        for k in conf.__dict__.keys():
            attr = getattr(conf, k)
            if not attr:
                continue

            if k == 'alignment_method':
                assert str(attr).lower() in ['magus', 'pasta', 'mafft'], \
                    'Backbone alignment method {} not implemented'.format(attr)
            elif k == 'backbone_size':
                assert int(attr) > 0, 'Backbone size needs to be > 0'
            elif k == 'selection_strategy':
                assert str(attr).lower() in ['median_length', 'random'], \
                    'Selection strategy {} not implemented'.format(attr)
            elif k == 'path':
                assert os.path.exists(os.path.abspath(str(attr))), \
                    '{} does not exist'.format(os.path.abspath(str(attr)))
    

# valid attribute check
def valid_attribute(k, v):
    assert isinstance(k, str)
    if isinstance(v, staticmethod):
        return False
    if not k.startswith('_'):
        return True
    return False

# print a list of all configurations
def getConfigs():
    print('\n********* Configurations **********')
    for k, v in Configs.__dict__.items():
        if valid_attribute(k, v):
            print('\tConfigs.{}: {}'.format(k, v))

'''
Read in from config file if it exists. Any cmd-line provided configs will
override the config file.

Original functionality comes from SEPP -> sepp/config.py
'''
def _read_config_file(filename, opts, expand=None):
    Configs.debug('Reading config from {}'.format(filename))
    config_defaults = []
    cparser = configparser.ConfigParser()
    cparser.optionxform = str
    cparser.read_file(filename)

    if cparser.has_section('commandline'):
        for k, v in cparser.items('commandline'):
            config_defaults.append('--{}'.format(k))
            config_defaults.append(v)

    for section in cparser.sections():
        if section == 'commandline':
            continue
        if getattr(opts, section, None):
            section_name_space = getattr(opts, section)
        else:
            section_name_space = Namespace()
        for k, v in cparser.items(section):
            if expand and k == 'path':
                v = os.path.join(expand, v)
            section_name_space.__setattr__(k, v)
        opts.__setattr__(section, section_name_space)
    return config_defaults

'''
Build configurations
'''
#def buildConfigs(parser, cmdline_args):
def buildConfigs(args):
    ## read in from config file first
    #global main_config_path
    #opts = Namespace()
    #main_cmd_defaults = []

    #if os.path.exists(main_config_path):
    #    #print('Found main configuration file at {}, '.format(
    #    #    main_config_path) + 'loading in...')
    #    with open(main_config_path, 'r') as cfile:
    #        main_cmd_defaults = _read_config_file(cfile, opts)
    #input_args = main_cmd_defaults + cmdline_args

    #args = parser.parse_args(input_args, namespace=opts)

    if args.input_path != None:
        Configs.input_path = os.path.abspath(args.input_path)
    if args.hmmdir != None:
        Configs.hmmdir = os.path.abspath(args.hmmdir)
    if args.backbone_path != None:
        Configs.backbone_path = os.path.abspath(args.backbone_path)
    if args.backbone_tree_path != None:
        Configs.backbone_tree_path = os.path.abspath(args.backbone_tree_path)
    if args.query_path != None:
        Configs.query_path = os.path.abspath(args.query_path)
    
    Configs.outdir = os.path.abspath(args.outdir)
    if not os.path.exists(Configs.outdir):
        os.makedirs(Configs.outdir)
    Configs.output_path = os.path.join(Configs.outdir, args.output_path)

    Configs.keeptemp = args.keeptemp
    Configs.keepsubalignment = args.keepsubalignment

    Configs.log_path = os.path.join(Configs.outdir, 'log.txt')
    Configs.error_path = os.path.join(Configs.outdir, 'error.txt')
    Configs.debug_path = os.path.join(Configs.outdir, 'debug.txt')
    Configs.runtime_path = os.path.join(Configs.outdir, 'runtime_breakdown.txt')

    if args.num_hmms > 0:
        Configs.num_hmms = args.num_hmms

    Configs.use_weight = args.use_weight == 1
    Configs.weight_adjust = args.weight_adjust 

    if args.subset_size > 0:
        if args.subset_size >= 25:
            Configs.warning('Subset size is recommended to be <= 25! Current: {}'.format(
                args.subset_size))
        Configs.subset_size = args.subset_size
    else:
        Configs.warning('Subset size given was < 0! Using default value: 1')

    if args.num_cpus > 0:
        Configs.num_cpus = args.num_cpus
    else:
        Configs.num_cpus = os.cpu_count()

    Configs.molecule = args.molecule

    # MAGUS/GCM options
    Configs.keepgcmtemp = args.keepgcmtemp
    Configs.inflation_factor = args.inflation_factor
    Configs.graphclustermethod = args.graphclustermethod
    Configs.graphtracemethod = args.graphtracemethod
    Configs.graphtraceoptimize = args.graphtraceoptimize
    
    # additional MAGUS/GCM failsafe option to timeout a MAGUS process
    # if [timeout] seconds are reached before the process finished
    Configs.timeout = args.timeout

    # add any additional arguments to Configs
    for k in args.__dict__.keys():
        if k not in Configs.__dict__:
            k_attr = getattr(args, k)

            # check whether the configuration is valid
            valid_configuration(k, k_attr)

            setattr(Configs, k, k_attr)

    ## backbone options
    #if args.backbone_size != None:
    #    setattr(Configs.backbone, 'backbone_size', args.backbone_size)
    #if args.selection_strategy != None:
    #    setattr(Configs.backbone, 'selection_strategy', args.selection_strategy)
    #if args.backbone_method != None:
    #    setattr(Configs.backbone, 'alignment_method', args.backbone_method)

    ## Bandaid way of logging main.config to log.txt
    #Configs.log('WITCH is running with: {}'.format(' '.join(cmdline_args)))
    #if os.path.exists(main_config_path):
    #    Configs.log('Main configuration loaded from {}'.format(
    #        main_config_path))
