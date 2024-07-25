'''
Created on 1.22.2022 by Chengze Shen

Global configuration class.
'''

import os
import time
try:
    import configparser
except ImportError:
    from ConfigParser import configparser
from argparse import ArgumentParser, Namespace
from platform import platform
from witch_msa.helpers.alignment_tools import inferDataType
from witch_msa.init_configs import init_config_file

# detect home.path location
homepath = os.path.dirname(__file__) + '/home.path'
_root_dir, main_config_path = init_config_file(homepath)

# default settings for tqdm progress bar style
tqdm_styles = {
        'desc': '\tRunning...', 'ascii': False,
        'ncols': 80, 
        #'disable': True,
        #'colour': 'green',
        'mininterval': 0.5
        }

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
    config_file = None     # Added @ 7.25.2024 

    chunksize = 1

    keeptemp = False
    keep_decomposition = True
    #keepsubalignment = False
    
    # WITCH configurations
    mode = 'witch-ng'
    num_hmms = 10
    use_weight = True 
    save_weight = False
    alignment_size = 10
    alignment_upper_bound = None
    #weight_adjust = 'none'
    #subset_size = 1

    num_cpus = -1
    max_concurrent_jobs = None
    molecule = 'dna'
    #collapse_singletons = True

    # hmmalign/hmmsearch/magus paths
    magus_path = None
    #gcm_path = os.path.join(_root_dir, 'tools/gcm137/gcm137')
    mafftpath = None
    fasttreepath = None
    mclpath = None
    hmmsearchpath = None
    hmmalignpath = None
    hmmbuildpath = None

    # miscellaneous option placeholder (not needed)
    bypass_setup = False

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
    timeout = 120

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

    @staticmethod
    def inferDataType(path):
        if Configs.molecule is None: 
            Configs.molecule = inferDataType(path)
            Configs.log('Molecule type was not specified. Inferred type: {}'.format(
                Configs.molecule))
        return Configs.molecule

# check for valid configurations and set them
def set_valid_configuration(name, conf):
    assert isinstance(conf, Namespace), \
            'Looking for Namespace object but find {}'.format(type(conf))

    # backbone alignment settings
    if name == 'Backbone':
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
                assert os.path.exists(os.path.realpath(str(attr))), \
                    '{} does not exist'.format(os.path.realpath(str(attr)))
        setattr(Configs, name, conf)
    # settings that change basic Configs class variables such as:
    # fasttreepath, hmmalignpath, etc.
    elif name == 'Basic':
        for k in conf.__dict__.keys():
            attr = getattr(conf, k)
            if not attr:
                continue
            # set variable [k] to [attr] if provided
            setattr(Configs, k, attr)
    elif name == 'MAGUS':
        setattr(Configs, name, conf)

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
    print('\thome.path: {}'.format(homepath))
    print('\tmain.config: {}\n'.format(main_config_path))
    for k, v in Configs.__dict__.items():
        if valid_attribute(k, v):
            print('\tConfigs.{}: {}'.format(k, v))

'''
Read in from config file if it exists. Any cmd-line provided configs will
override the config file.

Original functionality comes from SEPP -> sepp/config.py
'''
def _read_config_file(filename, cparser, opts, expand=None):
    Configs.debug('Reading config from {}'.format(filename))
    config_defaults = []
    #cparser = configparser.ConfigParser()
    #cparser.optionxform = str
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
#def buildConfigs(args):
def buildConfigs(parser, cmdline_args):
    # config parser, which first reads in main.config and later overrides
    # with user.config (if specified)
    cparser = configparser.ConfigParser()
    cparser.optionxform = str

    # load default_args from main.config
    default_args = Namespace()
    cmdline_default = []
    with open(main_config_path, 'r') as cfile:
        cmdline_default = _read_config_file(cfile, cparser, default_args)
    
    # load cmdline args first, then search for user.config if specified
    args = parser.parse_args(cmdline_args)
    cmdline_user = []
    if args.config_file != None:
        # override default_args
        Configs.config_file = args.config_file
        with open(Configs.config_file, 'r') as cfile:
            cmdline_user = _read_config_file(cfile, cparser, default_args)

    # finally, re-parse cmdline args in the order:
    #   [cmdline_default, cmd_user, cmdline_args] 
    args = parser.parse_args(cmdline_default + cmdline_user + cmdline_args,
            namespace=default_args)

    if args.input_path != None:
        Configs.input_path = os.path.realpath(args.input_path)
    if args.hmmdir != None:
        Configs.hmmdir = os.path.realpath(args.hmmdir)
    if args.backbone_path != None:
        Configs.backbone_path = os.path.realpath(args.backbone_path)
    if args.backbone_tree_path != None:
        Configs.backbone_tree_path = os.path.realpath(args.backbone_tree_path)
    if args.query_path != None:
        Configs.query_path = os.path.realpath(args.query_path)
    
    Configs.outdir = os.path.realpath(args.outdir)
    if not os.path.exists(Configs.outdir):
        os.makedirs(Configs.outdir)
    Configs.output_path = os.path.join(Configs.outdir, args.output_path)

    Configs.keeptemp = args.keeptemp
    Configs.keep_decomposition = args.keep_decomposition == 1
    #Configs.keepsubalignment = args.keepsubalignment

    Configs.log_path = os.path.join(Configs.outdir, 'log.txt')
    Configs.error_path = os.path.join(Configs.outdir, 'error.txt')
    Configs.debug_path = os.path.join(Configs.outdir, 'debug.txt')
    Configs.runtime_path = os.path.join(Configs.outdir, 'runtime_breakdown.txt')

    # query alignment mode. Default WITCH uses MAGUS/GCM, while I also
    # implemented the WITCH-ng way.
    # TODO: Implement a smart way that mixes the two for runtime, since
    # running Smith-Waterman can be slow on super long sequences
    Configs.mode = args.mode

    if args.num_hmms > 0:
        Configs.num_hmms = args.num_hmms

    Configs.use_weight = args.use_weight == 1
    Configs.save_weight = args.save_weight == 1
    Configs.alignment_size = args.alignment_size
    if args.alignment_upper_bound is not None:
        Configs.alignment_upper_bound = args.alignment_upper_bound \
                if args.alignment_upper_bound > 0 else None
    #Configs.weight_adjust = args.weight_adjust 

    #Configs.chunksize = args.chunksize

    if args.num_cpus > 0:
        Configs.num_cpus = args.num_cpus
    else:
        Configs.num_cpus = os.cpu_count()

    if args.max_concurrent_jobs:
        Configs.max_concurrent_jobs = args.max_concurrent_jobs
    else:
        Configs.max_concurrent_jobs = min(50, 10 * Configs.num_cpus)

    Configs.molecule = args.molecule

    #Configs.collapse_singletons = args.collapse_singletons == 1

    # MAGUS/GCM options
    #Configs.keepgcmtemp = args.keepgcmtemp
    #Configs.inflation_factor = args.inflation_factor
    #Configs.graphclustermethod = args.graphclustermethod
    #Configs.graphtracemethod = args.graphtracemethod
    #Configs.graphtraceoptimize = args.graphtraceoptimize
    
    # additional MAGUS/GCM failsafe option to timeout a MAGUS process
    # if [timeout] seconds are reached before the process finished
    Configs.timeout = args.timeout

    # add any additional arguments to Configs
    for k in args.__dict__.keys():
        if k not in Configs.__dict__:
            k_attr = getattr(args, k)

            # check whether the configuration is valid
            set_valid_configuration(k, k_attr)
            #setattr(Configs, k, k_attr)
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
