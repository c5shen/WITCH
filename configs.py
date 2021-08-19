import os
import time

'''
Configurations defined by users
'''
class Configs:
    hmmdir = None
    backbone_path = None
    query_path = None
    outdir = None
    keeptemp = False
    keepsubalignment = False
    
    num_hmms = 4
    use_weight = False
    weight_adjust = 'none'
    subset_size = 1
    num_cpus = -1

    # hmmalign/hmmsearch/magus paths
    hmmalign_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools/hmmer/hmmalign')
    hmmsearch_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools/hmmer/hmmsearch')
    magus_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools/magus/magus.py')

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

# print a list of all configurations
def getConfigs():
    print('Configs.hmmdir:', Configs.hmmdir)
    print('Configs.backbone_path:', Configs.backbone_path)
    print('Configs.outdir:', Configs.outdir)
    print('Configs.keeptemp:', Configs.keeptemp)
    print('Configs.keepsubalignment:', Configs.keepsubalignment)

    print('\nConfigs.num_hmms:', Configs.num_hmms)
    print('Configs.use_weight:', Configs.use_weight)
    print('Configs.weight_adjust:', Configs.weight_adjust)
    print('Configs.subset_size:', Configs.subset_size)
    print('Configs.num_cpus:', Configs.num_cpus)

    print('\nConfigs.hmmalign_path:', Configs.hmmalign_path)
    print('Configs.hmmsearch_path:', Configs.hmmsearch_path)
    print('Configs.magus_path:', Configs.magus_path)

    print('\nConfigs.log_path:', Configs.log_path)
    print('Configs.debug_path:', Configs.debug_path)
    print('Configs.error_path:', Configs.error_path)
    print('Configs.runtime_path:', Configs.runtime_path)
    
    print('\nConfigs.keepgcmtemp:', Configs.keepgcmtemp)
    print('Configs.inflation_factor:', Configs.inflation_factor)
    print('Configs.graphclustermethod:', Configs.graphclustermethod)
    print('Configs.graphtracemethod:', Configs.graphtracemethod)
    print('Configs.graphtraceoptimize:', Configs.graphtraceoptimize)

def buildConfigs(args):
    Configs.hmmdir = os.path.abspath(args.hmmdir)
    Configs.backbone_path = os.path.abspath(args.backbone_path)
    Configs.query_path = os.path.abspath(args.query_path)
    
    if args.outdir is not None:
        Configs.outdir = os.path.abspath(args.outdir)
    if not os.path.exists(Configs.outdir):
        os.makedirs(Configs.outdir)

    Configs.keeptemp = args.keeptemp
    Configs.keepsubalignment = args.keepsubalignment

    Configs.log_path = os.path.join(Configs.outdir, 'log.txt')
    Configs.error_path = os.path.join(Configs.outdir, 'error.txt')
    Configs.debug_path = os.path.join(Configs.outdir, 'debug.txt')
    Configs.runtime_path = os.path.join(Configs.outdir, 'runtime_breakdown.txt')

    if args.num_hmms > 0:
        Configs.num_hmms = args.num_hmms

    if args.use_weight:
        Configs.use_weight = True
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

    # MAGUS/GCM options
    Configs.keepgcmtemp = args.keepgcmtemp
    Configs.inflation_factor = args.inflation_factor
    Configs.graphclustermethod = args.graphclustermethod
    Configs.graphtracemethod = args.graphtracemethod
    Configs.graphtraceoptimize = args.graphtraceoptimize
