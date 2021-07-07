import os
import time

class Configs:
    hmmdir = None
    backbone_path = None
    outdir = None
    keeptemp = False
    
    num_hmms = 4
    use_weight = False
    subset_size = 10
    num_threads = 1

    # hmmalign/hmmsearch/magus paths
    hmmalign_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools/hmmer/hmmalign')
    hmmsearch_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools/hmmer/hmmsearch')
    magus_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools/magus/magus.py')

    log_path = None
    error_path = None

    @staticmethod
    def warning(msg, path=None):
        path = Configs.log_path if path is None else path
        Configs.write(msg, 'WARNING', path)

    @staticmethod
    def log(msg, path=None):
        path = Configs.log_path if path is None else path
        Configs.write(msg, 'LOG', path)

    @staticmethod
    def error(msg, path=None):
        path = Configs.error_path if path is None else path
        Configs.write(msg, 'ERROR', path)

    @staticmethod
    def write(msg, level, path):
        if path is not None:
            with open(path, 'a') as f:
                f.write('{}\t[{}]{}\n'.format(time.strftime('%Y-%m-%d %H:%M:%S'),
                    level, msg))

# print a list of all configurations
def getConfigs():
    print('Configs.hmmdir:', Configs.hmmdir)
    print('Configs.backbone_path:', Configs.backbone_path)
    print('Configs.outdir:', Configs.outdir)
    print('Configs.keeptemp:', Configs.keeptemp)

    print('\nConfigs.num_hmms:', Configs.num_hmms)
    print('Configs.use_weight:', Configs.use_weight)
    print('Configs.subset_size:', Configs.subset_size)
    print('Configs.num_threads:', Configs.num_threads)

    print('\nConfigs.hmmalign_path:', Configs.hmmalign_path)
    print('Configs.hmmsearch_path:', Configs.hmmsearch_path)
    print('Configs.magus_path:', Configs.magus_path)

    print('\nConfigs.log_path:', Configs.log_path)
    print('Configs.error_path:', Configs.error_path)

def buildConfigs(args):
    Configs.hmmdir = os.path.abspath(args.hmmdir)
    Configs.backbone_path = os.path.abspath(args.backbone_path)
    
    if args.outdir is not None:
        Configs.outdir = os.path.abspath(args.outdir)
    if not os.path.exists(Configs.outdir):
        os.makedirs(Configs.outdir)

    Configs.log_path = os.path.join(Configs.outdir, 'log.info')
    Configs.error_path = os.path.join(Configs.outdir, 'error.info')

    if args.num_hmms > 0:
        Configs.num_hmms = args.num_hmms

    if args.use_weight:
        Configs.use_weight = True

    if args.subset_size > 0:
        if args.subset_size >= 25:
            Configs.warning('Subset size is recommended to be <= 25! Current: {}'.format(
                args.subset_size))
        Configs.subset_size = args.subset_size

    if args.num_threads > 0:
        Configs.num_threads = args.num_threads
    else:
        Configs.num_threads = os.cpu_count()
