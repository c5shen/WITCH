import os, sys, time
import argparse
import logging

from configs import *

def main():
    start_time = time.time()
    args = parseArgs()
    buildConfigs(args)
    Configs.log('GCM+eHMMs is running with: {}'.format(' '.join(sys.argv)))

    getConfigs()


def parseArgs():
    parser = argparse.ArgumentParser()

    # requirements
    parser.add_argument('-p', '--hmmdir', type=str,
            help='Path to the HMMs directory generated from UPP', required=True)
    parser.add_argument('-b', '--backbone-path', type=str,
            help='Path to the backbone alignment from/for UPP', required=True)
    parser.add_argument('-o', '--outdir', type=str,
            help='Output directory, default: gcm+eHMMs_output/', required=False,
            default='gcm+eHMMs_output')

    # temporary files
    parser.add_argument('--keeptemp', action='store_const', const=True,
            help='Whether to keep temporary files for the run.',
            default=False)

    # options
    parser.add_argument('-k', '--num-hmms', type=int,
            help='The number of top HMMs used for aligning a query',
            required=False, default=4)
    parser.add_argument('-w', '--use-weight', action='store_const', const=True,
            help='Whether to use adjusted bitscore, default: 0', required=False,
            default=False)
    parser.add_argument('-s', '--subset-size', type=int,
            help='Number of queries in a single GCM run, default: 10',
            required=False, default=10)
    parser.add_argument('-t', '--num-threads', type=int,
            help='Number of threads for multi-threading (currently only supported for MAGUS/GCM)',
            required=False, default=1)

    return parser.parse_args()


if __name__ == '__main__':
    main()
