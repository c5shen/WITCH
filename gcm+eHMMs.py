import os, sys, time
import argparse

def main():
    start_time = time.time()
    args = parseArgs()
    print(args.weight)


def parseArgs():
    parser = argparse.ArgumentParser()

    parser.add_argument('-p', '--hmmdir', type=str,
            help='Path to the HMMs directory generated from UPP', required=True)
    parser.add_argument('-b', '--backbonepath', type=str,
            help='Path to the backbone alignment from/for UPP', required=True)
    parser.add_argument('-o' '--outdir', type=str,
            help='Output directory, default: gcm+eHMMs_output/', required=False,
            default='gcm+eHMMs_output')
    parser.add_argument('-k', '--numhmms', type=int,
            help='The number of top HMMs used for aligning a query',
            required=False, default=4)
    parser.add_argument('-w', '--weight', type=bool,
            help='Whether to use adjusted bitscore, default: 0', required=False,
            default=0)
    #parser.add_argument('--temp', type=bool)

    return parser.parse_args()


if __name__ == '__main__':
    main()
