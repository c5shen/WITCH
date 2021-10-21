-----------------------------
Summary
-----------------------------
### Overview
GCM+eHMMs is a new multiple sequence alignment (MSA) tool that combines techniques from [UPP](https://github.com/smirarab/sepp/blob/master/README.UPP.md) and [MAGUS](https://github.com/vlasmirnov/MAGUS). It aims to solve alignment problems when input sequences contain fragments or show high levels of sequence length hetereogeneity. The whole pipeline can be described as follows:
1. Given a set of unaligned sequences `S`, pick at most 1,000 full-length sequences to form a _backbone alignment_ `B` and a _backbone tree_ `T`. Full-length sequences refer to sequences of lengths that are within 25% of the median length.
2. Create an ensemble of HMMs (see [UPP](https://github.com/smirarab/sepp/blob/master/README.UPP.md) for more details) from `B` and `T` with decomposition subset size set to `z=2` (`z=10` in default UPP).
3. For each remaining unaligned sequence, align it to top `k` HMMs ranked by adjusted bit scores (i.e., weights) on the query sequence, and merge the `k` support alignments using Graph Clustering Merger with weights (GCM, an alignment merger technique introduced in MAGUS. We modified GCM to allow it to take a weight for each support alignment).
4. Transitively add the merged alignment of each query to `B`, and report the final alignment on `S`.

**Potential speedup on large datasets**: on large datasets, Step 2 considerably slows down because of the number of searches required between all query sequences and all HMMs. Instead of using `z=2`, users can choose to set `z=10`, which will reduce the number of HMMs created considerably (thus speed up searches) at a minor to none cost on alignment accuracy. Final alignment might also be longer (slightly higher under-alignment issue).

**(Oct 21st 2021) Current implementation of GCM+eHMMs only supports the following problem (we are working on intergrating UPP to GCM+eHMMs):**
* Input: a set of sequences potentially with fragmentary sequences.
* Intermediate output (**these are inputs to the current GCM+eHMMs**):
    1. a backbone alignment/tree on at most 1,000 full-length sequences (see the overview above);
    2. a set of HMMs generated using UPP tree decomposition.
* Output: an alignment on all sequences.

### Developers
Chengze Shen

### Publication
PLACEHOLDER

### Note and Acknowledgement
- GCM+eHMMs includes and uses:
    1. [MAGUS](https://github.com/vlasmirnov/MAGUS) (we use the Github version updated on April 5th 2021).
    2. [HMMER](http://hmmer.org/) suites (hmmalign v3.3.2, hmmsearch v3.1b2).
    3. (PLANNED) [UPP v4.5.1](https://github.com/smirarab/sepp/blob/master/README.UPP.md).
    4. [Dendropy](https://dendropy.org/) Python package.

---------------------------
Installation
---------------------------
This section lays out necessary steps to do to run GCM+eHMMs. We tested GCM+eHMMs on Linux. If you experience any difficulty running GCM+eHMMs, please contact Chengze Shen at chengze5@illinois.edu.

### Requirements
```
Python >=3.6
Dendropy >= 4.4.0
(PLANNED) UPP >= 4.5.1
```

### Installation Steps
1. Make sure all requirements are met and `run_upp.py` executable is on the PATH variable.
2. Execute the python file as `python gcm+eHMMs.py`.

----------------------------
Ruuning GCM+eHMMs
----------------------------
General command to run GCM+eHMMs is:
```
python gcm+eHMMs.py -p <HMM directory> -b <backbone alignment> -q <query sequences> -o <output directory> -w
```
`-p` denotes the directory where the ensemble of HMMs are stored (e.g., from a previous UPP run), `-b` denotes the path to the backbone alignment, `-q` denotes the path to unaligned sequences, `-w` tells the software that we want to use weights to rank HMMs for queries, and `-o` denotes the output directory. After GCM+eHMMs finishes, a file named `merged.fasta` that represents the final alignment will appear in the output directory.

#### Use regular bit scores
To use regular bit scores instead of weights to rank HMMs for queries, simply remove `-w` from the command.

#### Multi-processing
By default GCM+eHMMs uses all available cores to run. Users can choose the number of cores by the following command:
```
python gcm+eHMMs.py -t <number of cpus> [...other parameters...]
```


To obtain the full list of parameters and options, please use `python gcm+eHMMs.py -h` or `python gcm+eHMMs.py --help`:
```bash
REQUIRED PARAMETERS:
  These are required fields, assuming that an UPP output and its temporary
  files are available.

  -p HMMDIR, --hmmdir HMMDIR
                        Path to the HMMs directory generated from UPP
  -b BACKBONE_PATH, --backbone-path BACKBONE_PATH
                        Path to the backbone alignment from/for UPP
  -q QUERY_PATH, --query-path QUERY_PATH
                        Path to the queries file that we want to align
  -o OUTDIR, --outdir OUTDIR
                        Output directory, default: gcm+eHMMs_output/

GCM+EHMMS OPTIONS:
  These options are used to customize the GCM+eHMMs pipeline.

  --keeptemp            Keep ALL temporary files in the process (constraints,
                        backbones, HMMSearch results, GCM results, etc.)
  --keepsubalignment    Keep all subalignments by MAGUS/GCM
  -k NUM_HMMS, --num-hmms NUM_HMMS
                        The number of top HMMs used for aligning a query
  -w, --use-weight      Whether to use adjusted bitscore (weights), default: 0
  -s SUBSET_SIZE, --subset-size SUBSET_SIZE
                        Number of queries in a single GCM run, default: 1
  --weight-adjust {none,normalize,maxto1}
                        Optional adjustment of weights, default: none
  -t NUM_CPUS, --num-cpus NUM_CPUS
                        Number of cpus for multi-processing, default: -1 (all)

MAGUS/GCM OPTIONS:
  These options are used to customize MAGUS/GCM, for example the graph trace
  method.

  --keepgcmtemp         Keep temporary files generated from MAGUS/GCM
  --timeout TIMEOUT     Retry a MAGUS/GCM subtask after [timeout] seconds,
                        default: 60
  -f INFLATION_FACTOR, --inflation-factor INFLATION_FACTOR
                        Inflation factor for MCL, default: 4
  --graphclustermethod {mcl,mlrmcl,rg,none}
                        Method for initial clustering of the alignment graph,
                        default: mcl
  --graphtracemethod {minclusters,mwtgreedy,mwtsearch,fm,rg,rgfast}
                        Method for finding a trace from the alignment graph,
                        default: minclusters
  --graphtraceoptimize {true,false}
                        Run an optimization step on the graph trace, default:
                        false
```
