# WITCH - WeIghTed Consensus Hmm alignment
![visitors](https://visitor-badge.glitch.me/badge?page_id=c5shen.visitor-badge&left_color=blue&right_color=black)

(C) Chengze Shen

-----------------------------
Method Overview
-----------------------------
### Algorithm
WITCH is a new multiple sequence alignment (MSA) tool that combines techniques from [UPP](https://github.com/smirarab/sepp/blob/master/README.UPP.md) and [MAGUS](https://github.com/vlasmirnov/MAGUS). It aims to solve alignment problems particularly when input sequences contain fragments. The whole pipeline can be described as follows:
1. Given a set of unaligned sequences `S`, pick at most 1,000 "full-length" sequences to form a _backbone alignment_ `B` and a _backbone tree_ `T` (Full-length sequences refer to sequences of lengths that are within 25% of the median length).
2. Create an ensemble of HMMs (eHMM, see [UPP](https://github.com/smirarab/sepp/blob/master/README.UPP.md) for more details) from `B` and `T`.
3. For each remaining unaligned sequence, align it to high-ranked HMMs to obtain a set of weighted support alignments; then, merge the support alignments using Graph Clustering Merger (GCM, an alignment merger technique introduced in MAGUS).
4. Transitively add the merged alignment of each query to `B`, and report the final alignment on `S`.

For a more detailed explanation of the WITCH algorithm, please refer to the publication below:

#### Publication
1. Shen, Chengze, Minhyuk Park, and Tandy Warnow. “WITCH: Improved Multiple Sequence Alignment Through Weighted Consensus Hidden Markov Model Alignment.” Journal of Computational Biology, May 17, 2022. https://doi.org/10.1089/cmb.2021.0585.

### Note and Acknowledgement
- WITCH includes and uses:
    1. [MAGUS](https://github.com/vlasmirnov/MAGUS) (we use the Github version updated on April 5th 2021).
    2. [HMMER suites](http://hmmer.org/) (v3.3.2 - hmmbuild, hmmsearch, hmmalign).
    3. Partial functionalities from [UPP](https://github.com/smirarab/sepp/blob/master/README.UPP.md) (v4.5.1).
    4. [FastTreeMP](http://www.microbesonline.org/fasttree/FastTreeMP) (v2.1).


---------------------------
Installation
---------------------------
This section lays out necessary steps to do to run WITCH. We tested WITCH on the following system:
* Red Hat Enterprise Linux Server release 7.9 (Maipo) with **Python 3.7.0**
* Ubuntu 18.04.6 LTS with **Python 3.7.6**
* macOS Monterey 12.4 with **Python 3.9.13**

Now the program fully supports Linux and macOS systems (for at least the ones mentioned above). We provide necessary binary executables for both types of systems, but you can supplement your own by changing the paths in the `main.config` file.
If you experience any difficulty running WITCH, please contact Chengze Shen (chengze5@illinois.edu).

### Requirements
```
configparser==5.0.0
DendroPy==4.4.0
numpy==1.21.0
psutil==5.6.7
```

### Installation Steps
```bash
# 1. install via GitHub repo
git clone https://github.com/c5shen/WITCH.git

# 2. install all requirements
cd WITCH
pip3 install -r requirements.txt

# 3. run setup.py to set up main.config
python3 setup.py

# 4. execute the WITCH python script with -h
python3 witch.py -h
```

----------------------------
Usage
----------------------------
General command to run WITCH:
```
python3 witch.py -i <unaligned sequence file> -d <output directory> -o <output filename>
```
**Default behavior**: WITCH will pick at most 1,000 sequences from the input around the median length as the backbone sequences. Then, it uses MAGUS to align the backbone sequences and FastTree2 to estimate a tree. It uses UPP decomposition strategy to generate an eHMM, and uses HMMSearch to calculate bit scores between HMMs and unaligned sequences. Bit scores are used to calculate weights, and each unaligned sequence is aligned to top `k=4` HMMs ranked by weights.

#### Use regular bit scores
By default, WITCH uses HMMSearch to obtain bit scores, and then uses bit scores to calculate weights between unaligned sequences and HMMs. To use bit scores instead of weights, run WITCH by the following command:
```bash
python3 witch.py -w 0 [...other parameters...]
```

#### Multi-processing
By default, WITCH uses all available cores on the machine. Users can choose the number of cores by the following command:
```bash
python3 witch.py -t <number of cpus> [...other parameters...]
```

To obtain the full list of parameters and options, please use `python3 witch.py -h` or `python3 witch.py --help`.

-------------------------
Examples
-------------------------
All the following examples can be found in the **examples/run.sh** bash script.
### Scenario A - unaligned sequences only
```bash
python3 witch.py -i examples/data/unaligned_all.txt -d scenarioA_output -o aligned.txt
```

### Scenario B - unaligned sequences only; using bit scores; using 10 HMMs to align a sequence
```bash
python3 witch.py -i examples/data/unaligned_all.txt -d scenarioB_output -o aligned.txt -w 0 -k 10
```

### Scenario C - backbone alignment available; backbone tree missing; query sequences available
```bash
python3 witch.py -b examples/data/backbone.aln.fasta -q examples/data/unaligned_frag.txt -d scenarioC_output -o aligned.txt
```

-------------------------
To-do
-------------------------
1. Optimize how merging is done to reduce the output file size (need to do a similar task as UPP where insertion columns are marked as lower cases and squeezed together.
