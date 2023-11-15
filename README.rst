WITCH - WeIghTed Consensus Hmm alignment
========================================

|PyPI version fury.io| |PyPI pyversions| |PyPI license| |DOI|

Developer: Chengze Shen, Baqiao Liu

.. contents:: Table of contents
   :backlinks: top
   :local:

News
----
* (NEW) Automatically infer data type if None is specified (use ``--molecule`` to specify).
* (NEW) Checkpoint system set up for most steps except HMMSearch jobs (ongoing).
* (NEW) Added progress bar (python package ``tqdm``) to visualize the alignment progress at various stages.
* Implemented `WITCH-ng <https://github.com/RuneBlaze/WITCH-NG>`__’s way to align each query sequence with additional tweaks. Now the alignment process for query sequences is **fast and memory-efficient, particularly for short/fragmentary sequences**.


TODO list
---------
#. (Priority) Setting up a checkpoint system for HMMSearch jobs.


Method Overview
---------------
WITCH is a new multiple sequence alignment (MSA) tool that combines techniques from `UPP <https://github.com/smirarab/sepp/blob/master/README.UPP.md>`__ and `MAGUS <https://github.com/vlasmirnov/MAGUS>`__.
It aims to solve alignment problems, particularly when input sequences contain fragments. The whole pipeline can be described as follows:

#. Given a set of unaligned sequences ``S``, pick at most 1,000 “full-length” sequences to form a *backbone alignment* ``B`` and a *backbone tree* ``T``
   (Full-length sequences refer to sequences of lengths that are within 25% of the median length).
#. Create an ensemble of HMMs (eHMM, see `UPP <https://github.com/smirarab/sepp/blob/master/README.UPP.md>`__ for more details) from ``B`` and ``T``. 
#. For each remaining unaligned sequence, align it to high-ranked HMMs to obtain a set of weighted support alignments; then, merge the support alignments using Graph Clustering Merger
   (GCM, an alignment merger technique introduced in MAGUS). 4. Transitively add the merged alignment of  each query to ``B``, and report the final alignment on ``S``.

.. image:: https://chengzeshen.com/documents/gcm_ehmm/pipeline.png
   :alt: WITCH pipeline
   :width: 70%
   :align: center

For a more detailed explanation of the WITCH algorithm, please refer to the publication below:

+----------------------------------------+
| Publication                            |
+========================================+
| Shen,                                  |  
| Chengze, Minhyuk Park, and             |
| Tandy Warnow. “WITCH:                  |
| Improved Multiple Sequence             |
| Alignment Through Weighted             |
| Consensus Hidden Markov                |
| Model Alignment.” Journal              |
| of Computational Biology,              |
| May 17, 2022.                          |
| https://doi.org/10.1089/cmb.2021.0585. |
+----------------------------------------+

Note and Acknowledgement
~~~~~~~~~~~~~~~~~~~~~~~~
WITCH includes and uses:

#. `MAGUS <https://github.com/vlasmirnov/MAGUS>`__ (we use the Github version updated on April 5th 2021).
#. `HMMER suites <http://hmmer.org/>`__ (v3.1b2 - hmmbuild, hmmsearch, hmmalign).
#. `UPP <https://github.com/smirarab/sepp/blob/master/README.UPP.md>`__ (v4.5.1; we use only partial functionalities).
#. `FastTreeMP <http://www.microbesonline.org/fasttree/FastTreeMP>`__ (v2.1). 
#. `MAFFT <https://mafft.cbrc.jp/alignment/software/macportable.html>`__ (macOS v7.490).
#. `MCL <https://github.com/micans/mcl>`__ (linux version from MAGUS; macOS version 21-257).

Installation
------------

This section lays out the necessary steps to run WITCH. We tested
WITCH on the following systems: \* Red Hat Enterprise Linux Server
release 7.9 (Maipo) with **Python 3.7.0** \* Ubuntu 18.04.6 LTS with
**Python 3.7.6**, and Ubuntu 22.04 LTS with **Python 3.7.12** \* macOS
*(x86 chip)* Monterey 12.4 with **Python 3.9.13**.

Now, the program fully supports Linux and macOS systems (for at least the
ones mentioned above). We provide necessary binary executables for both
types of systems, but you can supplement your own by changing the paths
in the ``main.config`` file. In cases of conflicting installations
(e.g., different versions of MAFFT), please supplement with the version
on your system. If you experience any difficulty running WITCH, please
contact Chengze Shen (chengze5@illinois.edu).

   For the macOS system on the latest chips (e.g., M1/M2), you may need to compile and supply your own binaries for WITCH to run successfully.
   That is, change the paths of binaries in ``main.config`` to the ones on your system.

Python version (REQUIRED!)
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   python>=3.7

Requirements
~~~~~~~~~~~~

::

   cython>=0.29
   configparser>=5.0.0
   DendroPy>=4.4.0,<4.6.0
   numpy>=1.15
   psutil>=5.0
   tqdm>=4.0.0

Installation Steps
~~~~~~~~~~~~~~~~~~

.. code:: bash

   # 1. Install via GitHub repo
   git clone https://github.com/c5shen/WITCH.git

   # 2. Install all requirements
   # If you do not have root access, use "pip3 install -r requirements.txt --user"
   cd WITCH
   pip3 install -r requirements.txt

   # 3. Run setup.py to set up main.config. Please refer to default.config and use `-h` for additional information
   #    Additionally, software binaries that are available in the user's environment will be prioritized for usage.
   #    Use "-p false" to disable this priority.
   python3 setup.py [-h]

   # 4. Execute the WITCH python script with -h to see allowed commandline parameter settings
   #    When running WITCH normally, if step 3 is not run, WITCH will automatically generate a "main.config" file
   #    using the default "setup.py" settings.
   python3 witch.py [-h]

Usage
-----
The general command to run WITCH:

.. code:: bash

    python3 witch.py -i [unaligned sequence file] -d [output directory] -o [output filename]

**Default behavior**: WITCH will pick at most 1,000 sequences from the input around the median length as the backbone sequences.
Then, it uses MAGUS to align the backbone sequences and FastTree2 to estimate a tree. It uses UPP decomposition strategy to generate an eHMM,
and uses HMMSearch to calculate bit scores between HMMs and unaligned sequences. Bit scores are used to calculate weights, and each unaligned
sequence is aligned to top `k=10` HMMs ranked by weights.


Examples
~~~~~~~~

All the following examples can be found in the **examples/run.sh** bash
script.

Scenario A
++++++++++
Unaligned sequences only.

.. code:: bash

   python3 witch.py -i examples/data/unaligned_all.txt \
      -d scenarioA_output -o aligned.txt

Scenario B
++++++++++
Unaligned sequences only; using bit scores (instead of the default weighted bit scores); using 10 HMMs to align a sequence.

.. code:: bash

   python3 witch.py -i examples/data/unaligned_all.txt \
      -d scenarioB_output -o aligned.txt -w 0 -k 10

Scenario C
++++++++++
Backbone alignment available; backbone tree missing; query sequences available.

.. code:: bash

   python3 witch.py -b examples/data/backbone.aln.fasta \
      -q examples/data/unaligned_frag.txt -d scenarioC_output \
      -o aligned.txt

Scenario D - additional options
+++++++++++++++++++++++++++++++
Backbone alignment available; backbone tree available; query sequences available; saving weights to local; saving decomposition results for future usage (e.g., faster rerun).

.. code:: bash

   python3 witch.py -b examples/data/backbone.aln.fasta \
      -e examples/data/backbone.tre -q examples/data/unaligned_frag.txt \
      -d scenarioD_output -o aligned.txt \
      --save-weight 1 --keep-decomposition 1

.. |publication| image:: https://img.shields.io/badge/Publication-Journal_of_Computational_Biology-green?style=for-the-badge
   :target: https://doi.org/10.1089/cmb.2021.0585
.. |coverage| image:: https://coveralls.io/repos/github/c5shen/WITCH/badge.svg?branch=main
   :target: https://coveralls.io/github/c5shen/WITCH?branch=main
.. |PyPI version fury.io| image:: https://badge.fury.io/py/witch-msa.svg
   :target: https://pypi.python.org/pypi/witch-msa/
.. |PyPI license| image:: https://img.shields.io/pypi/l/witch-msa.svg
   :target: https://pypi.python.org/pypi/witch-msa/
.. |PyPI pyversions| image:: https://img.shields.io/pypi/pyversions/witch-msa.svg
   :target: https://pypi.python.org/pypi/witch-msa/
.. |DOI| image:: https://zenodo.org/badge/DOI/10.1089/cmb.2021.0585.svg
   :target: https://doi.org/10.1089/cmb.2021.0585
