WITCH v1.0.10
-------------
1. Added an adaptive schematic for inclusion of top HMMs for aligning each
   query sequence, when using adjusted bitscores (``--use-weights 1``).
   Previously, all top `k` HMMs will be used. Now, WITCH includes up to `k`
   HMMs, or until the sum of weights exceeds 0.999. This should keep the
   core design of WITCH, but this has not been tested with data yet.

WITCH v1.0.9
------------
1. Fixed the issue with feeding FastTree2 with gzipped alignment file for
   tree estimation. Now use piping to pipe the alignment file (gzipped or
   not) as stdin to the FastTree2 executable.
2. Enforced the start method for multiprocessing on macOS to be ``fork``.
   This ensures WITCH usability on a macOS environment.
3. Changed the default invocation of example codes from ``python3`` to
   just ``python``.

WITCH v1.0.8
------------
1. Supported reading gzipped alignment file. That is, when the user supplies
   with their own backbone alignment and adds some other query sequences).
   E.g., ``backbone.fasta.gz`` or ``backbone.fa.gzip``, etc.

WITCH v1.0.7
------------
1. Added example usages to ``witch.py --help``. Also changed the default
   formatter from a custom one to ``argparse.RawDescriptionHelpFormatter``.
2. Now by default will bypass the initial WITCH setup step. Previously, this
   was achieved by giving WITCH the parameter ``-y`` (``--bypass-setup``).
3. Changed the default behavior of ``examples/run.sh`` to running scenario D.
4. Fixed the runtime error when home.path exists but the pointed main.config is
   missing. Now will regenerate ``main.config`` at the pointed location.
5. Fixed ``--bypass-setup`` not working as True by default. Now will always
   create the config path at ``~/.witch/main.config``. 
6. Changed the default filename for the masked alignment output. Previously,
   it will be named as ``<name>.fasta.masked``, as ``<name>`` supplied by the user
   or default to ``aligned``. Now will be written as ``<name>.masked.fasta``.
   If the user gives ``-o <name>.fa`` or ``-o <name>.fasta``, the masked alignment
   will use the corresponding suffix (e.g., ``<name>.masked.fa``).

WITCH v1.0.6
------------
1. Added Software Output Explanation to the README to avoid confusion on what
   alignment file to use for downstream analyses.

WITCH v1.0.5
------------
1. Added compatibility to Dendropy with version >4.5.2 and removed its
   requirement from requirements.txt for pip.

WITCH v1.0.5b
-------------
1. Added a new parameter option allowing users to specify a customized config
   file to override the default ``main.config`` (usually can be found at
   ``~/.witch_msa/main.config``). Use ``-c <user config file>`` to do so.
   The priority for arguments: ``commandline > user.config > main.config``.

WITCH v1.0.5a
-------------
1. Added two sanity checks to HMMBuild and HMMSearch jobs: making sure all
   files are created correctly before proceeding.
2. Added a file number check utility function using the ``inspect`` package.

WITCH v1.0.4
------------
1. Added an additional parameter option to set an upper bound to the HMM
   subsets created (``-Z``, complementary to ``-A`` which is for lower bound),
   based on the number of sequences in a subset.
2. Changed the behavior for creating HMM subsets. Instead of reading in the
   backboen alignment at once, WITCH now reads line by line to avoid large
   memory consumption if the backbone is very large.
3. Changed the behavior for running HMMSearch. Now also uses my task manager
   to manage the maximum amount of concurrently running jobs (instead of
   submitting all at once).

WITCH v1.0.3
------------
1. Fixed an oversight in memory issue if the backbone alignment is large and
   we are creating a lot of subsets (``gcmm/algorithm.py``).
2. Added default behavior to function ``gcmm/task.py/handleFuture(...)`` so that
   return values are attached to ``success`` if no callback function is provided.
3. Moving towards using my generic task manager for all parallelization.

WITCH v1.0.2
------------
1. Added an option to bypass the initial config setup step. Use ``-y``, or
   ``--bypass-setup`` in your commandline running WITCH to avoid being asked where
   to generate the config directory (will default to ``~/.witch_msa``). Example:
   running ``witch.py -y -i [sequence files]`` for the first time will directly
   set up WITCH configuration file and start aligning the input sequences.

WITCH v1.0.2a1
--------------
1. Fixed not using absolute path when setting up the directory for
   ``main.config``.
2. Improved metavar naming for argument parameters.
3. Fixed a bug in the included MAGUS installation such that if the user has 
   a file or directory named ``fasttree/`` in the folder where WITCH is run,
   MAGUS will try to read it as the guide tree instead of creating a FastTree
   guide tree.
4. Changed the executable search from ``usr/bin/env python`` to
   ``usr/bin/env python3``.
5. Further added systemrecursionlimit from 10000 to 20000 to combat issues
   with large tree.

WITCH v1.0.1
------------
1. First working release across different platform and different python
   versions.
