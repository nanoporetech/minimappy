Welcome to minimappy's documentation!
=====================================

Python bindings to `minimap2` aligner; sufficient to load and index and perform
alignments of sequences to the index to obtain basic statistics.

These python bindings are licensed under Mozilla Public License 2.0, minimap is licenced
under The MIT License.

Installation
------------

The git source repository contains minimap2 as a submodule. The repository should
therefore be cloned using the recursive option.

The package `setup.py` script requires `libminimap2.a` to have been built in the
submodule directory before running. This can be performed via the `libminimap2.a`
target. To build and install the package one should therefore run:

.. code-block:: bash

    git clone --recursive https://github.com/nanoporetech/minimappy.git
    make minimap2/libminimap2.a
    python setup.py install


Performing Alignments
---------------------

The `MinimapAligner` class provides a pythonic interface to `minimap2` aligner.
It takes as input minimap index or reference .fasta on construction and can
then be used to find alignments of sequences given as strings:

.. code-block:: python

    from minimappy import MinimapAligner
    index = 'path/to/index' # the path given to minimap index
    seq = 'ACGATCGCGATCGA'

    aligner = MinimapAlinger(index)
    alignments = aligner.align_seq(seq)
    print('Found {} alignments.'.format(len(alignments))
    for aln in alignments:
        print('  ', aln)

The alignments are returned as a named tuple, e.g.:

.. code-block:: python

    Alignment(rname='yeast', orient='+', pos=0, mapq=60, cigar='915M3D29M3D27M3D13M', NM=12, flag=0)


Alignment parameters can be given as they are on the `minimap2` command line:

.. code-block:: python

    from minimappy import MinimapAligner
    index = 'path/to/index'
    options = '-x ont2d -A 1 -B 0'
    aligner = MinimapAligner(index, options=options)


Contents
--------

.. toctree::
   :maxdepth: 2


Full API reference
------------------

.. toctree::
   :maxdepth: 3
      
   minimappy

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

