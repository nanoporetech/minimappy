minimappy
=========

Status
------

Since the time of writing Heng Li, the author of minimap2, has provided a
python binding in the minimap2 source and at https://pypi.python.org/pypi/mappy.
Users are strongly encouraged to use Heng's python library rather that
the version presented here. Despite its short life, the latter will be no
longer supported out of deference to Heng's version.


[![Build Status](https://travis-ci.org/nanoporetech/minimappy.svg?branch=master)](https://travis-ci.org/nanoporetech/minimappy)

Python bindings to `minimap2` aligner; sufficient to load and index and perform
alignments of sequences to the index to obtain basic statistics.

These python bindings are licensed under Mozilla Public License 2.0, minimap is
licenced The MIT License.

Documentation can be found at https://nanoporetech.github.io/minimappy/.

Installation
------------

The git source repository contains minimap2 as a submodule. The repository
should therefore be cloned using the recursive option.

The package `setup.py` script requires `libminimap2.a` to have been built in the
submodule directory before running. This can be performed via the `libminimap2.a`
target. To build and install the package one should
therefore run:

    git clone --recursive https://github.com/nanoporetech/minimappy.git
    make minimap2/libminimap2.a 
    python setup.py install


Performing Alignments
---------------------

The `MinimapAligner` class provides a pythonic interface to `minimap2` aligner.
It takes as input a minimap index or .fasta reference on construction and can
then be used to find alignments of sequences given as strings:

```python
from minimappy import MinimapAligner
index = 'path/to/index' # the path given to minimap2 index
seq = 'ACGATCGCGATCGA'

aligner = MinimapAligner(index)
alignments = aligner.align_seq(seq)
print('Found {} alignments.'.format(len(alignments))
for aln in alignments:
    print('  ', aln)
```

The alignments are returned as a named tuple, e.g.:

```python
Alignment(rname='yeast', orient='+', pos=0, mapq=60, cigar='915M3D29M3D27M3D13M', NM=12, flags=0)
```

Alignment parameters can be given as they are on the `minimap2` command line:

```python
from minimappy import MinimapAligner
index = 'path/to/index'
options = '-ax map10k'
aligner = MinimapAligner(index, options=options)
```

