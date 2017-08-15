import argparse
from collections import namedtuple
import importlib
import imp
import os
import sys

from cffi import FFI
ffi = FFI()

"""High-level interface to minimap2 aligner."""


def get_shared_lib(name):
    """Cross-platform resolution of shared-object libraries, working
    around vagueries of setuptools.

    :param name: name of shared library to find.
    
    :returns: FFI shared library object.
    """
    try:
        # after 'python setup.py install' we should be able to do this
        lib_file = importlib.import_module(name).__file__
    except Exception as e:
        try:
            # after 'python setup.py develop' this should work
            lib_file = imp.find_module(name)[1]
        except Exception as e:
            raise ImportError('Cannot locate C library "{}".'.format(name))
        else:
            lib_file = os.path.abspath(lib_file)
    finally:
        library = ffi.dlopen(lib_file)
    return library


libmap = get_shared_lib('minimaplib')

ffi.cdef("""
  ////////////////
  // Mapping Index
  //
  typedef struct {
    uint64_t x, y;
  } mm128_t;
  
  typedef struct { size_t n, m; mm128_t *a; } mm128_v;
  
  typedef struct {
    mm128_v a;   // (minimizer, position) array
    int32_t n;   // size of the _p_ array
    uint64_t *p; // position array for minimizers appearing >1 times
    void *h;     // hash table indexing _p_ and minimizers appearing once
  } mm_idx_bucket_t;
  
  typedef struct {
    char *name;      // name of the db sequence
    uint64_t offset; // offset in mm_idx_t::S
    uint32_t len;    // length
  } mm_idx_seq_t;
  
  typedef struct {
    int32_t b, w, k, is_hpc;
    uint32_t n_seq;     // number of reference sequences
    mm_idx_seq_t *seq;  // sequence name, length and offset
    uint32_t *S;        // 4-bit packed sequence
    mm_idx_bucket_t *B; // index
    void *km;
  } mm_idx_t;

  mm_idx_t* get_index(const char* fname);
  void destroy_index(mm_idx_t* index);

  ////////////
  // Alignment
  //
  typedef struct {
    uint32_t capacity;
    int32_t dp_score, dp_max, dp_max2;
    uint32_t blen;
    uint32_t n_diff, n_ambi;
    uint32_t n_cigar;
    uint32_t cigar[];
  } mm_extra_t;

  typedef struct {
    int32_t id;
    uint32_t cnt:31, rev:1;
    uint32_t rid:31, inv:1;
    int32_t score;
    int32_t qs, qe, rs, re;
    int32_t parent, subsc;
    int32_t as;
    int32_t fuzzy_mlen, fuzzy_blen;
    uint32_t mapq:8, split:2, sam_pri:1, n_sub:21; // TODO: n_sub is not used for now
    mm_extra_t *p;
  } mm_reg1_t;

  typedef struct {size_t n; mm_reg1_t* reg;} mm_reg1_v;

  mm_reg1_v align(mm_idx_t* index, char* seq, char* seq_name);
""")


Alignment = namedtuple('Alignment', [
    'rname', 'orient', 'pos', 'mapq', 'cigar', 'NM', 'flag'
])

class MinimapAligner(object):
    def __init__(self, index:str, options:str=''):
        """Interface to minimap2 alignment.

        :param index: minimap index path.
        :param options: alignment options as would be given
            on the minimap command line.

        """
        self.index_file = index.encode()
        self._cigchar = "MIDSH"

        if options != '':
            raise NotImplementedError()

        self.index = libmap.get_index(self.index_file)
        if self.index == ffi.NULL:
            raise ValueError('Failed to load minimap index.')


    def __del__(self):
        if hasattr(self, 'index'):
            libmap.destroy_index(self.index)
        if hasattr(self, 'opts'):
            cffi.free(self.opts)


    def _build_alignment(self, reg, seq, reg0):
        flag = 0
        if reg.rev: flag |= 0x10
        if reg.parent != reg.id:
            flag |= 0x100
        elif not reg.sam_pri:
            flag |= 0x800

        l_seq = len(seq)
        # cigar clipping
        clip_char = 'H' if flag&0x800 else 'S'
        clip_start = l_seq - reg.qe if reg.rev else reg.qs
        if clip_start > 0:
            clip_start = '{}{}'.format(clip_start, clip_char)
        else:
           clip_start = ''
        clip_end = reg.qs if reg.rev else l_seq - reg.qe
        if clip_end > 0:
            clip_end = '{}{}'.format(clip_end, clip_char)
        else:
            clip_end = ''
        # cigar in aligned region
        cig = reg.p.cigar
        cigar = ''.join(
            # oplen + op
            str(cig[k]>>4) + self._cigchar[cig[k] & 0xf]
            for k in range(reg.p.n_cigar)
        )
        cigar = '{}{}{}'.format(clip_start, cigar, clip_end)

        return Alignment(
            ffi.string(self.index.seq[reg.rid].name).decode(),
            '+-'[reg.rev], reg.rs, reg.mapq, cigar, reg.p.n_diff, flag
        )

    def align_seq(self, seq:str, seq_name:str='no_name'):
        """Align a sequence to the index.

        :param seq: base sequence to align

        :returns: tuple of :class:`Alignment`
        """
        data = libmap.align(self.index, seq.encode(), seq_name.encode())
        res = [self._build_alignment(data.reg[i], seq, data.reg[0]) for i in range(data.n)]
        return res

def get_parser():
    parser = argparse.ArgumentParser('Align a sequence with minimap2.')
    parser.add_argument('index', help='minimap index path.')
    parser.add_argument('sequence', nargs='+', help='base sequence')
    return parser


def main():
    args, opts = get_parser().parse_known_args()
    options = ''
    if len(opts) > 0:
        options = ' '.join(opts)
    aligner = MinimapAligner(args.index, options=options)
    for i, seq in enumerate(args.sequence, 1):
        alignments = aligner.align_seq(seq)
        print('Found {} alignments for input {}.'.format(len(alignments), i))
        for aln in alignments:
            print('  ', aln)
 
