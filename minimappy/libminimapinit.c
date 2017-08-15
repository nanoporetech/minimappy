#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"


static struct option long_options[] = {
	{ "bucket-bits",    required_argument, 0, 0 },
	{ "mb-size",        required_argument, 0, 'K' },
	{ "int-rname",      no_argument,       0, 0 },
	{ "no-kalloc",      no_argument,       0, 0 },
	{ "print-qname",    no_argument,       0, 0 },
	{ "no-self",        no_argument,       0, 0 },
	{ "print-seed",     no_argument,       0, 0 },
	{ "max-chain-skip", required_argument, 0, 0 },
	{ "min-dp-len",     required_argument, 0, 0 },
	{ "print-aln-seq",  no_argument,       0, 0 },
	{ "version",        no_argument,       0, 'V' },
	{ "min-count",      required_argument, 0, 'n' },
	{ "min-chain-score",required_argument, 0, 'm' },
	{ "mask-level",     required_argument, 0, 'M' },
	{ "min-dp-score",   required_argument, 0, 's' },
	{ "sam",            no_argument,       0, 'a' },
	{ 0, 0, 0, 0}
};

typedef struct {mm_mapopt_t opt; mm_idx_t *mi;} mm_opt_ind; 

mm_opt_ind init_opt_ind(int argc, char *argv[])
{
	mm_mapopt_t opt;
	int c, k = 15, w = -1, bucket_bits = MM_IDX_DEF_B, n_threads = 3, keep_name = 1, is_idx, is_hpc = 0, long_idx, idx_par_set = 0;
	int minibatch_size = 200000000;
	uint64_t batch_size = 4000000000ULL;
	char *fnw = 0, *s;

	mm_mapopt_init(&opt);

	while ((c = getopt_long(argc, argv, "aw:k:K:t:r:f:Vv:g:I:d:XT:s:x:Hcp:M:n:z:A:B:O:E:m:N:Q", long_options, &long_idx)) >= 0) {
		if (c == 'w') w = atoi(optarg), idx_par_set = 1;
		else if (c == 'k') k = atoi(optarg), idx_par_set = 1;
		else if (c == 'H') is_hpc = 1, idx_par_set = 1;
		else if (c == 'd') fnw = optarg; // the above are indexing related options, except -I
		else if (c == 'r') opt.bw = atoi(optarg);
		else if (c == 'f') opt.mid_occ_frac = atof(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'v') mm_verbose = atoi(optarg);
		else if (c == 'g') opt.max_gap = atoi(optarg);
		else if (c == 'N') opt.best_n = atoi(optarg);
		else if (c == 'p') opt.pri_ratio = atof(optarg);
		else if (c == 'M') opt.mask_level = atof(optarg);
		else if (c == 'c') opt.flag |= MM_F_CIGAR;
		else if (c == 'X') opt.flag |= MM_F_AVA | MM_F_NO_SELF;
		else if (c == 'a') opt.flag |= MM_F_OUT_SAM | MM_F_CIGAR;
		else if (c == 'Q') opt.flag |= MM_F_NO_QUAL;
		else if (c == 'T') opt.sdust_thres = atoi(optarg);
		else if (c == 'n') opt.min_cnt = atoi(optarg);
		else if (c == 'm') opt.min_chain_score = atoi(optarg);
		else if (c == 'A') opt.a = atoi(optarg);
		else if (c == 'B') opt.b = atoi(optarg);
		else if (c == 'z') opt.zdrop = atoi(optarg);
		else if (c == 's') opt.min_dp_max = atoi(optarg);
		else if (c == 0 && long_idx == 0) bucket_bits = atoi(optarg); // --bucket-bits
		else if (c == 0 && long_idx == 2) keep_name = 0; // --int-rname
		else if (c == 0 && long_idx == 3) mm_dbg_flag |= MM_DBG_NO_KALLOC; // --no-kalloc
		else if (c == 0 && long_idx == 4) mm_dbg_flag |= MM_DBG_PRINT_QNAME; // --print-qname
		else if (c == 0 && long_idx == 5) opt.flag |= MM_F_NO_SELF; // --no-self
		else if (c == 0 && long_idx == 6) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_SEED; // --print-seed
		else if (c == 0 && long_idx == 7) opt.max_chain_skip = atoi(optarg); // --max-chain-skip
		else if (c == 0 && long_idx == 8) opt.min_ksw_len = atoi(optarg); // --min-dp-len
		else if (c == 0 && long_idx == 9) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_ALN_SEQ; // --print-aln-seq
		else if (c == 'V') {
			return (mm_opt_ind) {opt, NULL};
		} else if (c == 'O') {
			opt.q = opt.q2 = strtol(optarg, &s, 10);
			if (*s == ',') opt.q2 = strtol(s + 1, &s, 10);
		} else if (c == 'E') {
			opt.e = opt.e2 = strtol(optarg, &s, 10);
			if (*s == ',') opt.e2 = strtol(s + 1, &s, 10);
		} else if (c == 'I' || c == 'K') {
			double x;
			char *p;
			x = strtod(optarg, &p);
			if (*p == 'G' || *p == 'g') x *= 1e9;
			else if (*p == 'M' || *p == 'm') x *= 1e6;
			else if (*p == 'K' || *p == 'k') x *= 1e3;
			if (c == 'I') batch_size = (uint64_t)(x + .499);
			else minibatch_size = (uint64_t)(x + .499);
		} else if (c == 'x') {
			if (strcmp(optarg, "ava-ont") == 0) {
				opt.flag |= MM_F_AVA | MM_F_NO_SELF;
				opt.min_chain_score = 100, opt.pri_ratio = 0.0f, opt.max_gap = 10000, opt.max_chain_skip = 25;
				minibatch_size = 500000000;
				k = 15, w = 5;
			} else if (strcmp(optarg, "ava-pb") == 0) {
				opt.flag |= MM_F_AVA | MM_F_NO_SELF;
				opt.min_chain_score = 100, opt.pri_ratio = 0.0f, opt.max_gap = 10000, opt.max_chain_skip = 25;
				minibatch_size = 500000000;
				is_hpc = 1, k = 19, w = 5;
			} else if (strcmp(optarg, "map10k") == 0 || strcmp(optarg, "map-pb") == 0) {
				is_hpc = 1, k = 19;
			} else if (strcmp(optarg, "map-ont") == 0) {
				is_hpc = 0, k = 15;
			} else if (strcmp(optarg, "asm5") == 0) {
				k = 19, w = 19;
				opt.a = 1, opt.b = 19, opt.q = 39, opt.q2 = 81, opt.e = 3, opt.e2 = 1, opt.zdrop = 200;
				opt.min_dp_max = 200;
			} else if (strcmp(optarg, "asm10") == 0) {
				k = 19, w = 19;
				opt.a = 1, opt.b = 9, opt.q = 16, opt.q2 = 41, opt.e = 2, opt.e2 = 1, opt.zdrop = 200;
				opt.min_dp_max = 200;
			} else {
				return (mm_opt_ind) {opt, NULL};
			}
		}
	}
	if (w < 0) w = (int)(.6666667 * k + .499);

	is_idx = mm_idx_is_idx(argv[optind]);
	if (is_idx < 0) {
		return (mm_opt_ind) {opt, NULL};
	}

	mm_idx_t* mi = NULL;
	if (is_idx){
		FILE* fpr = fopen(argv[optind], "rb");
		mi = mm_idx_load(fpr);
		fclose(fpr);
	} else {
		mm_bseq_file_t* fp = mm_bseq_open(argv[optind]);
		if (!mm_bseq_eof(fp)) {
			mi = mm_idx_gen(fp, w, k, bucket_bits, is_hpc, minibatch_size, n_threads, batch_size, keep_name);
		}
		mm_bseq_close(fp);
	}

	mm_mapopt_update(&opt, mi);
	return (mm_opt_ind) {opt, mi};
}
