#include <Python.h>
#include <stdio.h>
#include <zlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include "bseq.h"
#include "minimap.h"

static PyMethodDef module_functions[] = {
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
  #define MOD_ERROR_VAL NULL
  #define MOD_SUCCESS_VAL(val) val
  #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          static struct PyModuleDef moduledef = { \
            PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
          ob = PyModule_Create(&moduledef);
#else
  #define MOD_ERROR_VAL
  #define MOD_SUCCESS_VAL(val)
  #define MOD_INIT(name) void init##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          ob = Py_InitModule3(name, methods, doc);
#endif

MOD_INIT(minimappy)
{
    PyObject *m;

    MOD_DEF(m, "minimaplib", "High-level binding to minimap2",
            module_functions)

    if (m == NULL)
        return MOD_ERROR_VAL;
    return MOD_SUCCESS_VAL(m);

}

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#  ifdef MODULE_API_EXPORTS
#    define MODULE_API __declspec(dllexport)
#    define restrict __restrict
#  else
#    define MODULE_API __declspec(dllimport)
#  endif
#else
#  define MODULE_API
#endif

MODULE_API int module_init();

#ifdef __cplusplus
}
#endif


typedef struct {size_t n; mm_reg1_t* reg;} mm_reg1_v;

MODULE_API mm_idx_t* get_index(const char* fname){
	int k = 15, w = -1, bucket_bits = MM_IDX_DEF_B, n_threads = 3, keep_name = 1, is_hpc = 0;
	int minibatch_size = 200000000;
	uint64_t batch_size = 4000000000ULL;

	int is_idx = mm_idx_is_idx(fname);
	if (is_idx < 0) {
		return NULL;
	}

	mm_idx_t* mi = NULL;
	if (is_idx){
		FILE* fpr = fopen(fname, "rb");
		mi = mm_idx_load(fpr);
		fclose(fpr);
	} else {
		mm_bseq_file_t* fp = mm_bseq_open(fname);
		if (!mm_bseq_eof(fp)) {
			mi = mm_idx_gen(fp, w, k, bucket_bits, is_hpc, minibatch_size, n_threads, batch_size, keep_name);
		}
		mm_bseq_close(fp);
	}
	return mi;
}

MODULE_API void destroy_index(mm_idx_t* index){
	mm_idx_destroy(index);
}

MODULE_API mm_reg1_v align(mm_idx_t* mi, char* seq, char* seq_name) {
	mm_verbose = 2;

	mm_mapopt_t opt;
	mm_mapopt_init(&opt);
	mm_mapopt_update(&opt, mi);
	opt.flag |= MM_F_CIGAR;

	mm_tbuf_t* tbuf = mm_tbuf_init();
	int n_reg;
	const mm_reg1_t* reg = mm_map(mi, strlen(seq), seq, &n_reg, tbuf, &opt, 0);
	mm_tbuf_destroy(tbuf);

	return (mm_reg1_v) {n_reg, reg};
}

