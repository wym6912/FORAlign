#ifndef __PARALLEL_SWG__
#define __PARALLEL_SWG__

#ifdef __cplusplus
extern "C"
{
#endif

#if (defined(_WIN32) || defined(_WIN64))
#ifdef SWG_LIB
#define DLL __declspec(dllexport)
#else
#define DLL __declspec(dllimport)
#endif // #ifdef SWG_LIB
#else
#define DLL
#endif

#ifndef AFFINE_MATRIX_H_ // defined in wfa2-system
/*
 * Affine Matrix
 */
typedef struct {
  int M; // Alignment matching/mismatching
  int I; // Alignment ends with a gap in the reference (insertion)
  int D; // Alignment ends with a gap in the read (deletion)
} affine_cell_t;
typedef struct {
  // Affine Matrix
  affine_cell_t** columns;
  int num_rows;
  int num_columns;
} affine_matrix_t;
/*
 * Affine penalties
 */
typedef struct {
  int match;              // (Penalty representation; usually M <= 0)
  int mismatch;           // (Penalty representation; usually X > 0)
  int gap_opening;        // (Penalty representation; usually O > 0)
  int gap_extension;      // (Penalty representation; usually E > 0)
} affine_penalties_t;

/*
 * Affine matrix-type (for traceback)
 */
typedef enum {
  affine_matrix_M,
  affine_matrix_I,
  affine_matrix_D,
} affine_matrix_type;
#endif


DLL void multithread_swg_init(int threads);
DLL void multithread_swg_free();
DLL void multithread_swg_compute_cv(int threads, affine_matrix_t* const affine_table, affine_penalties_t* const penalties, const char* const pattern, const int pattern_length, const char* const text, const int text_length);
DLL void multithread_swg_compute_barrier(int threads, affine_matrix_t* const affine_table, affine_penalties_t* const penalties, const char* const pattern, const int pattern_length, const char* const text, const int text_length);
DLL void multithread_swg_traceback(affine_matrix_t* const affine_matrix, affine_penalties_t* const penalties, const int pattern_length, const int text_length, char *seq1, char *seq2, char **oseq1, char **oseq2);
#ifdef __cplusplus
}
#endif
#endif