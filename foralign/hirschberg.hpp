#ifndef __HIRSCHBERG__
#define __HIRSCHBERG__

#include "russians.hpp"
#include "mtdp.hpp"
#include "mtset.hpp"
#include "../include/dpthreadpool/BS_thread_pool_light.hpp"
#include <queue>

using namespace BS;

void hirschberg_multi_init(ConcurrentSet_hirschberg_status_t *&S, F_affine_t *&F, int len, thread_pool_light *&pool1, thread_pool_light *&pool2, int dp_threads, std::mutex *&mtx1, std::mutex *&mtx2, bool multidp);
void hirschberg_multi_start(char *seq1, char *seq2, int len1, int len2, dp_t O, dp_t E, dp_t M, dp_t X, int russian_table_size, ConcurrentSet_hirschberg_status_t *&status_set, F_affine_t *F_, bool multidp, int dp_threads, 
                            thread_pool_light *pool1, thread_pool_light *pool2, std::mutex *mtx1, std::mutex *mtx2, int encode_table_size, int threads, bool simple_dp);
void hirschberg_multi_free(ConcurrentSet_hirschberg_status_t *&S, F_affine_t *&F, thread_pool_light *pool1, thread_pool_light *pool2, std::mutex *mtx1, std::mutex *mtx2);
void print_result_multi(char *FILE_OUT, int len1, int len2, char *seq1, char *seq2, char *name1, char *name2, ConcurrentSet_hirschberg_status_t *S);

#endif