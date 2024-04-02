#ifndef __DP_HPP__
#define __DP_HPP__

#include "mtxutl.hpp"
#include <algorithm>
#include <mutex>
#include <cstring>
#include <iostream>
#include <future>
#include <set>
#include <vector>
#include <tuple>
#include <cstdio>

#ifndef DP_LONG_LEN
#define DP_LONG_LEN 0
#endif

#if DP_LONG_LEN == 2
typedef int128_t dp_t;
#elif DP_LONG_LEN == 1
typedef long long dp_t;
#else // DP_LONG_LEN == 0
typedef int dp_t;
#endif

typedef std::vector <dp_t> dp_v_t;
typedef std::set <dp_v_t> dp_vs_t;
typedef std::set <dp_t> dp_s_t;

void print_vec(dp_t *v, size_t len);

char *strrev_(char *str);
dp_t* AllocateDPVec(size_t v);
dp_t readDPdata(FILE *F);
void printDPdata(dp_t x, FILE *F);
void FreeDPVec(dp_t *v);
void linear_dp_hirschberg(const char* seq1, const char* seq2, dp_t* CC, dp_t* DD, const int line1, const int len2, dp_t& M, dp_t& X, dp_t& O, dp_t& E, dp_t& tg);
void rev_linear_dp_hirschberg(const char* seq1, const char* seq2, dp_t* CC, dp_t* DD, const int line1, const int len1, const int len2, dp_t& M, dp_t& X, dp_t& O, dp_t& E, dp_t& tg);

void block_affine_dp_russian_final(const int block_size, dp_v_t &C1, dp_v_t &C2, dp_v_t &I1, dp_v_t &D2, const char *seq1, const char *seq2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t *delta);
void block_affine_dp_russian_seq2(const dp_t *C1, dp_t *C2, const dp_t *I1, dp_t *D2, const int block_size, int seq2len, const char *seq1, const char *seq2, dp_t &M, dp_t &X, dp_t &O, dp_t &E);
void rev_block_affine_dp_russian_final(const int block_size, dp_v_t &C1, dp_v_t &C2, dp_v_t &I1, dp_v_t &D2, const char *seq1, const char *seq2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t *delta);
void rev_block_affine_dp_russian_seq2(const dp_t *C1, dp_t *C2, const dp_t *I1, dp_t *D2, const int block_size, const int seq2len, const char *seq1, const char *seq2, dp_t &M, dp_t &X, dp_t &O, dp_t &E);

#endif
