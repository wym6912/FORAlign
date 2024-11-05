#ifndef __HIRSCHBERG_H__
#define __HIRSCHBERG_H__
#ifdef __cplusplus
#include <cstdio>
#include <cstdlib>
extern "C"
{
#else
#include <stdio.h>
#include <stdlib.h>
#endif

#if (defined(_WIN32) || defined(_WIN64))
#ifdef RUSSIAN_LIB
#define DLL __declspec(dllexport)
#else
#define DLL __declspec(dllimport)
#endif // #ifdef RUSSIAN_LIB
#else
#define DLL
#endif

#ifndef DP_LONG_LEN
#define DP_LONG_LEN 1
#endif

#if DP_LONG_LEN == 2
typedef int128_t dp_t;
#elif DP_LONG_LEN == 1
typedef long long dp_t;
#else // DP_LONG_LEN == 0
typedef int dp_t;
#endif
#ifndef sequence_t
typedef enum { DNA, Protein } sequence_t;
#endif

DLL void read_2_seqs(char *c_input, char **seq1, char **name1, char **seq2, char **name2, size_t *len1, size_t *len2);
DLL void print_seqs_to_file(char *seq1, char *comment1, char *seq2, char *comment2, char* c_output);
DLL void hirschberg_API(char *seq1, char *seq2, int len1, int len2, dp_t O, dp_t E, dp_t M, dp_t X, int table_size, sequence_t type_, int threads, char **out_seq1, char **out_seq2, short simple_dp, short multidp, int dp_threads);
DLL void hirschberg_cigar(char* seq1, char* seq2, int len1, int len2, dp_t O, dp_t E, dp_t M, dp_t X, int table_size, sequence_t type_, int threads, char* out_cigar, int *cigar_end, short simple_dp, short multidp, int dp_threads);
DLL void hirschberg_single_cigar(char *seq1, char *seq2, int len1, int len2, dp_t O, dp_t E, dp_t M, dp_t X, sequence_t type_, char* out_cigar, int* cigar_end);
#ifdef __cplusplus
}
#endif
#endif