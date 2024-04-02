#ifndef __MTXUTL__
#define __MTXUTL__

#include <cstdlib>
#include <cstdio>
#include "../include/boost/include/boost/multiprecision/cpp_int.hpp"

typedef long long ll;
using namespace boost::multiprecision;

char *AllocateCharVec(size_t len);
void FreeCharVec(char *v);

int *AllocateIntVec(size_t len);
void FreeIntVec(int *v);

int **AllocateIntMtx(size_t a, size_t b);
void FreeIntMtx(int **m);

ll *AllocateLongLongVec(size_t len);
void FreeLongLongVec(ll *v);

ll **AllocateLongLongMtx(size_t a, size_t b);
void FreeLongLongMtx(ll **m);

int128_t *AllocateInt128Vec(size_t a);
void FreeInt128Vec(int128_t* v);

int128_t **AllocateInt128Mtx(size_t a, size_t b);
void FreeInt128Mtx(int128_t **m);

#endif