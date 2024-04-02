#include "mtxutl.hpp"

char *AllocateCharVec(size_t len)
{
    char *vec;
	vec = (char *)calloc(len, sizeof(char));
	if(vec == nullptr)
	{	
		fprintf(stderr, "Error: Allocation error (char). Program will exit.\n");
		exit(1);
	}
	return vec;
}

void FreeCharVec(char *v)
{
    free(v);
}

int *AllocateIntVec(size_t len)
{
    int *vec;
	vec = (int *)calloc(len, sizeof(int));
	if(vec == nullptr)
	{	
		fprintf(stderr, "Error: Allocation error (vector). Program will exit.\n");
		exit(1);
	}
	return vec;
}

void FreeIntVec(int *v)
{
    free(v);
}

int **AllocateIntMtx(size_t a, size_t b)
{
    size_t i;
	int **mtx;

	mtx = (int **)calloc(a + 1, sizeof(int*));
	if(mtx == nullptr)
	{
		fprintf( stderr, "Error: Allocation error (matrix). Program will exit.\n");
		exit(1);
	}
	if(b)
	{
		for(i = 0; i < a; ++ i) mtx[i] = AllocateIntVec(b);
	}
	else
	{
		for(i = 0; i < a; ++ i) mtx[i] = nullptr;
	}
	mtx[a] = nullptr;

	return mtx;
}

void FreeIntMtx(int **m)
{
    size_t i;

	for( i = 0; m[i]; ++ i) 
		if(m[i]) { free(m[i]); m[i] = nullptr; }
	free(m);
}

ll *AllocateLongLongVec(size_t len)
{
    ll *vec;
	vec = (ll *)calloc(len, sizeof(ll));
	if(vec == nullptr)
	{	
		fprintf(stderr, "Error: Allocation error (long long vector). Program will exit.\n");
		exit(1);
	}
	return vec;
}

void FreeLongLongVec(ll *v)
{
    free(v);
}

ll **AllocateLongLongMtx(size_t a, size_t b)
{
    size_t i;
	ll **mtx;

	mtx = (ll **)calloc(a + 1, sizeof(ll*));
	if(mtx == nullptr)
	{
		fprintf( stderr, "Error: Allocation error (long long matrix). Program will exit.\n");
		exit(1);
	}
	if(b)
	{
		for(i = 0; i < a; ++ i) mtx[i] = AllocateLongLongVec(b);
	}
	else
	{
		for(i = 0; i < a; ++ i) mtx[i] = nullptr;
	}
	mtx[a] = nullptr;

	return mtx;
}

void FreeLongLongMtx(ll **m)
{
    size_t i;

	for( i = 0; m[i]; ++ i) 
		if(m[i]) { free(m[i]); m[i] = nullptr; }
	free(m);
}

int128_t *AllocateInt128Vec(size_t a)
{
    int128_t *vec;
	vec = (int128_t *)calloc(a, sizeof(int128_t));
	if(vec == nullptr)
	{	
		fprintf(stderr, "Error: Allocation error (vector). Program will exit.\n");
		exit(1);
	}
	return vec;

}

void FreeInt128Vec(int128_t* v)
{
	free(v);
}

int128_t **AllocateInt128Mtx(size_t a, size_t b)
{
    size_t i;
	int128_t **mtx;

	mtx = (int128_t **)calloc(a + 1, sizeof(int128_t*));
	if(mtx == nullptr)
	{
		fprintf( stderr, "Error: Allocation error (matrix). Program will exit.\n");
		exit(1);
	}
	if(b)
	{
		for(i = 0; i < a; ++ i) mtx[i] = AllocateInt128Vec(b);
	}
	else
	{
		for(i = 0; i < a; ++ i) mtx[i] = nullptr;
	}
	mtx[a] = nullptr;

	return mtx;
}

void FreeInt128Mtx(int128_t **m)
{
    size_t i;

	for( i = 0; m[i]; ++ i) 
		if(m[i]) { free(m[i]); m[i] = nullptr; }
	free(m);
}
