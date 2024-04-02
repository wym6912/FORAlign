#include "dp.hpp"
#define DP_DEBUG 0

#define MATCH_SCORE(a, b, M, X) ((a) == (b) ? (M) : (X))

inline dp_t** AllocateDPMtx(size_t a, size_t b)
{
#if DP_LONG_LEN == 2
    return AllocateInt128Mtx(a, b);
#elif DP_LONG_LEN == 1
    return AllocateLongLongMtx(a, b);
#else
    return AllocateIntMtx(a, b);
#endif
}

inline void FreeDPMtx(dp_t **m)
{
#if DP_LONG_LEN == 2
    FreeInt128Mtx(m);
#elif DP_LONG_LEN == 1
    FreeLongLongMtx(m);
#else
    FreeIntMtx(m);
#endif
}

dp_t* AllocateDPVec(size_t v)
{
#if DP_LONG_LEN == 2
    return AllocateInt128Vec(v);
#elif DP_LONG_LEN == 1
    return AllocateLongLongVec(v);
#else
    return AllocateIntVec(v);
#endif
}

dp_t readDPdata(FILE *F)
{
    dp_t x;
#if DP_LONG_LEN == 2
    x = 0;
    dp_t f = 1;
    char c = fgetc(F);
    while(! isdigit(c)){ if(c == '-') f = -1;  c = fgetc(F); }
    while(isdigit(c))  { x = x * 10 + c - '0'; c = fgetc(F); }
#elif DP_LONG_LEN == 1
    fscanf(F, "%lld", &x);
#else
    fscanf(F, "%d", &x);
#endif
    return x;
}

void printDPdata(dp_t x, FILE *F)
{
#if DP_LONG_LEN == 2
    fprintf(F, "%s", x.str().c_str());
#elif DP_LONG_LEN == 1
    fprintf(F, "%lld", x);
#else
    fprintf(F, "%d", x);
#endif
}

void FreeDPVec(dp_t *v)
{
#if DP_LONG_LEN == 2
    FreeInt128Vec(v);
#elif DP_LONG_LEN == 1
    FreeLongLongVec(v);
#else
    FreeIntVec(v);
#endif

}

void print_matrix(dp_t **mat, size_t len1, size_t len2)
{
    for(size_t i = 0; i <= len1; ++ i, std::cerr << "\n")
        for(size_t j = 0; j <= len2; ++ j, std::cerr << ' ')
            std::cerr << mat[i][j];
    std::cerr << std::endl;
}

void print_vec(dp_t *v, size_t len)
{
    for(size_t i = 0; i <= len; ++ i, std::cerr << ' ') std::cerr << v[i];
    std::cerr << std::endl;
}

char *strrev_(char *str)
{
    // code from https://stackoverflow.com/a/8534275 - it works!
    char *p1, *p2;
    if (! str || ! *str)
        return str;
    for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2)
    {
        *p1 ^= *p2;
        *p2 ^= *p1;
        *p1 ^= *p2;
    }
    return str;
}

void block_affine_dp_russian_seq2(const dp_t *C1, dp_t *C2, const dp_t *I1, dp_t *D2, const int block_size, const int seq2len, const char *seq1, const char *seq2, dp_t &M, dp_t &X, dp_t &O, dp_t &E)
{
    dp_t *C = AllocateDPVec(seq2len + 1), *D = AllocateDPVec(seq2len + 1), e, c, s;
    // round 0
    C[0] = 0;
    for(int j = 1; j <= seq2len; ++ j)
    {
        C[j] = C[j - 1] + C2[j - 1];
        D[j] = C[j] + D2[j - 1] + (O + E);
    }
    // other_rounds
    for(int i = 1; i <= block_size; ++ i)
    {
        // copy C1 and I1
        s = C[0];                       // s = C[i - 1][0]
        C[0] = c = s + C1[i - 1];       // c = C[i][0]
        e = C[0] + I1[i - 1] + (O + E); // e = I[i][0]
        for(int j = 1; j <= seq2len; ++ j)
        {
            e = std::min(e, c + O) + E;          // e = I[i][j]
            D[j] = std::min(D[j], C[j] + O) + E; // DD[j] = D[i][j]
            c = std::min(std::min(D[j], e), s + MATCH_SCORE(seq1[i - 1], seq2[j - 1], M, X));
            s = C[j];                            // s = C[i - 1][j]
            C[j] = c;                            // c = C[i][j]
        }
    }
    for(int j = 1; j <= seq2len; ++ j)
    {
        C2[j - 1] = C[j] - C[j - 1];
        D2[j - 1] = D[j] - C[j] - (O + E);
    }
    FreeDPVec(C); FreeDPVec(D);
}

void block_affine_dp_russian_final(const int block_size, dp_v_t &C1, dp_v_t &C2, dp_v_t &I1, dp_v_t &D2, const char *seq1, const char *seq2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t *delta)
{
    dp_t *C = AllocateDPVec(block_size + 1), *D = AllocateDPVec(block_size + 1), e, c, s;
    // round 0
    C[0] = 0;
    for(int j = 1; j <= block_size; ++ j)
    {
        C[j] = C[j - 1] + C2[j - 1];
        D[j] = C[j] + D2[j - 1] + (O + E);
    }
    // other_rounds
    for(int i = 1; i <= block_size; ++ i)
    {
        // copy C1 and I1
        s = C[0];                       // s = C[i - 1][0]
        C[0] = c = s + C1[i - 1];       // c = C[i][0]
        e = C[0] + I1[i - 1] + (O + E); // e = I[i][0]
        for(int j = 1; j <= block_size; ++ j)
        {
            e = std::min(e, c + O) + E;          // e = I[i][j]
            D[j] = std::min(D[j], C[j] + O) + E; // DD[j] = D[i][j]
            c = std::min(std::min(D[j], e), s + MATCH_SCORE(seq1[i - 1], seq2[j - 1], M, X));
            s = C[j];                            // s = C[i - 1][j]
            C[j] = c;                            // c = C[i][j]
        }
        C1[i - 1] = c - s;           // s = C[i - 1][block_size]; c = C[i][block_size]
        I1[i - 1] = e - c - (O + E); // e = I[i][block_size]
    }

    *delta = C[block_size];
    for(int j = 1; j <= block_size; ++ j)
    {
        C2[j - 1] = C[j] - C[j - 1];
        D2[j - 1] = D[j] - C[j] - (O + E);
    }
    FreeDPVec(C); FreeDPVec(D);
}

void rev_block_affine_dp_russian_final(const int block_size, dp_v_t &C1, dp_v_t &C2, dp_v_t &I1, dp_v_t &D2, const char *seq1, const char *seq2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t *delta)
{
    dp_t *C = AllocateDPVec(block_size + 1), *D = AllocateDPVec(block_size + 1), e, c, s;
    // round 0
    C[0] = 0;
    for(int j = 1; j <= block_size; ++ j)
    {
        C[j] = C[j - 1] + C2[j - 1];
        D[j] = C[j] + D2[j - 1] + (O + E);
    }
    // other_rounds
    for(int i = 1; i <= block_size; ++ i)
    {
        // copy C1 and I1
        s = C[0];                       // s = C[i - 1][0]
        C[0] = c = s + C1[i - 1];       // c = C[i][0]
        e = C[0] + I1[i - 1] + (O + E); // e = I[i][0]
        for(int j = 1; j <= block_size; ++ j)
        {
            e = std::min(e, c + O) + E;          // e = I[i][j]
            D[j] = std::min(D[j], C[j] + O) + E; // DD[j] = D[i][j]
            c = std::min(std::min(D[j], e), s + MATCH_SCORE(seq1[-(i - 1)], seq2[-(j - 1)], M, X));
            s = C[j];                            // s = C[i - 1][j]
            C[j] = c;                            // c = C[i][j]
        }
        C1[i - 1] = c - s;           // s = C[i - 1][block_size]; c = C[i][block_size]
        I1[i - 1] = e - c - (O + E); // e = I[i][block_size]
    }

    *delta = C[block_size];
    for(int j = 1; j <= block_size; ++ j)
    {
        C2[j - 1] = C[j] - C[j - 1];
        D2[j - 1] = D[j] - C[j] - (O + E);
    }
    FreeDPVec(C); FreeDPVec(D);
}

void rev_block_affine_dp_russian_seq2(const dp_t *C1, dp_t *C2, const dp_t *I1, dp_t *D2, const int block_size, const int seq2len, const char *seq1, const char *seq2, dp_t &M, dp_t &X, dp_t &O, dp_t &E)
{
    dp_t *C = AllocateDPVec(seq2len + 1), *D = AllocateDPVec(seq2len + 1), e, c, s;
    // round 0
    C[0] = 0;
    for(int j = 1; j <= seq2len; ++ j)
    {
        C[j] = C[j - 1] + C2[j - 1];
        D[j] = C[j] + D2[j - 1] + (O + E);
    }
    // other_rounds
    for(int i = 1; i <= block_size; ++ i)
    {
        // copy C1 and I1
        s = C[0];                       // s = C[i - 1][0]
        C[0] = c = s + C1[i - 1];       // c = C[i][0]
        e = C[0] + I1[i - 1] + (O + E); // e = I[i][0]
        for(int j = 1; j <= seq2len; ++ j)
        {
            e = std::min(e, c + O) + E;          // e = I[i][j]
            D[j] = std::min(D[j], C[j] + O) + E; // DD[j] = D[i][j]
            c = std::min(std::min(D[j], e), s + MATCH_SCORE(seq1[-(i - 1)], seq2[-(j - 1)], M, X));
            s = C[j];                            // s = C[i - 1][j]
            C[j] = c;                            // c = C[i][j]
        }
    }

    for(int j = 1; j <= seq2len; ++ j)
    {
        C2[j - 1] = C[j] - C[j - 1];
        D2[j - 1] = D[j] - C[j] - (O + E);
    }
    FreeDPVec(C); FreeDPVec(D);
}

void linear_dp_hirschberg(const char* seq1, const char* seq2, dp_t *CC, dp_t *DD, const int line1, const int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t &tg)
{
    dp_t t, s, e, c;
    CC[0] = 0;
    t = O;
    for(int j = 1; j <= len2; ++ j)
    {
        t = t + E;
        CC[j] = t;
        DD[j] = t + O;
    }
#define DP_DEBUG 0
#if DP_DEBUG
    fprintf(stderr, "CC[0] = "); print_vec(CC, len2);
    fprintf(stderr, "DD[0] = "); print_vec(DD, len2);
#endif
    t = tg; // [*]
    for(int i = 1; i <= line1; ++ i)
    {
        s = CC[0];
        t = t + E;
        c = t;
        CC[0] = c; DD[0] = CC[0];
        e = t + O;
        for(int j = 1; j <= len2; ++ j)
        {
            e = std::min(e, c + O) + E;
            DD[j] = std::min(DD[j], CC[j] + O) + E;
            c = std::min(std::min(DD[j], e), s + MATCH_SCORE(seq1[i - 1], seq2[j - 1], M, X));
            s = CC[j];
            CC[j] = c;
        }
#if DP_DEBUG
        fprintf(stderr, "CC[%zu] = ", i); print_vec(CC, len2);
        fprintf(stderr, "DD[%zu] = ", i); print_vec(DD, len2);
#endif
    }
}

void rev_linear_dp_hirschberg(const char* seq1, const char* seq2, dp_t *CC, dp_t *DD, const int line1, const int len1, const int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t &tg)
{
    dp_t t, s, e, c;
    CC[0] = 0;
    t = O;
    for(int j = 1; j <= len2; ++ j)
    {
        t = t + E;
        CC[j] = t;
        DD[j] = t + O;
    }
#if DP_DEBUG
    fprintf(stderr, "CC[0] = "); print_vec(CC, len2);
    fprintf(stderr, "DD[0] = "); print_vec(DD, len2);
#endif
    t = tg; // [*]
    for(int i = 1; i <= line1; ++ i)
    {
        s = CC[0];
        t = t + E;
        c = t;
        CC[0] = c; DD[0] = CC[0];
        e = t + O;
        for(int j = 1; j <= len2; ++ j)
        {
            e = std::min(e, c + O) + E;
            DD[j] = std::min(DD[j], CC[j] + O) + E;
            c = std::min(std::min(DD[j], e), s + MATCH_SCORE(seq1[len1 - i], seq2[len2 - j], M, X));
            s = CC[j];
            CC[j] = c;
        }
#if DP_DEBUG
        fprintf(stderr, "CC[%zu] = ", i); print_vec(CC, len2);
        fprintf(stderr, "DD[%zu] = ", i); print_vec(DD, len2);
#endif
    }
}
