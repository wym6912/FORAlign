#include "hirschberg_single.hpp"

#define MATCH_SCORE(a, b, M, X) ((a) == (b) ? (M) : (X))
#define HIRSCHBERG_GAP(N, O, E) ((N) <= 0 ? 0 : (O) + (E) * (N))


inline void del_single(int s, std::deque <int> &result)
{
    if (result.size() && result.back() < 0)
        result.back() -= s;
    else result.emplace_back(-s);
}

inline void add_single(int s, std::deque <int> &result)
{
    if (result.size() && result.back() < 0)
    {
        auto tmp = result.back();
        result.back() = s;
        result.emplace_back(std::move(tmp));
    }
    else result.emplace_back(s);
}

inline void match_single(int s, std::deque <int> &result)
{
    while (s--) result.emplace_back(0);
}


void hirschberg_single(const char* seq1, const char* seq2, int b1, int l1, int b2, int l2, dp_t tb, dp_t te, dp_t O, dp_t E, dp_t M, dp_t X, std::deque <int> &result)
{
    /* if sequence B is empty.... */
    if(l2 <= 0)
    {
        /* if sequence A is not empty.... */
        if(l1 > 0)
        {
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "del l1: b1 = %zu, l1 = %zu (final)\n", b1, l1);
#endif
            del_single(l1, result);
        }
        return;
    }
    /* if sequence A is empty.... */
    else if(l1 <= 1)
    {
        if(l1 <= 0)
        {
            /* insert residues from B[1] to B[N] */
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "add l2: b2 = %zu, l2 = %zu (final)\n", b2, l2);
#endif
            add_single(l2, result);
            return;
        }
        /* if sequence A has just one residue.... */
        int midj = 0;
        dp_t min_score, cur;
        min_score = (tb + E) + HIRSCHBERG_GAP(l2, te, E);
        cur       = (te + E) + HIRSCHBERG_GAP(l2, tb, E);
        if(min_score > cur) min_score = cur;
        for(int j = 1; j <= l2; ++ j)
        {
            cur = HIRSCHBERG_GAP(j - 1, tb, E) + MATCH_SCORE(seq1[b1 + l1 - 1], seq2[b2 + j - 1], M, X) + HIRSCHBERG_GAP(l2 - j, te, E);
            if(cur < min_score)
            {
                min_score = cur;
                midj = j;
            }
        }

        if(! midj)
        {
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "del 1: b1 = %zu, l1 = %zu\n", b1, l1);
#endif
            del_single(1, result);
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "add l2: b2 = %zu, l2 = %zu\n", b2, l2);
#endif
            add_single(l2, result);
        }
        else
        {
            if(midj > 1)
            {
#if HIRSCHBERG_DEBUG
                fprintf(stderr, "add l2: b2 = %zu, midj - 1 = %zu (before)\n", b2, midj - 1);
#endif
                add_single(midj - 1, result);
            }
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "match 1: midi = %zu, midj = %zu\n", l1 + b1, midj + b2);
#endif
            match_single(1, result);
            if(midj < l2)
            {
#if HIRSCHBERG_DEBUG
                fprintf(stderr, "add l2: b2 = %zu, l2 - midj = %zu (after)\n", b2, l2 - midj);
#endif
                add_single(l2 - midj, result);
            }
        }
        return;
    }
    else
    {
        size_t imid = l1 >> 1, jmid;
        dp_t *CC = AllocateDPVec(l2 + 1), *DD = AllocateDPVec(l2 + 1), *RR = AllocateDPVec(l2 + 1), *SS = AllocateDPVec(l2 + 1);
        linear_dp_hirschberg    (seq1 + b1, seq2 + b2, CC, DD, imid,          l2, M, X, O, E, tb);
        rev_linear_dp_hirschberg(seq1 + b1, seq2 + b2, RR, SS, l1 - imid, l1, l2, M, X, O, E, te);
        bool type; dp_t min_val, cur_val;
        /* find midj, such that CC[j] + RR[j] or DD[j] + SS[j] + gap is the max */
        min_val = CC[0] + RR[l2]; jmid = 0; type = false;
        // type 1
        for(size_t j = 0; j <= l2; ++ j)
        {
            cur_val = CC[j] + RR[l2 - j];
            if(cur_val <= min_val)
            {
                if(cur_val < min_val || (CC[j] != DD[j] && RR[l2 - j] == SS[l2 - j]))
                {
                    min_val = cur_val;
                    jmid = j;
                }
            }
        }
        // type 2
        for(int j = l2; j >= 0; -- j)
        {
            cur_val = DD[j] + SS[l2 - j] - O;
            if(cur_val < min_val)
            {
                type = true;
                jmid = j;
                min_val = cur_val;
            }
        }
#if HIRSCHBERG_DEBUG
        fprintf(stderr, "CC    DD    RR    SS\n");
        for(size_t i = 0; i <= l2; ++ i)
        {
            fprintf(stderr, "%-5d %-5d %-5d %d\n", CC[i], DD[i], RR[i], SS[i]);
        }
#endif
        FreeDPVec(CC); FreeDPVec(DD); FreeDPVec(RR); FreeDPVec(SS);
        /* Conquer recursively around midpoint */
        if(! type)
        {
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "This function will choose type 1 with (%zu, %zu) => (%zu, %zu) and (%zu, %zu) => (%zu, %zu); max = %d, (i*, j*) = (%zu, %zu)\n",
                            b1, imid, b2, jmid, b1 + imid, l1 - imid, b2 + jmid, l2 - jmid, max_val, imid, jmid);
#endif
            hirschberg_single(seq1, seq2, b1,        imid,      b2,        jmid,      tb, O,  O, E, M, X, result);
            hirschberg_single(seq1, seq2, b1 + imid, l1 - imid, b2 + jmid, l2 - jmid, O,  te, O, E, M, X, result);
        }
        else
        {
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "This function will choose type 2 with (%zu, %zu) => (%zu, %zu), with 2 deletions, and (%zu, %zu) => (%zu, %zu); max = %d, (i*, j*) = (%zu, %zu)\n",
                            b1, imid - 1, b2, jmid, b1 + imid + 1, l1 - imid - 1, b2 + jmid, l2 - jmid, max_val, imid, jmid);
#endif
            hirschberg_single(seq1, seq2, b1,            imid - 1,      b2,        jmid,      tb, 0,  O, E, M, X, result);
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "del 2: midi = %zu\n", imid + b1);
#endif
            del_single(2, result);
            hirschberg_single(seq1, seq2, b1 + imid + 1, l1 - imid - 1, b2 + jmid, l2 - jmid, 0,  te, O, E, M, X, result);
        }
    }
}

void write_final_result(size_t len1, size_t len2, char* seq1, char* seq2, char* cigar, int* final_cigar, std::deque <int> &result)
{
    while (!result.empty())
    {
        auto top = result.front();
        if (!top)
        {
            if (*seq1++ == *seq2++) *cigar++ = 'M';
            else *cigar++ = 'X';
            (*final_cigar)++;
        }
        else if (top > 0)
        {
            while (top--)
            {
                *cigar++ = 'I';
                seq2++;
                (*final_cigar)++;
            }
        }
        else
        {
            auto k = -top;
            while(k--)
            {
                *cigar++ = 'D';
                seq1++;
                (*final_cigar)++;
            }
        }
        result.pop_front();
    }
}


DLL void hirschberg_single_cigar(char *seq1, char *seq2, int len1, int len2, dp_t O, dp_t E, dp_t M, dp_t X, sequence_t type_, char* out_cigar, int* cigar_end)
{
    int encode_size;
    char *russian_table = nullptr, *comp_seq1, *comp_seq2;
    std::deque <int> result;
    russian_table = russian_table_init(type_, russian_table, encode_size);
    str_compress(seq1, comp_seq1, len1, russian_table);
    str_compress(seq2, comp_seq2, len2, russian_table);
    hirschberg_single(comp_seq1, comp_seq2, 0, len1, 0, len2, O, O, O, E, M, X, result);
    FreeCharVec(comp_seq1); FreeCharVec(comp_seq2); // save memory
    write_final_result(len1, len2, seq1, seq2, out_cigar, cigar_end, result);
    FreeCharVec(russian_table);
}
