#include "russians.hpp"
#undef min

#define MATCH_SCORE(a, b, M, X) ((a) == (b) ? (M) : (X))

char* russian_table_init(sequence_t sequence_type, char* table, int& table_size)
{
    table = table ? table : AllocateCharVec(128);
    if(sequence_type == DNA)
    {
        for(char i = 0; i < 127; ++ i)
        {
            switch(i)
            {
            case 'A':
            case 'a':
                table[i] = 0;
                break;
            case 'C':
            case 'c':
                table[i] = 1;
                break;
            case 'G':
            case 'g':
                table[i] = 2;
                break;
            case 'T':
            case 't':
            case 'U':
            case 'u':
                table[i] = 3;
                break;
            default:
                table[i] = 4;
                break;
            }
        }
        table_size = 5;
    }
    else if(sequence_type == Protein)
    {
        for(char i = 0; i < 127; ++ i)
        {
            switch(i)
            {
            case 'A':
            case 'a':
                table[i] = 0;
                break;
            case 'C':
            case 'c':
                table[i] = 1;
                break;
            case 'D':
            case 'd':
                table[i] = 2;
                break;
            case 'E':
            case 'e':
                table[i] = 3;
                break;
            case 'F':
            case 'f':
                table[i] = 4;
                break;
            case 'G':
            case 'g':
                table[i] = 5;
                break;
            case 'H':
            case 'h':
                table[i] = 6;
                break;
            case 'I':
            case 'i':
                table[i] = 7;
                break;
            case 'K':
            case 'k':
                table[i] = 8;
                break;
            case 'L':
            case 'l':
                table[i] = 9;
                break;
            case 'M':
            case 'm':
                table[i] = 10;
                break;
            case 'N':
            case 'n':
                table[i] = 11;
                break;
            case 'P':
            case 'p':
                table[i] = 12;
                break;
            case 'Q':
            case 'q':
                table[i] = 13;
                break;
            case 'R':
            case 'r':
                table[i] = 14;
                break;
            case 'S':
            case 's':
                table[i] = 15;
                break;
            case 'T':
            case 't':
                table[i] = 16;
                break;
            case 'V':
            case 'v':
                table[i] = 17;
                break;
            case 'W':
            case 'w':
                table[i] = 18;
                break;
            case 'Y':
            case 'y':
                table[i] = 19;
                break;
            default:
                table[i] = 20;
                break;
            }
        }
        table_size = 21;
    }
    return table;
}

void compress_without_table(const char *seq1, const char *seq2, int len, int table_size, std::string &ans)
{
    std::vector <short> v(table_size, 0);
    int count = 1;
    size_t block_len = len;
    auto it = ans.begin();
    while(block_len --)
    {
        if(! v[(*seq1) - 1]) v[(*seq1) - 1] = count ++;
        *it = (v[(*seq1) - 1] - 1);
        ++ it;
        ++ seq1;
    }
    while(len --)
    {
        if(! v[(*seq2) - 1]) v[(*seq2) - 1] = count ++;
        *it = (v[(*seq2) - 1] - 1);
        ++ it;
        ++ seq2;
    }
}

void rev_compress_without_table(const char *seq1, const char *seq2, int len, int table_size, std::string &ans)
{
    std::vector <short> v(table_size, 0);
    int count = 1;
    size_t block_len = len;
    auto it = ans.begin();
    while(block_len --)
    {
        if(! v[(*seq1) - 1]) v[(*seq1) - 1] = count ++;
        *it = (v[(*seq1) - 1] - 1);
        ++ it;
        -- seq1;
    }
    while(len --)
    {
        if(! v[(*seq2) - 1]) v[(*seq2) - 1] = count ++;
        *it = (v[(*seq2) - 1] - 1);
        ++ it;
        -- seq2;
    }
}

size_t compress_with_table(char *seq1, char *seq2, int len, int table_size, char* table)
{
    std::vector <short> v(table_size, 0);
    short count = 1;
    size_t ans = 0, block_len = len;
    while(block_len --)
    {
        if(! v[*(table + *seq1)]) v[*(table + *seq1)] = count ++;
        ans = ans * table_size + v[*(table + *seq1)] - 1;
        ++ seq1;
    }
    while(len --)
    {
        if(! v[*(table + *seq2)]) v[*(table + *seq2)] = count ++;
        ans = ans * table_size + v[*(table + *seq2)] - 1;
        ++ seq2;
    }
    return ans;
}

void str_compress(char *seq, char *&encoded_seq, int len, char *table)
{
    encoded_seq = AllocateCharVec(len + 1);
    char *tmp = encoded_seq;
    while(len --) *tmp ++ = table[*seq ++] + 1;
    *tmp = 0;
}

void linear_block_affine_main_dp(char *real_seq1, int line1, char *real_seq2, int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t tb, int &block_size, F_affine_t *F, dp_t* final_C, dp_t* final_D, int &table_size)
{
    if(line1 < block_size || len2 < block_size)
    {
        // do simple dp, no need to do russians
        linear_dp_hirschberg((const char*)real_seq1, (const char*)real_seq2, final_C, final_D, (const int)line1, (const int)len2, M, X, O, E, tb);
        return;
    }
    int seq1_remain_length = 0, seq2_remain_length = 0;
    // seq1 will roll in last rounds
    if(line1 % block_size != 0) seq1_remain_length = int(line1 % block_size);
    // seq2 will roll as same
    if(len2 % block_size != 0)  seq2_remain_length = int(len2 % block_size);
    int block1 = line1 / block_size, block2 = len2 / block_size, block2_remain_start = block2 * block_size;

    dp_v_t b1(block_size), I1(block_size);
    dp_b_t b2(block2, block_size), D2(block2, block_size);
    dp_v_t bleft(block2 + 1);
    bleft[0] = 0; b1[0] = tb + E; I1[0] = O - (O + E);
    for(int i = 1; i < block_size; ++ i) { b1[i] = E; I1[i] = O - (O + E); }
    for(int block2_ = 0; block2_ < block2; ++ block2_)
    {
        auto &C_block = b2[block2_], &D_block = D2[block2_];
        for(int i = 0; i < block_size; ++ i)
        {
            C_block[i] = E; D_block[i] = O - (O + E);
        }
        bleft[block2_ + 1] = O + E * (block2_ + 1) * block_size;
    }
    b2[0][0] = O + E;
    if(seq2_remain_length)
    {
        for(int i = block2_remain_start + 1; i <= len2; ++ i)
        {
            final_C[i] = E; final_D[i] = O - (O + E);
        }
    }
// #define RUSSIAN_DP_DEBUG 0
#ifdef RUSSIAN_DP_DEBUG
    fprintf(stderr, "b1 = "); print_vec(b1, block_size);
    fprintf(stderr, "b2 = "); print_vec(final_C, block2_remain_start);
    print_vec(bleft, block2 + 1);
#endif
    affine_input_t now_in;
    affine_output_t now_out;
    std::string this_block_s((2 * block_size), 0);
    dp_t dp_before, dp_next;
    for(int block1_ = 0; block1_ < block1; ++ block1_)
    {
        dp_before = bleft[0];
        dp_next = bleft[0];
        for(int block2_ = 0; block2_ < block2; ++ block2_)
        {
            size_t b1_first = block_size * block1_, b2_first = block_size * block2_;
            char *seq1_ = real_seq1 + b1_first, *seq2_ = real_seq2 + b2_first;
            dp_t dp_tmp;
            compress_without_table(seq1_, seq2_, block_size, table_size, this_block_s);
            now_in = std::tie(this_block_s, b1, b2[block2_], I1, D2[block2_]);
#ifdef RUSSIAN_DP_DEBUG
            fprintf(stderr, "block (%zu, %zu):\n", block1_, block2_);
            fprintf(stderr, "b1 = "); print_vec(b1, block_size - 1);
            fprintf(stderr, "b2 = "); print_vec(b2[block2_], block_size - 1);
            fprintf(stderr, "I1 = "); print_vec(I1, block_size - 1);
            fprintf(stderr, "D2 = "); print_vec(D2[block2_], block_size - 1);
            fprintf(stderr, "bleft = %lld\n", bleft[block2_]);
            fprintf(stderr, "val = %s\n", this_name.c_str());
            for(int q = 0; q < block_size; ++ q) fputc(seq1_[q] + '0', stderr); fputc('\n', stderr);
            for(int q = 0; q < block_size; ++ q) fputc(seq2_[q] + '0', stderr); fputc('\n', stderr);
#endif
            if(F -> find(now_in, now_out))
            {
                // found in table F
                std::tie(b1, b2[block2_], I1, D2[block2_], dp_tmp) = now_out;
            }
            else
            {
                block_affine_dp_russian_final(block_size, b1, b2[block2_], I1, D2[block2_], seq1_, seq2_, M, X, O, E, &dp_tmp);
                now_out = std::tie(b1, b2[block2_], I1, D2[block2_], dp_tmp);
                F -> insert(now_in, now_out);
            }
            dp_next = bleft[block2_] + dp_tmp;
            bleft[block2_] = dp_before;
            dp_before = dp_next;
#ifdef RUSSIAN_DP_DEBUG
            fprintf(stderr, "=>\nbleft = %lld\n", dp_next);
            fprintf(stderr, "b1 = "); print_vec(b1, block_size - 1);
            fprintf(stderr, "b2 = "); print_vec(b2[block2_], block_size - 1);
            fprintf(stderr, "I1 = "); print_vec(I1, block_size - 1);
            fprintf(stderr, "D2 = "); print_vec(D2[block2_], block_size - 1);
            fputs("---------\n", stderr);
#endif
        }
        bleft[block2] = dp_next;
#ifdef RUSSIAN_DP_DEBUG
        print_vec(bleft, block2);
#endif
        if(seq2_remain_length)
        {
            block_affine_dp_russian_seq2(b1.data(), final_C + block2_remain_start + 1, I1.data(), final_D + block2_remain_start + 1, block_size, seq2_remain_length, real_seq1 + block_size * block1_, real_seq2 + block2_remain_start, M, X, O, E);
        }
        for(int i = 0; i < block_size; ++ i) { b1[i] = E; I1[i] = O - (O + E); } // next round
        bleft[0] = tb + E * block_size * (block1_ + 1);
    }
    final_C[0] = bleft[0];
    // save to C and D
    int now_block = 0, b2table = 0;
    for(int i = 1; i <= block2_remain_start; ++ i)
    {
        final_C[i] = final_C[i - 1] + b2[b2table][i - now_block - 1];
        final_D[i] = final_C[i] + (O + E) + D2[b2table][i - now_block - 1];
        if(i % block_size == 0) now_block += block_size, ++ b2table;
    }
    for(int i = block2_remain_start + 1; i <= len2; ++ i)
    {
        final_C[i] += final_C[i - 1];
        final_D[i] += final_C[i] + (O + E);
    }
    decltype(bleft)().swap(bleft);
    decltype(b1)().swap(b1);
    decltype(I1)().swap(I1);
    decltype(b2.dp_block)().swap(b2.dp_block);
    decltype(D2.dp_block)().swap(D2.dp_block);
#ifdef RUSSIAN_DP_DEBUG
    fprintf(stderr, "C = "); for(size_t i = 0; i <= len2; ++ i) fprintf(stderr, "%lld ", final_C[i]); fputc('\n', stderr);
    fprintf(stderr, "D = "); for(size_t i = 0; i <= len2; ++ i) fprintf(stderr, "%lld ", final_D[i]); fputc('\n', stderr);
#endif
    // process remain in seq1
    // dp_before defined before
    int now_seq_1 = line1 - seq1_remain_length;
    dp_t final_I1;
    while(seq1_remain_length --)
    {
        dp_t c, s;
        s = final_C[0];
        final_C[0] = c = final_C[0] + E; // must have block, so no need to add tb
        final_I1 = final_C[0] + O;
        for(int i = 1; i <= len2; ++ i)
        {
            final_I1 = std::min(final_I1, c + O) + E;
            final_D[i] = std::min(final_D[i], final_C[i] + O) + E;
            c = std::min(std::min(final_D[i], final_I1), s + MATCH_SCORE(real_seq1[now_seq_1], real_seq2[i - 1], M, X));
            s = final_C[i];
            final_C[i] = c;
        }
        ++ now_seq_1;
#ifdef RUSSIAN_DP_DEBUG
        fprintf(stderr, "C = "); for(size_t i = 0; i <= len2; ++ i) fprintf(stderr, "%lld ", final_C[i]); fputc('\n', stderr);
        fprintf(stderr, "D = "); for(size_t i = 0; i <= len2; ++ i) fprintf(stderr, "%lld ", final_D[i]); fputc('\n', stderr);
#endif
    }
    final_D[0] = final_C[0];
}

void rev_linear_block_affine_main_dp(char *real_seq1, int line1, int len1, char *real_seq2, int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t &tb, int &block_size, F_affine_t *F, dp_t* final_C, dp_t* final_D, int &table_size)
{
    if(line1 < block_size || len2 < block_size)
    {
        // do simple dp, no need to do russians
        rev_linear_dp_hirschberg((const char*)real_seq1, (const char*)real_seq2, final_C, final_D, (const int)line1, (const int)len1, (const int)len2, M, X, O, E, tb);
        return;
    }
    real_seq1 += len1 - 1; real_seq2 += len2 - 1;
    int seq1_remain_length = 0, seq2_remain_length = 0;
    // seq1 will roll in last rounds
    if(line1 % block_size != 0) seq1_remain_length = int(line1 % block_size);
    // seq2 will roll as same
    if(len2 % block_size != 0)  seq2_remain_length = int(len2 % block_size);
    int block1 = line1 / block_size, block2 = len2 / block_size, block2_remain_start = block2 * block_size;

    dp_v_t b1(block_size), I1(block_size);
    dp_b_t b2(block2, block_size), D2(block2, block_size);
    dp_v_t bleft(block2 + 1);
    bleft[0] = 0; b1[0] = tb + E; I1[0] = O - (O + E);
    for(int i = 1; i < block_size; ++ i) { b1[i] = E; I1[i] = O - (O + E); }
    for(int block2_ = 0; block2_ < block2; ++ block2_)
    {
        auto &C_block = b2[block2_], &D_block = D2[block2_];
        for(int i = 0; i < block_size; ++ i)
        {
            C_block[i] = E; D_block[i] = O - (O + E);
        }
        bleft[block2_ + 1] = O + E * (block2_ + 1) * block_size;
    }
    b2[0][0] = O + E;
    if(seq2_remain_length)
    {
        for(int i = block2_remain_start + 1; i <= len2; ++ i)
        {
            final_C[i] = E; final_D[i] = O - (O + E);
        }
    }
#ifdef RUSSIAN_DP_DEBUG
    fprintf(stderr, "b1 = "); print_vec(b1.data(), block_size - 1);
    fprintf(stderr, "b2 = "); print_vec(final_C, block2_remain_start);
    print_vec(bleft.data(), block2 + 1);
#endif
    affine_input_t now_in;
    affine_output_t now_out;
    std::string this_block_s((2 * block_size), 0);
    dp_t dp_before, dp_next;
    for(int block1_ = 0; block1_ < block1; ++ block1_)
    {
        dp_before = bleft[0];
        dp_next = bleft[0];
        for(int block2_ = 0; block2_ < block2; ++ block2_)
        {
            size_t b1_first = block_size * block1_, b2_first = block_size * block2_;
            char *seq1_ = real_seq1 - b1_first, *seq2_ = real_seq2 - b2_first;
            dp_t dp_tmp;
            rev_compress_without_table(seq1_, seq2_, block_size, table_size, this_block_s);
            now_in = std::tie(this_block_s, b1, b2[block2_], I1, D2[block2_]);
#ifdef RUSSIAN_DP_DEBUG
            fprintf(stderr, "block (%zu, %zu):\n", block1_, block2_);
            fprintf(stderr, "b1 = "); print_vec(b1.data(), block_size - 1);
            fprintf(stderr, "b2 = "); print_vec(b2[block2_].data(), block_size - 1);
            fprintf(stderr, "I1 = "); print_vec(I1.data(), block_size - 1);
            fprintf(stderr, "D2 = "); print_vec(D2[block2_].data(), block_size - 1);
            fprintf(stderr, "bleft = %lld\n", bleft[block2_]);
            for(int q = 0; -q < block_size; -- q) fputc(seq1_[q] + '0', stderr); fputc('\n', stderr);
            for(int q = 0; -q < block_size; -- q) fputc(seq2_[q] + '0', stderr); fputc('\n', stderr);
#endif
            if(F -> find(now_in, now_out))
            {
                // found in table F
                std::tie(b1, b2[block2_], I1, D2[block2_], dp_tmp) = now_out;
            }
            else
            {
                rev_block_affine_dp_russian_final(block_size, b1, b2[block2_], I1, D2[block2_], seq1_, seq2_, M, X, O, E, &dp_tmp);
                now_out = std::tie(b1, b2[block2_], I1, D2[block2_], dp_tmp);
                F -> insert(now_in, now_out);
            }
            dp_next = bleft[block2_] + dp_tmp;
            bleft[block2_] = dp_before;
            dp_before = dp_next;
#ifdef RUSSIAN_DP_DEBUG
            fprintf(stderr, "=>\nbleft = %lld\n", dp_next);
            fprintf(stderr, "b1 = "); print_vec(b1.data(), block_size - 1);
            fprintf(stderr, "b2 = "); print_vec(b2[block2_].data(), block_size - 1);
            fprintf(stderr, "I1 = "); print_vec(I1.data(), block_size - 1);
            fprintf(stderr, "D2 = "); print_vec(D2[block2_].data(), block_size - 1);
            fputs("---------\n", stderr);
#endif
        }
        bleft[block2] = dp_next;
#ifdef RUSSIAN_DP_DEBUG
        print_vec(bleft.data(), block2);
#endif
        if(seq2_remain_length)
        {
            rev_block_affine_dp_russian_seq2(b1.data(), final_C + block2_remain_start + 1, I1.data(), final_D + block2_remain_start + 1, block_size, seq2_remain_length, real_seq1 - block_size * block1_, real_seq2 - block2_remain_start, M, X, O, E);
        }
        for(int i = 0; i < block_size; ++ i) { b1[i] = E; I1[i] = O - (O + E); } // next round
        bleft[0] = tb + E * block_size * (block1_ + 1);
    }
    final_C[0] = bleft[0];
    // save to C and D
    int now_block = 0, b2table = 0;
    for(int i = 1; i <= block2_remain_start; ++ i)
    {
        final_C[i] = final_C[i - 1] + b2[b2table][i - now_block - 1];
        final_D[i] = final_C[i] + (O + E) + D2[b2table][i - now_block - 1];
        if(i % block_size == 0) now_block += block_size, ++ b2table;
    }
    for(int i = block2_remain_start + 1; i <= len2; ++ i)
    {
        final_C[i] += final_C[i - 1];
        final_D[i] += final_C[i] + (O + E);
    }
    decltype(bleft)().swap(bleft);
    decltype(b1)().swap(b1);
    decltype(I1)().swap(I1);
    decltype(b2.dp_block)().swap(b2.dp_block);
    decltype(D2.dp_block)().swap(D2.dp_block);
#ifdef RUSSIAN_DP_DEBUG
    fprintf(stderr, "C = "); for(size_t i = 0; i <= len2; ++ i) fprintf(stderr, "%lld ", final_C[i]); fputc('\n', stderr);
    fprintf(stderr, "D = "); for(size_t i = 0; i <= len2; ++ i) fprintf(stderr, "%lld ", final_D[i]); fputc('\n', stderr);
#endif
    // process remain in seq1
    // dp_before defined before
    int now_seq_1 = line1 - seq1_remain_length;
    dp_t final_I1;
    while(seq1_remain_length --)
    {
        dp_t c, s;
        s = final_C[0];
        final_C[0] = c = final_C[0] + E; // must have block, so no need to add tb
        final_I1 = final_C[0] + O;
        for(int i = 1; i <= len2; ++ i)
        {
            final_I1 = std::min(final_I1, c + O) + E;
            final_D[i] = std::min(final_D[i], final_C[i] + O) + E;
            c = std::min(std::min(final_D[i], final_I1), s + MATCH_SCORE(*(real_seq1 - now_seq_1), *(real_seq2 - i + 1), M, X));
            s = final_C[i];
            final_C[i] = c;
        }
        ++ now_seq_1;
#ifdef RUSSIAN_DP_DEBUG
        fprintf(stderr, "C = "); for(size_t i = 0; i <= len2; ++ i) fprintf(stderr, "%lld ", final_C[i]); fputc('\n', stderr);
        fprintf(stderr, "D = "); for(size_t i = 0; i <= len2; ++ i) fprintf(stderr, "%lld ", final_D[i]); fputc('\n', stderr);
#endif
    }
    final_D[0] = final_C[0];
}
