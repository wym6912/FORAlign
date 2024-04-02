#ifndef __RUSSIANS_ALL__
#define __RUSSIANS_ALL__

#include "dp.hpp"
#include "mtmap.hpp"
#include "hirschberg.h"
#include "mtdp.hpp"
#include "../include/dpthreadpool/BS_thread_pool_light.hpp"

typedef struct dp_b_t
{
    dp_b_t() {}
    dp_b_t(size_t dp_sz) { dp_block.resize(dp_sz); }
    dp_b_t(size_t dp_sz, size_t block_sz) { dp_block.resize(dp_sz); for (auto& every : dp_block) every.resize(block_sz); }
    std::vector <dp_v_t> dp_block;
    dp_v_t &operator[](size_t sz) { return dp_block[sz]; }
} dp_b_t;

typedef std::tuple <size_t, dp_v_t, dp_v_t> simple_input_t; // block_str, C1, C2
typedef std::tuple <dp_v_t, dp_v_t, dp_t>   simple_output_t; // C1', C2', delta_ans
typedef std::tuple <std::string, dp_v_t, dp_v_t, dp_v_t, dp_v_t> affine_input_t; // block_str, C1, C2, I1, D2
typedef std::tuple <dp_v_t, dp_v_t, dp_v_t, dp_v_t, dp_t>   affine_output_t; // C1', C2', I1', D2', delta_ans

// define hash function
namespace std
{
    template<>
    struct hash<simple_input_t>
    {
        typedef simple_input_t argument_type;
        typedef std::size_t result_type;

        result_type operator()(argument_type const& h) const
        {
            result_type const h1 ( std::hash<std::size_t>()(std::get<0>(h)) );
            return h1;
        }
    };

    template<>
    struct hash<affine_input_t>
    {
        typedef affine_input_t argument_type;
        typedef std::size_t result_type;

        result_type operator()(argument_type const& h) const
        {
            result_type const h1 ( std::hash<std::string>()(std::get<0>(h)) );
            result_type h2 = 0;
            std::hash<dp_t> hasher;
            for (auto i : std::get<1>(h))
            {
                h2 ^= hasher(i) + (h2 << 6) + (h2 >> 2);
            }
            for (auto i : std::get<2>(h))
            {
                h2 ^= hasher(i) + (h2 << 6) + (h2 >> 2);
            }
            return h1 ^ h2;
        }
    };

}

typedef ConcurrentHashMap <simple_input_t, simple_output_t> F_simple_t;
typedef ConcurrentHashMap <affine_input_t, affine_output_t> F_affine_t;

char* russian_table_init(sequence_t sequence_type, char* table, int &table_size);

void linear_block_affine_main_dp(char *real_seq1, int line1, char *real_seq2, int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t tb, int &block_size, F_affine_t *F, dp_t* final_C, dp_t* final_D, int &table_size);
void rev_linear_block_affine_main_dp(char *real_seq1, int line1, int len1, char *real_seq2, int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t &tb, int &block_size, F_affine_t *F, dp_t* final_C, dp_t* final_D, int &table_size);
size_t compress_without_table(const char *seq1, const char *seq2, int len, int table_size);
size_t rev_compress_without_table(const char *seq1, const char *seq2, int len, int table_size);
void compress_without_table(const char *seq1, const char *seq2, int len, int table_size, std::string &ans);
void rev_compress_without_table(const char *seq1, const char *seq2, int len, int table_size, std::string &ans);

// in mtdp.hpp
void linear_russians_affine_multithread(int threads, const char* real_seq1, const int line1, const char* real_seq2, const int len2, dp_t& M, dp_t& X, dp_t& O, dp_t& E, dp_t& tb, int& block_size, F_affine_t* F, dp_t* final_C, dp_t* final_D, const int& table_size, std::mutex* mtx, BS::thread_pool_light* pool);
void rev_linear_russians_affine_multithread(int threads, const char* real_seq1, const int line1, const int len1, const char* real_seq2, const int len2, dp_t& M, dp_t& X, dp_t& O, dp_t& E, dp_t tb, int& block_size, F_affine_t* F, dp_t* final_C, dp_t* final_D, const int& table_size, std::mutex *mtx, BS::thread_pool_light *pool);


size_t compress_with_table(char *seq1, char *seq2, int len, int table_size, char* table);
void str_compress(char *seq, char *&encoded_seq, int len, char *table);

#endif
