#include "parallel_swg.h"
#include <vector>
#include <thread>
#include <barrier>
#include <condition_variable>
#include <cstdio>

#include "../include/dpthreadpool/BS_thread_pool_light.hpp"

#define MIN(a,b) (((a)<=(b))?(a):(b))
#define AFFINE_SCORE_MAX (10000000)

using namespace BS;

thread_pool_light *swg_pool;

DLL void multithread_swg_init(int threads)
{
    swg_pool = new thread_pool_light(threads);
    if(swg_pool == nullptr) { fprintf(stderr, "FATAL ERROR: can not start thread pool. Program will exit.\n"); exit(1); }
}

DLL void multithread_swg_free()
{
    delete swg_pool;
    swg_pool = nullptr;
}

#if __cplusplus >= 202002L
DLL void multithread_swg_compute_barrier(int threads, affine_matrix_t* const affine_table, affine_penalties_t* const penalties, const char* const pattern, const int pattern_length, const char* const text, const int text_length)
{
    // Parameters
    affine_cell_t** const dp = affine_table->columns;
    int h, v;
    // Init DP
    dp[0][0].D = AFFINE_SCORE_MAX;
    dp[0][0].I = AFFINE_SCORE_MAX;
    dp[0][0].M = 0;
    for (v = 1; v <= pattern_length; ++ v)
    {   // Init first column
        dp[0][v].D = penalties->gap_opening + v * penalties->gap_extension;
        dp[0][v].I = AFFINE_SCORE_MAX;
        dp[0][v].M = dp[0][v].D;
    }
    for (h = 1; h <= text_length; ++ h)
    {   // Init first row
        dp[h][0].D = AFFINE_SCORE_MAX;
        dp[h][0].I = penalties->gap_opening + h*penalties->gap_extension;
        dp[h][0].M = dp[h][0].I;
    }
    std::barrier sync(threads);

    auto swg_thread = [&threads, &sync, &dp, &penalties, &pattern, &pattern_length, &text, &text_length]
                      (const int &thread_id)
    {
        auto dp_routine = [&](const int &h, const int &v)
        {
            // Update DP.D
            const int del_new = dp[h][v-1].M + penalties->gap_opening + penalties->gap_extension;
            const int del_ext = dp[h][v-1].D + penalties->gap_extension;
            const int del = MIN(del_new, del_ext);
            dp[h][v].D = del;
            // Update DP.I
            const int ins_new = dp[h-1][v].M + penalties->gap_opening + penalties->gap_extension;
            const int ins_ext = dp[h-1][v].I + penalties->gap_extension;
            const int ins = MIN(ins_new, ins_ext);
            dp[h][v].I = ins;
            // Update DP.M
            const int m_match = dp[h-1][v-1].M + ((pattern[v - 1] == text[h - 1]) ? penalties->match : penalties->mismatch);
            dp[h][v].M = MIN(m_match, MIN(ins, del));
        };
        if(text_length >= pattern_length)
        {
            // part 1
            for (int round = 2; round <= pattern_length; ++ round)
            {
                for(int h = thread_id; h < round; h += threads)
                {
                    const int v = round - h;
                    dp_routine(h, v);
                }
                sync.arrive_and_wait();
            }
            // part 2
            for (int round = pattern_length + 1; round <= text_length; ++ round)
            {
                for(int v = thread_id; v <= pattern_length; v += threads)
                {
                    const int h = round - v;
                    dp_routine(h, v);
                }
                sync.arrive_and_wait();
            }
            // part 3
            for (int round = 1; round <= pattern_length; ++ round)
            {
                for(int v = round + thread_id - 1; v <= pattern_length; v += threads)
                {
                    const int h = text_length + round - v;
                    dp_routine(h, v);
                }
                sync.arrive_and_wait();
            }
        }
        else
        {
            // part 1
            for (int round = 2; round <= text_length; ++ round)
            {
                for(int v = thread_id; v < round; v += threads)
                {
                    const int h = round - v;
                    dp_routine(h, v);
                }
                sync.arrive_and_wait();
            }
            // part 2
            for (int round = text_length + 1; round <= pattern_length; ++ round)
            {
                for(int h = thread_id; h <= text_length; h += threads)
                {
                    const int v = round - h;
                    dp_routine(h, v);
                }
                sync.arrive_and_wait();
            }
            // part 3
            for (int round = 1; round <= text_length; ++ round)
            {
                for(int h = round + thread_id - 1; h <= text_length; h += threads)
                {
                    const int v = pattern_length + round - h;
                    dp_routine(h, v);
                }
                sync.arrive_and_wait();
            }
        }
    };
    // Compute DP
    for (int i = 0; i < threads; ++ i)
        swg_pool->push_task(swg_thread, i + 1);
    swg_pool->wait_for_tasks();
}
#else
DLL void multithread_swg_compute_barrier(int threads, affine_matrix_t* const affine_table, affine_penalties_t* const penalties, const char* const pattern, const int pattern_length, const char* const text, const int text_length)
{
    fprintf(stderr, "Warning: not supported C++20. This function can not use. \nTry use a C++20 compiler and run again.\n");
    return;
}
#endif

DLL void multithread_swg_compute_cv(int threads, affine_matrix_t* const affine_table, affine_penalties_t* const penalties, const char* const pattern, const int pattern_length, const char* const text, const int text_length)
{
    // Parameters
    affine_cell_t** const dp = affine_table->columns;
    int *now_line = new int[pattern_length + 1]();
    int h, v;
    // Init DP
    dp[0][0].D = AFFINE_SCORE_MAX;
    dp[0][0].I = AFFINE_SCORE_MAX;
    dp[0][0].M = 0;
    for (v = 1; v <= pattern_length; ++ v)
    {   // Init first column
        dp[0][v].D = penalties->gap_opening + v * penalties->gap_extension;
        dp[0][v].I = AFFINE_SCORE_MAX;
        dp[0][v].M = dp[0][v].D;
    }
    for (h = 1; h <= text_length; ++ h)
    {   // Init first row
        dp[h][0].D = AFFINE_SCORE_MAX;
        dp[h][0].I = penalties->gap_opening + h*penalties->gap_extension;
        dp[h][0].M = dp[h][0].I;
    }
    now_line[0] = text_length; // all perpared
    std::mutex *cv_mutex = new std::mutex[threads];
    std::condition_variable *cv = new std::condition_variable[threads];

    auto swg_thread = [&threads, &dp, &now_line, &penalties, &pattern, &pattern_length, &text, &text_length]
                      (const int &thread_id, std::condition_variable &cv_before, std::mutex &cv_mutex_before, std::condition_variable &cv_this, std::mutex &cv_mutex_this)
    {
        std::unique_lock<std::mutex> lock_before(cv_mutex_before, std::defer_lock), lock_this(cv_mutex_this, std::defer_lock);
        for (int h = thread_id; h <= text_length; h += threads)
        {
            for (int v = 1; v <= pattern_length; ++ v)
            {
                while(now_line[v] != h - 1)
                {
                    lock_before.lock();
                    cv_before.wait_for(lock_before, std::chrono::nanoseconds(1), [&now_line, &v, &h]() { return now_line[v] == h - 1; });  // wait pre thread
                    lock_before.unlock();
                }
                // Update DP.D
                const int del_new = dp[h][v-1].M + penalties->gap_opening + penalties->gap_extension;
                const int del_ext = dp[h][v-1].D + penalties->gap_extension;
                const int del = MIN(del_new, del_ext);
                dp[h][v].D = del;
                // Update DP.I
                const int ins_new = dp[h-1][v].M + penalties->gap_opening + penalties->gap_extension;
                const int ins_ext = dp[h-1][v].I + penalties->gap_extension;
                const int ins = MIN(ins_new, ins_ext);
                dp[h][v].I = ins;
                // Update DP.M
                const int m_match = dp[h-1][v-1].M + ((pattern[v - 1] == text[h - 1]) ? penalties->match : penalties->mismatch);
                dp[h][v].M = MIN(m_match, MIN(ins, del));
                // lock_this.lock();
                now_line[v] = h; // dp[h][v] is ready, use wait_for to avoid lock
                // lock_this.unlock();
                cv_this.notify_all();
            }
        }
    };
    // Compute DP
    for (int i = 0; i < threads - 1; ++ i)
        swg_pool->push_task(swg_thread, i + 1, std::ref(cv[i]), std::ref(cv_mutex[i]), std::ref(cv[i + 1]), std::ref(cv_mutex[i + 1]));
    swg_pool->push_task(swg_thread, threads, std::ref(cv[threads - 1]), std::ref(cv_mutex[threads - 1]), std::ref(cv[0]), std::ref(cv_mutex[0]));

    cv[0].notify_one();  // 启动第一个线程

    swg_pool->wait_for_tasks();

    delete[] cv_mutex; delete[] cv;
    delete[] now_line;
}

DLL void multithread_swg_traceback(affine_matrix_t* const affine_matrix, affine_penalties_t* const penalties, const int pattern_length, const int text_length, char *seq1, char *seq2, char **oseq1, char **oseq2)
{ 
  // Parameters
  affine_cell_t** const dp = affine_matrix->columns;
  *oseq1 = (char*)malloc(pattern_length + text_length + 1);
  *oseq2 = (char*)malloc(pattern_length + text_length + 1);
  char *s1 = *oseq1, *s2 = *oseq2;
  if(*oseq1 == nullptr || *oseq2 == nullptr)
  {
    fprintf(stderr, "Error: can not traceback.\n");
    exit(1);
  }
  // Add final insertions/deletions
  int i;
  // Compute traceback
  affine_matrix_type matrix_type = affine_matrix_M;
  int h = pattern_length - 1;
  int v = text_length - 1;
  while (h>0 && v>0) {
    switch (matrix_type) {
      case affine_matrix_D:
        // Traceback D-matrix
        *s1 = *seq1; s1 ++; seq1 ++;
        *s2 = '-';   s2 ++;
        if (dp[h][v].D != dp[h][v-1].D + penalties->gap_extension) {
          matrix_type = affine_matrix_M;
        }
        --v;
        break;
      case affine_matrix_I:
        // Traceback I-matrix
        *s1 = '-';   s1 ++;
        *s2 = *seq2; s2 ++; seq2 ++;
        if (dp[h][v].I != dp[h-1][v].I + penalties->gap_extension) {
          matrix_type = affine_matrix_M;
        }
        --h;
        break;
      case affine_matrix_M:
        // Traceback M-matrix
        if (dp[h][v].M == dp[h-1][v-1].M + penalties->mismatch) {
          *s1 = *seq1; s1 ++;
          *s2 = *seq2; s2 ++;
          --h; --v;
        } else if (dp[h][v].M == dp[h][v].D) {
          matrix_type = affine_matrix_D;
        } else if (dp[h][v].M == dp[h][v].I) {
          matrix_type = affine_matrix_I;
        } else if (dp[h][v].M == dp[h-1][v-1].M + penalties->match) {
          *s1 = *seq1; s1 ++;
          *s2 = *seq2; s2 ++; 
          --h; --v;
        } else {
          fprintf(stderr,"SWG backtrace. No backtrace operation found\n");
          exit(1);
        }
        break;
    }
  }
  // Add initial deletions/insertions
  while (v>0) { *s1 = *seq1; s1 ++; seq1 ++; *s2 = '-'; s2 ++; --v; }
  while (h>0) { *s1 = '-'; s1 ++; *s2 = *seq2; s2 ++; seq2 ++; --h; }
}
