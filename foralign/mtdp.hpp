#ifndef __MULTIDP__
#define __MULTIDP__

#include "dp.hpp"
#include "../include/dpthreadpool/BS_thread_pool_light.hpp"

#include <thread>
#include <mutex>
#include <vector>
#include <condition_variable>
#include <atomic>

void linear_dp_hirschberg_multithread(int threads, const char* seq1, const char* seq2, dp_t *CC, dp_t *DD, const int line1, const int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t &tg, std::mutex *mtx, BS::thread_pool_light *pool);
void rev_linear_dp_hirschberg_multithread(int threads, const char* seq1, const char* seq2, dp_t *CC, dp_t *DD, const int line1, const int len1, const int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t &tg, std::mutex *mtx, BS::thread_pool_light *pool);

#endif