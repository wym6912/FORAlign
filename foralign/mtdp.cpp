#include "mtdp.hpp"
#include "russians.hpp"


#define MATCH_SCORE(a, b, M, X) ((a) == (b) ? (M) : (X))

void linear_dp_hirschberg_multithread(int threads, const char* seq1, const char* seq2, dp_t *CC, dp_t *DD, const int line1, const int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t &tg, std::mutex *mtx, BS::thread_pool_light *pool)
{
    if(threads <= 1)
    {
        linear_dp_hirschberg(seq1, seq2, CC, DD, line1, len2, M, X, O, E, tg);
        return;
    }
    // Parameters
    int *now_line = new int[len2 + 1]();
    // Init DP
    dp_t t0;
    CC[0] = 0;
    t0 = O;
    for(int j = 1; j <= len2; ++ j)
    {
        t0 = t0 + E;
        CC[j] = t0;
        DD[j] = t0 + O;
    }
    t0 = tg; // [*]

    std::mutex *cv_mutex = new std::mutex[threads];
    std::condition_variable *cv = new std::condition_variable[threads];

    auto every_thread =
    [&line1, &len2, &threads, &t0, &E, &O, &CC, &DD, &seq1, &seq2, &M, &X, &now_line]
    (const int &thread_id, std::condition_variable &cv_before, std::mutex &cv_mutex_before, std::condition_variable &cv_this, std::mutex& cv_mutex_this)
    {
        dp_t t, s, e, c;
        std::unique_lock<std::mutex> lock_before(cv_mutex_before, std::defer_lock), lock_this(cv_mutex_this, std::defer_lock);
        for(int i = thread_id; i <= line1; i += threads)
        {
            if(now_line[0] != i - 1)
            {
                lock_before.lock();
                cv_before.wait(lock_before, [&now_line, &i](){ return now_line[0] == i - 1; }); // wait pre thread
                lock_before.unlock();
            }
            // calc dp[i][0]
            s = CC[0];
            CC[0] = c = t = t0 + E * i;
            e = t + O;
            // calc dp[i][0] finish

            lock_this.lock();
            now_line[0] = i; // dp[i][0] is ready
            lock_this.unlock();

            cv_this.notify_one();
            for(int j = 1; j <= len2; ++ j)
            {
                if(now_line[j] != i - 1)
                {
                    lock_before.lock();
                    cv_before.wait(lock_before, [&now_line, &i, &j]() { return now_line[j] == i - 1; });  // wait pre thread
                    lock_before.unlock();
                }
                // calc dp[i][j]
                e = std::min(e, c + O) + E;
                DD[j] = std::min(DD[j], CC[j] + O) + E;
                c = std::min(std::min(DD[j], e), s + (seq1[i - 1] == seq2[j - 1] ? M : X));
                s = CC[j];
                CC[j] = c;
                // calc dp[i][0] finish

                lock_this.lock();
                now_line[j] = i; // dp[i][j] is ready
                lock_this.unlock();

                cv_this.notify_one();
            }
        }
    };
    {
        std::lock_guard <std::mutex> lock(*mtx);
        // Compute DP
        for (int i = 0; i < threads - 1; ++ i)
            pool->push_task(every_thread, i + 1, std::ref(cv[i]), std::ref(cv_mutex[i]), std::ref(cv[i + 1]), std::ref(cv_mutex[i + 1]));
        pool->push_task(every_thread, threads, std::ref(cv[threads - 1]), std::ref(cv_mutex[threads - 1]), std::ref(cv[0]), std::ref(cv_mutex[0]));
        cv[0].notify_one();  // start first thread
        pool->wait_for_tasks();
    }

    DD[0] = CC[0];
    delete[] cv_mutex; delete[] cv;
    delete[] now_line;
}

void rev_linear_dp_hirschberg_multithread(int threads, const char* seq1, const char* seq2, dp_t *CC, dp_t *DD, const int line1, const int len1, const int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t &tg, std::mutex *mtx, BS::thread_pool_light *pool)
{
    if(threads <= 1)
    {
        rev_linear_dp_hirschberg(seq1, seq2, CC, DD, line1, len1, len2, M, X, O, E, tg);
        return;
    }
    // Parameters
    int *now_line = new int[len2 + 1]();
    // Init DP
    dp_t t0;
    CC[0] = 0;
    t0 = O;
    for(int j = 1; j <= len2; ++ j)
    {
        t0 = t0 + E;
        CC[j] = t0;
        DD[j] = t0 + O;
    }
    t0 = tg; // [*]

    std::mutex *cv_mutex = new std::mutex[threads];
    std::condition_variable *cv = new std::condition_variable[threads];

    auto every_thread =
    [&line1, &len1, &len2, &threads, &t0, &E, &O, &CC, &DD, &seq1, &seq2, &M, &X, &now_line]
    (const int& thread_id, std::condition_variable& cv_before, std::mutex& cv_mutex_before, std::condition_variable& cv_this, std::mutex& cv_mutex_this)
    {
        dp_t t, s, e, c;
        std::unique_lock<std::mutex> lock_before(cv_mutex_before, std::defer_lock), lock_this(cv_mutex_this, std::defer_lock);
        for(int i = thread_id; i <= line1; i += threads)
        {
            if(now_line[0] != i - 1)
            {
                lock_before.lock();
                cv_before.wait(lock_before, [&now_line, &i](){ return now_line[0] == i - 1; }); // wait pre thread
                lock_before.unlock();
            }
            // calc dp[i][0]
            s = CC[0];
            CC[0] = c = t = t0 + E * i;
            e = t + O;
            // calc dp[i][0] finish

            lock_this.lock();
            now_line[0] = i; // dp[i][0] is ready
            lock_this.unlock();

            cv_this.notify_one();
            for(int j = 1; j <= len2; ++ j)
            {
                if(now_line[j] != i - 1)
                {
                    lock_before.lock();
                    cv_before.wait(lock_before, [&now_line, &i, &j]() { return now_line[j] == i - 1; });  // wait pre thread
                    lock_before.unlock();
                }
                // calc dp[i][j]
                e = std::min(e, c + O) + E;
                DD[j] = std::min(DD[j], CC[j] + O) + E;
                c = std::min(std::min(DD[j], e), s + (seq1[len1 - i] == seq2[len2 - j] ? M : X));
                s = CC[j];
                CC[j] = c;
                // calc dp[i][0] finish

                lock_this.lock();
                now_line[j] = i; // dp[i][j] is ready
                lock_this.unlock();

                cv_this.notify_one();
            }
        }
    };
    {
        std::lock_guard <std::mutex> lock(*mtx);
        // Compute DP
        for (int i = 0; i < threads - 1; ++ i)
            pool->push_task(every_thread, i + 1, std::ref(cv[i]), std::ref(cv_mutex[i]), std::ref(cv[i + 1]), std::ref(cv_mutex[i + 1]));
        pool->push_task(every_thread, threads, std::ref(cv[threads - 1]), std::ref(cv_mutex[threads - 1]), std::ref(cv[0]), std::ref(cv_mutex[0]));
        cv[0].notify_one();  // start first thread
        pool->wait_for_tasks();
    }

    DD[0] = CC[0];
    delete[] cv_mutex; delete[] cv;
    delete[] now_line;
}

void linear_russians_affine_multithread(int threads, const char *real_seq1, const int line1, const char *real_seq2, const int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t &tb, int &block_size, F_affine_t *F, dp_t* final_C, dp_t* final_D, const int &table_size, std::mutex *mtx, BS::thread_pool_light *pool)
{
    if(line1 < block_size || len2 < block_size)
    {
        // do simple dp, no need to do russians
        linear_dp_hirschberg_multithread(threads, real_seq1, real_seq2, final_C, final_D, line1, len2, M, X, O, E, tb, mtx, pool);
        return;
    }
    int seq1_remain_length = 0, seq2_remain_length = 0;
    // seq1 will roll in last rounds
    if(line1 % block_size != 0) seq1_remain_length = int(line1 % block_size);
    // seq2 will roll as same
    if(len2 % block_size != 0)  seq2_remain_length = int(len2 % block_size);
    int block1 = line1 / block_size, block2 = len2 / block_size, block2_remain_start = block2 * block_size;
    int *now_line = new int[block2 + 2]();
    // init russians
    dp_b_t b1(threads, block_size), I1(threads, block_size);
    dp_b_t b2(block2, block_size), D2(block2, block_size);
    dp_v_t bleft(block2 + 1);
    bleft[0] = 0; b1[0][0] = tb + E; I1[0][0] = O - (O + E);
    for(int j = 1; j < block_size; ++ j) { b1[0][j] = E; I1[0][j] = O - (O + E); }
    for(int i = 1; i < threads; ++ i)
        for(int j = 0; j < block_size; ++ j)
        {
            b1[i][j] = E;
            I1[i][j] = O - (O + E);
        }
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
    // init threads
    std::mutex *cv_mutex = new std::mutex[threads];
    std::condition_variable *cv = new std::condition_variable[threads];
    auto every_thread =
    [&](const int &thread_id, std::condition_variable &cv_before, std::mutex &cv_mutex_before, std::condition_variable &cv_this, std::mutex& cv_mutex_this,
        dp_v_t &b1, dp_v_t &I1)
    {
        affine_input_t now_in;
        affine_output_t now_out;
        std::string this_block_s((2 * block_size), 0);
        dp_t dp_before, dp_next, dp_tmp;
        std::unique_lock<std::mutex> lock_before(cv_mutex_before, std::defer_lock), lock_this(cv_mutex_this, std::defer_lock);
        for(int block1_ = thread_id; block1_ < block1; block1_ += threads)
        {
            if(now_line[0] != block1_)
            {
                lock_before.lock();
                cv_before.wait(lock_before, [&](){ return now_line[0] == block1_; });
                lock_before.unlock();
            }
            bleft[0] = tb + E * block_size * block1_;
            dp_before = bleft[0];
            dp_next = bleft[0];

            lock_this.lock();
            now_line[0] = block1_ + 1;
            lock_this.unlock();

            cv_this.notify_one();
            for(int block2_ = 0; block2_ < block2; ++ block2_)
            {
                if(now_line[block2_ + 1] != block1_)
                {
                    lock_before.lock();
                    cv_before.wait(lock_before, [&]() { return now_line[block2_ + 1] == block1_; });  // wait pre thread
                    lock_before.unlock();
                }
                size_t b1_first = block_size * block1_, b2_first = block_size * block2_;
                const char *seq1_ = real_seq1 + b1_first, *seq2_ = real_seq2 + b2_first;
                compress_without_table(seq1_, seq2_, block_size, table_size, this_block_s);
                now_in = std::tie(this_block_s, b1, b2[block2_], I1, D2[block2_]);
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

                lock_this.lock();
                now_line[block2_ + 1] = block1_ + 1; // dp[i][j] is ready
                lock_this.unlock();

                cv_this.notify_one();
            }
            if(now_line[block2 + 1] != block1_)
            {
                lock_before.lock();
                cv_before.wait(lock_before, [&]() { return now_line[block2 + 1] == block1_; });  // wait pre thread
                lock_before.unlock();
            }
            bleft[block2] = dp_next;
            if(seq2_remain_length)
            {
                block_affine_dp_russian_seq2(b1.data(), final_C + block2_remain_start + 1, I1.data(), final_D + block2_remain_start + 1, block_size, seq2_remain_length, real_seq1 + block_size * block1_, real_seq2 + block2_remain_start, M, X, O, E);
            }

            lock_this.lock();
            now_line[block2 + 1] = block1_ + 1; // dp[i][j] is ready
            lock_this.unlock();

            cv_this.notify_one();
            for(int i = 0; i < block_size; ++ i) { b1[i] = E; I1[i] = O - (O + E); } // next round
        }
    };

    {
        std::lock_guard <std::mutex> lock(*mtx);
        // Compute DP
        for (int i = 0; i < threads - 1; ++ i)
            pool->push_task(every_thread, i, std::ref(cv[i]), std::ref(cv_mutex[i]), std::ref(cv[i + 1]), std::ref(cv_mutex[i + 1]), std::ref(b1[i]), std::ref(I1[i]));
        pool->push_task(every_thread, threads - 1, std::ref(cv[threads - 1]), std::ref(cv_mutex[threads - 1]), std::ref(cv[0]), std::ref(cv_mutex[0]), std::ref(b1[threads - 1]), std::ref(I1[threads - 1]));
        cv[0].notify_one();  // start first thread
        pool->wait_for_tasks();
    }
    final_C[0] = tb + E * block_size * block1;
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
    decltype(b1.dp_block)().swap(b1.dp_block);
    decltype(I1.dp_block)().swap(I1.dp_block);
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
    delete[] cv; delete[] cv_mutex;
    delete[] now_line;
}

void rev_linear_russians_affine_multithread(int threads, const char *real_seq1, const int line1, const int len1, const char *real_seq2, const int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t tb, int &block_size, F_affine_t *F, dp_t* final_C, dp_t* final_D, const int &table_size, std::mutex *mtx, BS::thread_pool_light *pool)
{
    if(line1 < block_size || len2 < block_size)
    {
        // do simple dp, no need to do russians
        rev_linear_dp_hirschberg_multithread(threads, real_seq1, real_seq2, final_C, final_D, line1, len1, len2, M, X, O, E, tb, mtx, pool);
        return;
    }
    real_seq1 += len1 - 1; real_seq2 += len2 - 1;
    int seq1_remain_length = 0, seq2_remain_length = 0;
    // seq1 will roll in last rounds
    if(line1 % block_size != 0) seq1_remain_length = int(line1 % block_size);
    // seq2 will roll as same
    if(len2 % block_size != 0)  seq2_remain_length = int(len2 % block_size);
    int block1 = line1 / block_size, block2 = len2 / block_size, block2_remain_start = block2 * block_size;
    int *now_line = new int[block2 + 2]();
    // init russians
    dp_b_t b1(threads, block_size), I1(threads, block_size);
    dp_b_t b2(block2, block_size), D2(block2, block_size);
    dp_v_t bleft(block2 + 1);
    bleft[0] = 0; b1[0][0] = tb + E; I1[0][0] = O - (O + E);
    for(int j = 1; j < block_size; ++ j) { b1[0][j] = E; I1[0][j] = O - (O + E); }
    for(int i = 1; i < threads; ++ i)
        for(int j = 0; j < block_size; ++ j)
        {
            b1[i][j] = E;
            I1[i][j] = O - (O + E);
        }
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
    // init threads
    std::mutex *cv_mutex = new std::mutex[threads];
    std::condition_variable *cv = new std::condition_variable[threads];
    auto every_thread =
    [&](const int &thread_id, std::condition_variable &cv_before, std::mutex &cv_mutex_before, std::condition_variable &cv_this, std::mutex& cv_mutex_this,
        dp_v_t &b1, dp_v_t &I1)
    {
        affine_input_t now_in;
        affine_output_t now_out;
        std::string this_block_s((2 * block_size), 0);
        dp_t dp_before, dp_next, dp_tmp;
        std::unique_lock<std::mutex> lock_before(cv_mutex_before, std::defer_lock), lock_this(cv_mutex_this, std::defer_lock);
        for(int block1_ = thread_id; block1_ < block1; block1_ += threads)
        {
            if(now_line[0] != block1_)
            {
                lock_before.lock();
                cv_before.wait(lock_before, [&](){ return now_line[0] == block1_; });
                lock_before.unlock();
            }
            bleft[0] = tb + E * block_size * block1_;
            dp_before = bleft[0];
            dp_next = bleft[0];

            lock_this.lock();
            now_line[0] = block1_ + 1;
            lock_this.unlock();

            cv_this.notify_one();
            for(int block2_ = 0; block2_ < block2; ++ block2_)
            {
                if(now_line[block2_ + 1] != block1_)
                {
                    lock_before.lock();
                    cv_before.wait(lock_before, [&]() { return now_line[block2_ + 1] == block1_; });  // wait pre thread
                    lock_before.unlock();
                }
                size_t b1_first = block_size * block1_, b2_first = block_size * block2_;
                const char *seq1_ = real_seq1 - b1_first, *seq2_ = real_seq2 - b2_first;
                rev_compress_without_table(seq1_, seq2_, block_size, table_size, this_block_s);
                now_in = std::tie(this_block_s, b1, b2[block2_], I1, D2[block2_]);
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

                lock_this.lock();
                now_line[block2_ + 1] = block1_ + 1; // dp[i][j] is ready
                lock_this.unlock();

                cv_this.notify_one();
            }
            if(now_line[block2 + 1] != block1_)
            {
                lock_before.lock();
                cv_before.wait(lock_before, [&]() { return now_line[block2 + 1] == block1_; });  // wait pre thread
                lock_before.unlock();
            }
            bleft[block2] = dp_next;
            if(seq2_remain_length)
            {
                rev_block_affine_dp_russian_seq2(b1.data(), final_C + block2_remain_start + 1, I1.data(), final_D + block2_remain_start + 1, block_size, seq2_remain_length, real_seq1 - block_size * block1_, real_seq2 - block2_remain_start, M, X, O, E);
            }

            lock_this.lock();
            now_line[block2 + 1] = block1_ + 1; // dp[i][j] is ready
            lock_this.unlock();

            cv_this.notify_one();
            for(int i = 0; i < block_size; ++ i) { b1[i] = E; I1[i] = O - (O + E); } // next round
        }
    };

    {
        std::lock_guard <std::mutex> lock(*mtx);
        // Compute DP
        for (int i = 0; i < threads - 1; ++ i)
            pool->push_task(every_thread, i, std::ref(cv[i]), std::ref(cv_mutex[i]), std::ref(cv[i + 1]), std::ref(cv_mutex[i + 1]), std::ref(b1[i]), std::ref(I1[i]));
        pool->push_task(every_thread, threads - 1, std::ref(cv[threads - 1]), std::ref(cv_mutex[threads - 1]), std::ref(cv[0]), std::ref(cv_mutex[0]), std::ref(b1[threads - 1]), std::ref(I1[threads - 1]));
        cv[0].notify_one();  // start first thread
        pool->wait_for_tasks();
    }

    final_C[0] = tb + E * block_size * block1;
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
    decltype(b1.dp_block)().swap(b1.dp_block);
    decltype(I1.dp_block)().swap(I1.dp_block);
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
    delete[] cv; delete[] cv_mutex;
    delete[] now_line;
}
