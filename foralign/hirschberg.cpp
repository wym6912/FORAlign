#include "hirschberg.hpp"
#include "../include/threadpool/include/scheduler.hpp"
#include "hirschberg.h"
#define HIRSCHBERG_DEBUG 0
#define HIRSCHBERG_DEBUG_PRINT 0
#define HIRSCHBERG_ASSERT 0

#define MATCH_SCORE(a, b, M, X) ((a) == (b) ? (M) : (X))
#define HIRSCHBERG_GAP(N, O, E) ((N) <= 0 ? 0 : (O) + (E) * (N))

inline dp_t hirschberg_gap(dp_t &N, dp_t &O, dp_t &E)
{
    return N <= 0 ? 0 : O + E * N;
}

using namespace staccato;

class HirschbergTask
{
public:
    HirschbergTask(char* seq1, int b1, int l1, char* seq2, int b2, int l2, dp_t tb, dp_t te, dp_t &O, dp_t &E, dp_t &M, dp_t &X, bool multidp, int dp_threads,
                   thread_pool_light *dp_pool_fwd, thread_pool_light *dp_pool_rev, std::mutex *mtx_fwd, std::mutex *mtx_rev,
                   ConcurrentSet_hirschberg_status_t *&status_set, F_affine_t *F_, int &russian_size, int &encode_size, bool simple_dp) :
                   seq1(seq1), b1(b1), l1(l1), seq2(seq2), b2(b2), l2(l2), tb(tb), te(te), O(O), E(E), M(M), X(X), multidp(multidp), dp_threads(dp_threads),
                   dp_pool_fwd(dp_pool_fwd), dp_pool_rev(dp_pool_rev), mtx_fwd(mtx_fwd), mtx_rev(mtx_rev),
                   status_set(status_set), F_(F_), russian_size(russian_size), encode_size(encode_size), simple_dp(simple_dp) {}
    char *seq1, *seq2;
    int b1, l1, b2, l2, dp_threads, russian_size, encode_size;
    dp_t tb, te, O, E, M, X;
    thread_pool_light *dp_pool_fwd, *dp_pool_rev;
    std::mutex *mtx_fwd, *mtx_rev;
    ConcurrentSet_hirschberg_status_t *status_set;
    F_affine_t *F_;
    bool simple_dp, multidp;
};

class linear_dp_Task
{
public:
    linear_dp_Task(dp_t *CC, dp_t *DD, char *seq1, int line1, int b1, int l1, char *seq2, int b2, int len2, dp_t &M, dp_t &X, dp_t &O, dp_t &E, dp_t tg,
                   bool rev, bool multidp, int &dp_threads, thread_pool_light *pool, std::mutex *mtx, F_affine_t *F_, int &russian_size, int &encode_size, bool simple_dp) :
                   CC(CC), DD(DD), seq1(seq1), line1(line1), b1(b1), l1(l1), seq2(seq2), b2(b2), len2(len2), M(M), X(X), O(O), E(E), tg(tg),
                   rev(rev), multidp(multidp), dp_threads(dp_threads), pool(pool), mtx(mtx),
                   F_(F_), russian_size(russian_size), encode_size(encode_size), simple_dp(simple_dp) {}
    dp_t *CC, *DD;
    char *seq1, *seq2;
    int line1, b1, l1, b2, len2, dp_threads, russian_size, encode_size;
    dp_t M, X, O, E, tg;
    bool rev, simple_dp, multidp;
    F_affine_t *F_;
    std::mutex *mtx;
    thread_pool_light *pool;
};

union task_info_t
{
    linear_dp_Task dpTask;
    HirschbergTask HirTask;
    task_info_t(linear_dp_Task dpTask):  dpTask(dpTask)   {}
    task_info_t(HirschbergTask HirTask): HirTask(HirTask) {}
};

class TaskRunner : public staccato::task<TaskRunner>
{
public:
    TaskRunner(linear_dp_Task dpTask) : Info(dpTask)  { type_ = 1; }
    TaskRunner(HirschbergTask HirTask): Info(HirTask) { type_ = 2; }
    void execute();
    void execute_linear_dp();
    void execute_hirschberg();
private:
    task_info_t Info;
    int type_;
};

void hirschberg_multi_init(ConcurrentSet_hirschberg_status_t *&S, F_affine_t *&F, int len, thread_pool_light *&pool1, thread_pool_light *&pool2, int dp_threads, std::mutex *&mtx1, std::mutex *&mtx2, bool multidp)
{
    S = new ConcurrentSet_hirschberg_status_t(len <= 100007 ? len : 100007); // TODO: less len to find error on sorting
    F = new F_affine_t;
    if(multidp)
    {
        pool1 = new thread_pool_light(dp_threads);
        pool2 = new thread_pool_light(dp_threads);
        mtx1 = new std::mutex;
        mtx2 = new std::mutex;
    }
    else pool1 = pool2 = nullptr, mtx1 = mtx2 = nullptr;
}

void hirschberg_multi_free(ConcurrentSet_hirschberg_status_t *&S, F_affine_t *&F, thread_pool_light *pool1, thread_pool_light *pool2, std::mutex *mtx1, std::mutex *mtx2)
{
    delete S;
    delete F;
    delete pool1;
    delete pool2;
    delete mtx1;
    delete mtx2;
}

inline void del_multi(int b1, int l1, int b2, ConcurrentSet_hirschberg_status_t &S)
{
#if HIRSCHBERG_DEBUG
    fprintf(stderr, "%s seq1 = %zu %zd seq2 = %zu %zd\n", "del  ", b1, l1, b2, 0);
#endif
    S.insert(b1, l1, b2, 0);
}

inline void add_multi(int b2, int l2, int b1, ConcurrentSet_hirschberg_status_t &S)
{
#if HIRSCHBERG_DEBUG
    fprintf(stderr, "%s seq1 = %zu %zd seq2 = %zu %zd\n", "ins  ", b1, 0, b2, l2);
#endif
    S.insert(b1, 0, b2, l2);
}

inline void match_multi(int b1, int l1, int b2, int l2, ConcurrentSet_hirschberg_status_t &S)
{
#if HIRSCHBERG_DEBUG
    fprintf(stderr, "%s seq1 = %zu %zd seq2 = %zu %zd\n", "match", b1, l1, b2, l2);
#endif
    S.insert(b1, l1, b2, l2);
}

void TaskRunner::execute_linear_dp()
{
#define F(X) this->Info.dpTask.X
    dp_t *CC = F(CC), *DD = F(DD);
    char *seq1 = F(seq1), *seq2 = F(seq2);
    int line1 = F(line1), b1 = F(b1), l1 = F(l1), b2 = F(b2), len2 = F(len2), russian_size = F(russian_size), encode_size = F(encode_size);
    dp_t M = F(M), X = F(X), O = F(O), E = F(E), tg = F(tg);
    bool rev = F(rev), simple_dp = F(simple_dp), multidp = F(multidp);
    F_affine_t* F_ = F(F_);
    int dp_threads = F(dp_threads);
    thread_pool_light *pool = F(pool);
    std::mutex *mtx = F(mtx);
#undef F
    if(rev)
    {
        if (simple_dp)
        {
            if (multidp && dp_threads > 1 && (line1 >= 1000 && len2 >= 1000)) rev_linear_dp_hirschberg_multithread(dp_threads, seq1 + b1, seq2 + b2, CC, DD, line1, l1, len2, M, X, O, E, tg, mtx, pool);
            else rev_linear_dp_hirschberg(seq1 + b1, seq2 + b2, CC, DD, line1, l1, len2, M, X, O, E, tg);
        }
        else
        {
            if (multidp && dp_threads > 1 && (line1 >= 1000 && len2 >= 1000)) rev_linear_russians_affine_multithread(dp_threads, seq1 + b1, line1, l1, seq2 + b2, len2, M, X, O, E, tg, russian_size, F_, CC, DD, encode_size, mtx, pool);
            else rev_linear_block_affine_main_dp(seq1 + b1, line1, l1, seq2 + b2, len2, M, X, O, E, tg, russian_size, F_, CC, DD, encode_size);
        }
    }
    else
    {
        if(simple_dp) 
        {
            if(multidp && dp_threads > 1 && (line1 >= 1000 && len2 >= 1000)) linear_dp_hirschberg_multithread(dp_threads, seq1 + b1, seq2 + b2, CC, DD, line1, len2, M, X, O, E, tg, mtx, pool);
            else linear_dp_hirschberg(seq1 + b1, seq2 + b2, CC, DD, line1, len2, M, X, O, E, tg);
        }
        else
        {
            if (multidp && dp_threads > 1 && (line1 >= 1000 && len2 >= 1000)) linear_russians_affine_multithread(dp_threads, seq1 + b1, line1, seq2 + b2, len2, M, X, O, E, tg, russian_size, F_, CC, DD, encode_size, mtx, pool);
            else linear_block_affine_main_dp(seq1 + b1, line1, seq2 + b2, len2, M, X, O, E, tg, russian_size, F_, CC, DD, encode_size);
        }
    }
}

void TaskRunner::execute_hirschberg()
{
#define F(X) this->Info.HirTask.X
    auto b1 = F(b1), l1 = F(l1), b2 = F(b2), l2 = F(l2);
    auto tb = F(tb), te = F(te), O = F(O), E = F(E), M = F(M), X = F(X);
    auto status_set = F(status_set);
    auto F_ = F(F_);
    auto seq1 = F(seq1), seq2 = F(seq2);
    auto russian_size = F(russian_size);
    auto encode_size = F(encode_size);
    auto simple_dp = F(simple_dp), multidp = F(multidp);
    auto dp_threads = F(dp_threads);
    auto *dp_pool_fwd = F(dp_pool_fwd), *dp_pool_rev = F(dp_pool_rev);
    auto *mtx_fwd = F(mtx_fwd), *mtx_rev = F(mtx_rev);
#undef F
    /* if sequence B is empty.... */
    if(l2 <= 0)
    {
        /* if sequence A is not empty.... */
        if(l1 > 0)
        {
#define HIRSCHBERG_DEBUG 0
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "del l1: b1 = %zu, l1 = %zu (final)\n", b1, l1);
#endif
            del_multi(b1, l1, b2, *status_set);
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
            add_multi(b2, l2, b1, *status_set);
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
            del_multi(b1, 1, b2, *status_set);
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "add l2: b2 = %zu, l2 = %zu\n", b2, l2);
#endif
            add_multi(b2, l2, b1 + 1, *status_set);
        }
        else
        {
            if(midj > 1)
            {
#if HIRSCHBERG_DEBUG
                fprintf(stderr, "add l2: b2 = %zu, midj - 1 = %zu (before)\n", b2, midj - 1);
#endif
                add_multi(b2, midj - 1, b1 + l1 - 1, *status_set);
            }
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "match 1: midi = %zu, midj = %zu (in)\n", l1 + b1, midj + b2);
#endif
            match_multi(b1 + l1 - 1, 1, b2 + midj - 1, 1, *status_set);
            if(midj < l2)
            {
#if HIRSCHBERG_DEBUG
                fprintf(stderr, "add l2: b2 = %zu, l2 - midj = %zu (after)\n", b2, l2 - midj);
#endif
                add_multi(b2 + midj, l2 - midj, b1 + l1, *status_set);
            }
        }
        return;
    }
    else
    {
        int imid = l1 >> 1, jmid;
        dp_t *CC = AllocateDPVec(l2 + 1), *DD = AllocateDPVec(l2 + 1), *RR = AllocateDPVec(l2 + 1), *SS = AllocateDPVec(l2 + 1);
        linear_dp_Task dp_pre(CC, DD, seq1, imid,      b1, l1, seq2, b2, l2, M, X, O, E, tb, false, multidp, dp_threads, dp_pool_fwd, mtx_fwd, F_, russian_size, encode_size, simple_dp),
                       dp_rev(RR, SS, seq1, l1 - imid, b1, l1, seq2, b2, l2, M, X, O, E, te, true,  multidp, dp_threads, dp_pool_rev, mtx_rev, F_, russian_size, encode_size, simple_dp);
        spawn(new(child()) TaskRunner(dp_pre));
        spawn(new(child()) TaskRunner(dp_rev));
        wait();

        bool type; dp_t min_val, cur_val;
        /* find midj, such that CC[j] + RR[j] or DD[j] + SS[j] + gap is the max */
        min_val = CC[0] + RR[l2]; jmid = 0; type = false;
        // type 1
        for(int j = 0; j <= l2; ++ j)
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
        for (size_t i = 0; i < l1; ++i) fputc(seq1[i + b1] + '0', stderr); fputc('\n', stderr);
        for (size_t i = 0; i < l2; ++i) fputc(seq2[i + b2] + '0', stderr); fputc('\n', stderr);
        fprintf(stderr, "tb = "); printDPdata(tb, stderr); fputc('\n', stderr);
        fprintf(stderr, "te = "); printDPdata(te, stderr); fputc('\n', stderr);
        fprintf(stderr, "CC    DD    RR    SS\n");
        for(size_t i = 0; i <= l2; ++ i)
        {
            fprintf(stderr, "%-5lld %-5lld %-5lld %-5lld\n", CC[i], DD[i], RR[i], SS[i]);
        }
#endif
        FreeDPVec(CC); FreeDPVec(DD); FreeDPVec(RR); FreeDPVec(SS);
        //fprintf(stderr, "Use [%d, %d] space, next will use [%d, %d] and [%d, %d] place\n", b2, b2 + l2, b2, b2 + jmid, b2 + l2 + jmid + 1, b2 + l2 + l2 + jmid + 1);
        /* Conquer recursively around midpoint */
        if(! type)
        {
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "This function will choose type 1 with (%zu, %zu) => (%zu, %zu) and (%zu, %zu) => (%zu, %zu); max = %d, (i*, j*) = (%zu, %zu)\n",
                            b1, imid, b2, jmid, b1 + imid, l1 - imid, b2 + jmid, l2 - jmid, min_val, imid, jmid);
#endif
            HirschbergTask left_task (seq1, b1,        imid,      seq2, b2,        jmid,      tb, O,  O, E, M, X, multidp, dp_threads, dp_pool_fwd, dp_pool_rev, mtx_fwd, mtx_rev, status_set, F_, russian_size, encode_size, simple_dp),
                           right_task(seq1, b1 + imid, l1 - imid, seq2, b2 + jmid, l2 - jmid, O,  te, O, E, M, X, multidp, dp_threads, dp_pool_fwd, dp_pool_rev, mtx_fwd, mtx_rev, status_set, F_, russian_size, encode_size, simple_dp);
            spawn(new(child()) TaskRunner(left_task));
            spawn(new(child()) TaskRunner(right_task));
            wait();
        }
        else
        {
            del_multi(imid + b1 - 1, 2, jmid + b2, *status_set);
#if HIRSCHBERG_DEBUG
            fprintf(stderr, "This function will choose type 2 with (%zu, %zu) => (%zu, %zu), with 2 deletions, and (%zu, %zu) => (%zu, %zu); max = %d, (i*, j*) = (%zu, %zu)\n",
                            b1, imid - 1, b2, jmid, b1 + imid + 1, l1 - imid - 1, b2 + jmid, l2 - jmid, min_val, imid, jmid);
            fprintf(stderr, "del 2: midi = %zu\n", imid + b1);
#endif
#undef HIRSCHBERG_DEBUG
            HirschbergTask left_task (seq1, b1,            imid - 1,      seq2, b2,        jmid,      tb, 0,  O, E, M, X, multidp, dp_threads, dp_pool_fwd, dp_pool_rev, mtx_fwd, mtx_rev, status_set, F_, russian_size, encode_size, simple_dp),
                           right_task(seq1, b1 + imid + 1, l1 - imid - 1, seq2, b2 + jmid, l2 - jmid, 0,  te, O, E, M, X, multidp, dp_threads, dp_pool_fwd, dp_pool_rev, mtx_fwd, mtx_rev, status_set, F_, russian_size, encode_size, simple_dp);
            spawn(new(child()) TaskRunner(left_task));
            spawn(new(child()) TaskRunner(right_task));
            wait();
        }
    }
}

void TaskRunner::execute()
{
    switch (this->type_)
    {
    case 1:
        this->execute_linear_dp();
        break;
    case 2:
        this->execute_hirschberg();
        break;
    default:
        break;
    }
}

void hirschberg_multi_start(char *seq1, char *seq2, int len1, int len2, dp_t O, dp_t E, dp_t M, dp_t X, int russian_table_size,
                            ConcurrentSet_hirschberg_status_t *&status_set, F_affine_t *F_, bool multidp, int dp_threads, 
                            thread_pool_light *pool1, thread_pool_light *pool2, std::mutex *mtx1, std::mutex *mtx2,
                            int encode_table_size, int threads, bool simple_dp)
{
    scheduler<TaskRunner> ____(4, threads);
    ____.spawn(new(____.root()) TaskRunner(HirschbergTask(seq1, 0, len1, seq2, 0, len2, O, O, O, E, M, X, multidp, dp_threads, pool1, pool2, mtx1, mtx2, status_set, F_, russian_table_size, encode_table_size, simple_dp)));
    ____.wait();
}

inline void del_single(int &s, std::deque <int> &q)
{
    if(q.size() && q.back() < 0)
        q.back() -= s;
    else q.emplace_back(-s);
}

inline void add_single(int &s, std::deque <int> &q)
{
    if(q.size() && q.back() < 0)
    {
        auto tmp = q.back();
        q.back() = s;
        q.emplace_back(std::move(tmp));
    }
    else q.emplace_back(s);
}

inline void match_single(int &s, std::deque <int> &q)
{
    while(s --) q.emplace_back(0);
}


void print_result_multi(char *FILE_OUT, int len1, int len2, char *seq1, char *seq2, char *name1, char *name2, ConcurrentSet_hirschberg_status_t *S)
{
    int now1 = 0, now2 = 0, t1p = 0, t2p = 0;
#if HIRSCHBERG_ASSERT
    int sum1 = 0, sum2 = 0;
#endif
    char *tmpseq1 = AllocateCharVec(len1 + len2 + 1), *tmpseq2 = AllocateCharVec(len1 + len2 + 1);
#if HIRSCHBERG_DEBUG_PRINT
    fprintf(stderr, "multi mode\n");
#endif
    // if not priority_queue, sort the data
    t1p = t2p = 0;
    now1 = now2 = 0;
    std::vector <hirschberg_status_t*> v;
    S->getFinalSnapshotAndFree(v);
    std::deque <int> result_multi;
    for(auto &item: v)
    {
#if HIRSCHBERG_ASSERT
        assert(sum1 == item->seq1.first && sum2 == item->seq2.first);
        sum1 = item->seq1.first + item->seq1.second;
        sum2 = item->seq2.first + item->seq2.second;
#endif
#if HIRSCHBERG_PRINT
        fprintf(stderr, "%s seq1 = %zu %zu seq2 = %zu %zu\n",
                        item->seq1.second == item->seq2.second ? "match" : item->seq1.second ? "del  " : "ins  ",
                        item->seq1.first, item->seq1.second, item->seq2.first, item->seq2.second);
#endif
        if(item->seq1.second && item->seq2.second) match_single(item->seq1.second, result_multi);
        else if(item->seq1.second) del_single(item->seq1.second, result_multi);
        else add_single(item->seq2.second, result_multi);
        delete item;
    }
    decltype(v)().swap(v);
    while(! result_multi.empty())
    {
        auto &top = result_multi.front();
#if HIRSCHBERG_PRINT
        fprintf(stderr, "%d ", top);
#endif
        if(! top)
        {
            tmpseq1[t1p] = *seq1 ++;
            tmpseq2[t1p] = *seq2 ++;
            ++ t1p;
        }
        else if(top > 0)
        {
            for(int j = 0; j < top; ++ j)
            {
                tmpseq1[t1p + j] = '-';
                tmpseq2[t1p + j] = *seq2 ++;
            }
            t1p += top;
        }
        else
        {
            auto k = -top;
            for(int j = 0; j < k; ++ j)
            {
                tmpseq1[t1p + j] = *seq1 ++;
                tmpseq2[t1p + j] = '-';
            }
            t1p += k;
        }
        result_multi.pop_front();
    }
#if HIRSCHBERG_DEBUG_PRINT
    fprintf(stderr, "multi mode:\nseq1 = %s\nseq2 = %s\n", tmpseq1, tmpseq2);
#endif
    FILE *final_o = fopen(FILE_OUT, "wb+");
    if(final_o == nullptr) { fprintf(stderr, "Error: can not open file %s. Program will exit.\n", FILE_OUT); exit(1); }
    fprintf(final_o, ">%s\n%s\n>%s\n%s\n", name1, tmpseq1, name2, tmpseq2);
    FreeCharVec(tmpseq1); FreeCharVec(tmpseq2);
}

void write_final_result(int len1, int len2, char *seq1, char *seq2, ConcurrentSet_hirschberg_status_t *S, char **final_seq1, char **final_seq2)
{
    size_t now1 = 0, now2 = 0, t1p = 0;
    char *tmpseq1 = AllocateCharVec(len1 + len2 + 1), *tmpseq2 = AllocateCharVec(len1 + len2 + 1);
    // if not priority_queue, sort the data
    t1p = 0;
    now1 = now2 = 0;
    std::vector <hirschberg_status_t*> v;
    S->getFinalSnapshotAndFree(v);
    std::deque <int> result_multi;
    for(auto &item: v)
    {
        if(item->seq1.second && item->seq2.second) match_single(item->seq1.second, result_multi);
        else if(item->seq1.second) del_single(item->seq1.second, result_multi);
        else add_single(item->seq2.second, result_multi);
        delete item;
    }
    decltype(v)().swap(v);
    while(! result_multi.empty())
    {
        auto &top = result_multi.front();
        if(! top)
        {
            (tmpseq1[t1p]) = *seq1 ++;
            (tmpseq2[t1p]) = *seq2 ++;
            ++ t1p;
        }
        else if(top > 0)
        {
            for(size_t j = 0; j < top; ++ j)
            {
                (tmpseq1[t1p + j]) = '-';
                (tmpseq2[t1p + j]) = *seq2++;
            }
            t1p += top;
        }
        else
        {
            auto k = -top;
            for(size_t j = 0; j < k; ++ j)
            {
                (tmpseq1[t1p + j]) = *seq1 ++;
                (tmpseq2[t1p + j]) = '-';
            }
            t1p += k;
        }
        result_multi.pop_front();
    }
    (tmpseq1[t1p]) = (tmpseq2[t1p]) = 0;
    *final_seq1 = tmpseq1; *final_seq2 = tmpseq2;
}

void write_final_result(int len1, int len2, char* seq1, char* seq2, ConcurrentSet_hirschberg_status_t* S, char* cigar, int *final_cigar)
{
    *final_cigar = 0;
    std::vector <hirschberg_status_t*> v;
    S->getFinalSnapshotAndFree(v);
    std::deque <int> result_multi;
    for (auto& item : v)
    {
        if (item->seq1.second && item->seq2.second) match_single(item->seq1.second, result_multi);
        else if (item->seq1.second) del_single(item->seq1.second, result_multi);
        else add_single(item->seq2.second, result_multi);
        delete item;
    }
    decltype(v)().swap(v);
    while (!result_multi.empty())
    {
        auto top = result_multi.front();
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

        result_multi.pop_front();
    }
}

DLL void hirschberg_API(char *seq1, char *seq2, int len1, int len2, dp_t O, dp_t E, dp_t M, dp_t X, int table_size, sequence_t type_, int threads, char **out_seq1, char **out_seq2, short simple_dp, short multidp, int dp_threads)
{
    int encode_size;
    char *russian_table = nullptr, *comp_seq1, *comp_seq2;
    russian_table = russian_table_init(type_, russian_table, encode_size);
    ConcurrentSet_hirschberg_status_t *s;
    thread_pool_light *pool1, *pool2;
    std::mutex *mtx1, *mtx2;
    F_affine_t *F;
    str_compress(seq1, comp_seq1, len1, russian_table);
    str_compress(seq2, comp_seq2, len2, russian_table);
    hirschberg_multi_init(s, F, len1 + len2, pool1, pool2, dp_threads, mtx1, mtx2, multidp ? true : false);
    hirschberg_multi_start(comp_seq1, comp_seq2, len1, len2, O, E, M, X, table_size, s, F, multidp, dp_threads, pool1, pool2, mtx1, mtx2, encode_size, threads, simple_dp);
    //fprintf(stderr, "Set size = %zu, map size = %zu\n", s->getSize(), F->getSize());
    FreeCharVec(comp_seq1); FreeCharVec(comp_seq2); // save memory
    write_final_result(len1, len2, seq1, seq2, s, out_seq1, out_seq2);
    hirschberg_multi_free(s, F, pool1, pool2, mtx1, mtx2);
    FreeCharVec(russian_table);
}

DLL void hirschberg_cigar(char* seq1, char* seq2, int len1, int len2, dp_t O, dp_t E, dp_t M, dp_t X, int table_size, sequence_t type_, int threads, char* out_cigar, int* cigar_end, short simple_dp, short multidp, int dp_threads)
{
    int encode_size;
    char* russian_table = nullptr, * comp_seq1, * comp_seq2;
    russian_table = russian_table_init(type_, russian_table, encode_size);
    ConcurrentSet_hirschberg_status_t* s;
    thread_pool_light *pool1, *pool2;
    std::mutex *mtx1, *mtx2;
    F_affine_t* F;
    str_compress(seq1, comp_seq1, len1, russian_table);
    str_compress(seq2, comp_seq2, len2, russian_table);
    hirschberg_multi_init(s, F, len1 + len2, pool1, pool2, dp_threads, mtx1, mtx2, multidp ? true : false);
    hirschberg_multi_start(comp_seq1, comp_seq2, len1, len2, O, E, M, X, table_size, s, F, multidp, dp_threads, pool1, pool2, mtx1, mtx2, encode_size, threads, simple_dp);
    //fprintf(stderr, "Set size = %zu, map size = %zu\n", s->getSize(), F->getSize());
    FreeCharVec(comp_seq1); FreeCharVec(comp_seq2); // save memory
    write_final_result(len1, len2, seq1, seq2, s, out_cigar, cigar_end);
    hirschberg_multi_free(s, F, pool1, pool2, mtx1, mtx2);
    FreeCharVec(russian_table);
}

