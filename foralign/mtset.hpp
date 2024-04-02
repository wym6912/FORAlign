#ifndef __MUTLITHREAD_SET__
#define __MUTLITHREAD_SET__

#include <atomic>
#include <thread>
#include <shared_mutex>
#include <mutex>
#include <vector>
#include <cassert>
#include <algorithm>
#include <cstddef>

typedef struct hirschberg_status_t
{
    std::pair <int, int> seq1, seq2; // first is ID, second is __len__
    bool operator < (const hirschberg_status_t& b) const
    {
        return seq1.first != b.seq1.first ? seq1.first < b.seq1.first : (seq2.first != b.seq2.first ? seq2.first < b.seq2.first : (assert(false), false));
    }
    bool operator == (const hirschberg_status_t& b) const
    {
        return seq1.first == b.seq1.first && seq2.first == b.seq2.first;
    }
    bool operator > (const hirschberg_status_t& b) const
    {
        return !(*this == b || *this < b);
    }
    hirschberg_status_t() {}
    hirschberg_status_t(int& a, int ad, int& b, int bd) : seq1(a, ad), seq2(b, bd) {}
} hirschberg_status_t;

class ConcurrentSet_hirschberg_status_t
{
private:
    struct Node
    {
        hirschberg_status_t *value;
        std::atomic<Node*> next;

        Node(hirschberg_status_t *val) : value(val), next(nullptr) {}
    };

    std::vector<std::atomic<Node*>> *table;
    int capacity;
    std::atomic<size_t> set_size;  // 添加一个原子计数器
    std::mutex* mutexes;  // 读写锁数组
    std::atomic<bool> need_sorted;

public:
    ConcurrentSet_hirschberg_status_t(int cap = 23) : capacity(cap)
    {
        table = new std::vector<std::atomic<Node*>>(capacity);
        for (int i = 0; i < capacity; ++i) {
            std::atomic_init(&(*table)[i], nullptr);
        }
        mutexes = new std::mutex[capacity];  // 初始化读写锁数组
        set_size.store(0);
        need_sorted.store(false);
    }

    ~ConcurrentSet_hirschberg_status_t() {
        delete table;      // 在析构函数中释放内存
        delete[] mutexes;  // 删除读写锁数组
    }

    void insert(int& a, int ad, int& b, int bd);
    bool getFinalSnapshotAndFree(std::vector<hirschberg_status_t*> &v);
    
    size_t getSize();
};


#endif