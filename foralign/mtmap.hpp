#ifndef __MUTLITHREAD_MAP__
#define __MUTLITHREAD_MAP__

#include <atomic>
#include <thread>
#include <shared_mutex>
#include <mutex>
#include <vector>
#include <functional>
#include <cassert>

template <typename K, typename V>
class ConcurrentHashMap
{
private:
    struct Node
    {
        K key;
        V value;
        std::atomic<Node*> next;
        Node(const K& k, const V& v) : key(k), value(v), next(nullptr) {}
    };

    int num_buckets;
    std::vector<std::atomic<Node*>> *buckets;
    std::shared_mutex* mutexes;  // 读写锁数组
    std::atomic<size_t> map_size;  // 添加一个原子计数器

    int hash(const K& key)
    {
        return std::hash<K>{}(key) % num_buckets;
    }

public:
    ConcurrentHashMap(int num_buckets = 100007) : num_buckets(num_buckets)
    {
        buckets = new std::vector<std::atomic<Node*>>(num_buckets);
        for (int i = 0; i < num_buckets; ++i)
        {
            std::atomic_init(&(*buckets)[i], nullptr);
        }
        mutexes = new std::shared_mutex[num_buckets];  // 初始化读写锁数组
        map_size.store(0);
    }

    ~ConcurrentHashMap()
    {
        FreeTable();
        delete buckets;   // 在析构函数中释放内存
        delete[] mutexes; // 删除读写锁数组
    }

    void insert(const K& key, const V& value)
    {
        int idx = hash(key);

        std::lock_guard<std::shared_mutex> lock(mutexes[idx]);

        // 检查是否已经存在相同的键
        Node* curr = (*buckets)[idx].load(), *pre = curr;
        if(curr == nullptr)
        {
            // 插入到桶首
            (*buckets)[idx].store(new Node(key, value));
            map_size.fetch_add(1);
            return;
        }
        while (curr)
        {
            if (curr->key == key)
            {
                // No need to update
                return;
            }
            pre = curr;
            curr = curr->next.load();
        }

        // 不存在相同的键，将新节点插入链表尾部
        pre->next.store(new Node(key, value)); 
        map_size.fetch_add(1);
    }

    bool find(const K& key, V& value)
    {
        int idx = hash(key);

        std::shared_lock<std::shared_mutex> lock(mutexes[idx]);
        Node* curr = (*buckets)[idx].load();
        while (curr)
        {
            if (curr->key == key)
            {
                value = curr->value;
                return true;
            }
            curr = curr->next.load();
        }
        return false;
    }

    size_t getSize()
    {
        return map_size.load();
    }

    void getFinalSnapshotAndFree(std::vector<std::pair<K, V> > &v)
    {
        std::vector<std::unique_lock<std::shared_mutex>> locks(mutexes, mutexes + num_buckets);  // 获取所有哈希表的读锁
        Node *curr, *pre;
        v.reserve(map_size.load());
        for(int i = 0; i < num_buckets; ++ i)
        {
            curr = (*buckets)[i].load();
            while(curr != nullptr)
            {
                v.emplace_back(curr -> key, curr -> value);
                pre = curr;
                curr = curr->next;
                delete pre;
            }
        }
    }

    void FreeTable()
    {
        Node *curr, *pre;
        for(int i = 0; i < num_buckets; ++ i)
        {
            curr = (*buckets)[i].load();
            while(curr != nullptr)
            {
                pre = curr;
                curr = curr->next;
                delete pre;
            }
        }
    }
};


#endif
