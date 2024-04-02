#include "mtset.hpp"

bool less_hirschberg_t(const hirschberg_status_t* a, const hirschberg_status_t* b)
{
    return *a < *b;
}

void ConcurrentSet_hirschberg_status_t::insert(int& a, int ad, int& b, int bd)
{
    int index = (a + b) % capacity;

    hirschberg_status_t *n = new hirschberg_status_t(a, ad, b, bd);
    Node* newNode = new Node(n), *curr;
    {
        std::lock_guard<std::mutex> lock(mutexes[index]);  // 获取写锁
        curr = (*table)[index].load();
        newNode->next.store(curr);
        (*table)[index].store(newNode);
        set_size.fetch_add(1);
    }
    if (curr) need_sorted.store(true);
}

size_t ConcurrentSet_hirschberg_status_t::getSize()
{
    return set_size.load();
}

bool ConcurrentSet_hirschberg_status_t::getFinalSnapshotAndFree(std::vector<hirschberg_status_t*> &v)
{
    std::vector<std::lock_guard<std::mutex>> locks(mutexes, mutexes + capacity);
    Node *curr, *pre;
    v.reserve(set_size.load());
    for(int i = 0; i < capacity; ++ i)
    {
        curr = (*table)[i];
        while(curr != nullptr)
        {
            v.emplace_back(curr -> value);
            pre = curr;
            curr = curr->next.load();
            delete pre;
        }
    }
    if(need_sorted.load() == true) std::sort(v.begin(), v.end(), less_hirschberg_t);
    return true;
}
