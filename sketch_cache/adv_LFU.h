#ifndef _adv_LFU_H
#define _adv_LFU_H

#include <queue>
#include "LFU.h"
#include "../Algorithm/adv.h"
using namespace std;

class tower_lfu : public lfu
{
private:
    uint64_t last_timestamp;
    Ours_tower *o;
    unordered_map<uint32_t, uint64_t> periodic_item;
    typedef unordered_multimap<uint32_t, uint64_t> mymap;
    mymap periodic_interval;
    priority_queue<ItemPair, vector<ItemPair>, lfu_cmp> heap;
    uint64_t T, advance;
    int ours_length;

public:
    int count = 0;
    int periodic_count = 0;
    tower_lfu(int len, uint64_t _T, uint64_t _advance, int _ours_length)
        : lfu(len),  // 关键：调用父类lfu的构造函数，传入缓存大小len
          T(_T), advance(_advance), ours_length(_ours_length) 
    {
        o = new Ours_tower(ours_length,0.1,{0.5,0.5},{32,16});
        last_timestamp = 0;
    }
    ~tower_lfu()
    {
        delete o;
    }
    void access(uint32_t x, uint64_t t)
    {
        if (t > last_timestamp + T)
        {
            last_timestamp = t;
            periodic_item.clear();
            periodic_interval.clear();
            HashMap ans = o->Report1(3);
            for (HashMap::iterator it = ans.begin(); it != ans.end(); ++it)
            {
                if (it->first.time == 0)
                {
                    continue;
                }
                //printf("%lu %lu\n", it->first.item, it->first.time);
                periodic_interval.insert(make_pair(it->first.item, it->first.time));
                periodic_item[it->first.item] = 0;
            }
            delete o;
            o = new Ours_tower(ours_length,0.1,{0.5,0.5},{16,32});
        }

        while (!heap.empty() && heap.top().time <= t)
        {
            ItemPair p = heap.top();
            heap.pop();
            lfu::access(p.item, 5, false);
        }

        if (periodic_item.find(x) != periodic_item.end())
        {
            ++periodic_count;
            periodic_item[x] = t;
            pair<mymap::iterator, mymap::iterator> pr = periodic_interval.equal_range(x);
            for (; pr.first != pr.second; ++pr.first)
                heap.push(ItemPair(t + pr.first->second - advance, x));
        }

        o->Insert(ItemPair(t, x));

        lfu::access(x);
    }
};

#endif