#ifndef _ours_LRU_H
#define _ours_LRU_H

#include <queue>
#include "LRU.h"
#include "Periodic/Algorithm/Ours.h"
using namespace std;

struct cmp
{
    bool operator()(ItemPair x, ItemPair y)
    {
        return x.time > y.time;
    }
};

class ours_lru : public lru
{
private:
    uint64_t last_timestamp;
    Ours<2> *o;
    unordered_map<uint32_t, uint64_t> periodic_item;
    typedef unordered_multimap<uint32_t, uint64_t> mymap;
    mymap periodic_interval;
    priority_queue<ItemPair, vector<ItemPair>, cmp> heap;
    uint64_t T, advance;
    int ours_length;

public:
    int count = 0;
    ours_lru(int len, uint64_t _T, uint64_t _advance, int _ours_length) : T(_T), advance(_advance), ours_length(_ours_length)
    {
        l = list<uint32_t>(len);
        o = new Ours<2>(ours_length);
        last_timestamp = 0;
    }

    void access(uint32_t x, uint64_t t)
    {
        if (t > last_timestamp + T)
        {
            last_timestamp = t;
            periodic_item.clear();
            periodic_interval.clear();
            HashMap ans = o->Report(5);
            for (HashMap::iterator it = ans.begin(); it != ans.end(); ++it)
            {
                if (it->first.time == 0)
                {
                    continue;
                }
                //printf("%u %lu\n", it->first.item, it->first.time);
                periodic_interval.insert(make_pair(it->first.item, it->first.time));
                periodic_item[it->first.item] = 0;
            }
            delete o;
            o = new Ours<2>(ours_length);
        }

        while (!heap.empty() && heap.top().time <= t)
        {
            ItemPair p = heap.top();
            heap.pop();
            if (m.find(p.item) == m.end())
            {
                ++count;
                if (l.get_head()->v)
                    m.erase(l.get_head()->value);
                l.del(l.get_head());
                m[p.item] = l.insert(p.item, true);
            }
            else
            {
                l.del(m[p.item]);
                m[p.item] = l.insert(p.item, true);
            }
        }

        if (periodic_item.find(x) != periodic_item.end())
        {
            periodic_item[x] = t;
            pair<mymap::iterator, mymap::iterator> pr = periodic_interval.equal_range(x);
            for (; pr.first != pr.second; ++pr.first)
                heap.push(ItemPair(t + pr.first->second - advance, x));
        }

        o->Insert(ItemPair(t, x));

        ++tot_access;
        if (m.find(x) == m.end())
        {
            if (l.get_head()->v)
                m.erase(l.get_head()->value);
            l.del(l.get_head());
            m[x] = l.insert(x, true);
            ++failed_access;
        }
        else
        {
            l.del(m[x]);
            m[x] = l.insert(x, true);
            ++success_access;
        }
    }
};

#endif