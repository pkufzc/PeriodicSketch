#ifndef _tower_LRU_H
#define _tower_LRU_H

#include <queue>
#include "LRU.h"
#include "../Algorithm/tower.h"
#include "list.h"
using namespace std;


class tower_lru : public lru
{
private:
    uint64_t last_timestamp;
    Ours_tower *o;
    unordered_map<uint32_t, uint64_t> periodic_item;
    typedef unordered_multimap<uint32_t, uint64_t> mymap;
    mymap periodic_interval;
    priority_queue<ItemPair, vector<ItemPair>, cmp> heap;
    uint64_t T, advance;
    int ours_length;

public:
    int count = 0;
    tower_lru(int len, uint64_t _T, uint64_t _advance, int _ours_length)
        : lru(len),  // 关键：调用父类lru的构造函数，传入缓存大小len
          T(_T), advance(_advance), ours_length(_ours_length)
    {
        // l = ::list<uint32_t>(len);
        o = new Ours_tower(ours_length,0.1,{0.5,0.5},{32,16});
        last_timestamp = 0;
        //std::cout << "tower_lru 实例创建成功，缓存大小=" << len << std::endl;
    }

    void access(uint32_t x, uint64_t t)
    {
        //std::cout << "准备调用access：key=" << x << ", ts=" << t << std::endl; 
        if (t > last_timestamp + T)
        {
            last_timestamp = t;
            periodic_item.clear();
            periodic_interval.clear();
            HashMap ans = o->Report1(3);
            static int insert_count = 0;  // 静态变量，统计总插入次数
            insert_count++;
            if (insert_count % 10000 == 0) {  // 每插入1000条数据打印一次
                //printf("Ours 已插入数据量：%d\n", insert_count);
            }
            //printf("Report 返回的周期性数据数量：%d\n", ans.size());
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
            //std::cout << "步骤1完成：Ours_tower已重置" << std::endl;
            o = new Ours_tower(ours_length,0.1,{0.5,0.5},{16,32});
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
        //std::cout << "步骤3：检查periodic_item" << std::endl;
        if (periodic_item.find(x) != periodic_item.end())
        {
            periodic_item[x] = t;
            pair<mymap::iterator, mymap::iterator> pr = periodic_interval.equal_range(x);
            for (; pr.first != pr.second; ++pr.first)
                heap.push(ItemPair(t + pr.first->second - advance, x));
        }
        //std::cout << "tower_lru::access: 准备插入 key=" << x << ", timestamp=" << t << std::endl;
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
