#ifndef _LRU_H
#define _LRU_H

#include <unordered_map>
#include "list.h"
using namespace std;

class lru
{
protected:
    int tot_access = 0, failed_access = 0, success_access = 0;
    list<uint32_t> l;
    unordered_map<uint32_t, list_node<uint32_t> *> m;

public:
    lru(int len = 0)
    {
        l = list<uint32_t>(len);
        l.check();
    }
    void access(uint32_t x)
    {
        ++tot_access;
        if (m.find(x) == m.end())
        {
            if (l.get_head()->v){
                m.erase(l.get_head()->value);
            }
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
    double get_hit_rate() { return 1.0 * success_access / tot_access; }
    int get_tot_access() { return tot_access; }
    int get_success_access() { return success_access; }
    int get_failed_access() { return failed_access; }
    void clear_access()
    {
        tot_access = success_access = failed_access = 0;
    }
};

#endif