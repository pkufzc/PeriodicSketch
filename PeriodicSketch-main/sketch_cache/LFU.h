#ifndef _LFU_H
#define _LFU_H

#include <unordered_map>
#include "list.h"
using namespace std;

class lfu
{
protected:
    int max_number, total_number = 0;
    int tot_access = 0, failed_access = 0, success_access = 0;
    list<list<pair<uint32_t, uint32_t>>> l;
    unordered_map<uint32_t, list_node<list<pair<uint32_t, uint32_t>>> *> h1;
    unordered_map<uint32_t, list_node<pair<uint32_t, uint32_t>> *> h2;

public:
    lfu(int x) : max_number(x) {}
    virtual ~lfu()
    {
    }
    void access(uint32_t x, uint32_t y = 1, bool op = true)
    {
        if (op)
            ++tot_access;
        if (h2.find(x) != h2.end())
        {
            if (op)
                ++success_access;
            list_node<pair<uint32_t, uint32_t>> *p = h2[x];
            uint32_t old_value = p->value.second;
            uint32_t new_value = old_value + y;
            list_node<list<pair<uint32_t, uint32_t>>> *old_list = h1[old_value];
            list_node<list<pair<uint32_t, uint32_t>>> *new_list;
            if (h1.find(new_value) != h1.end())
                new_list = h1[new_value];
            else
            {
                list_node<list<pair<uint32_t, uint32_t>>> *tmp = new list_node<list<pair<uint32_t, uint32_t>>>;
                tmp->pre = old_list;
                while (tmp->pre->next->next != NULL && tmp->pre->next->value.get_head()->value.second < new_value)
                    tmp->pre = tmp->pre->next;
                tmp->next = tmp->pre->next;
                tmp->pre->next = tmp;
                tmp->next->pre = tmp;
                h1[new_value] = tmp;
                new_list = tmp;
            }

            old_list->value.del(p);
            //delete p;
            if (old_list->value.empty())
            {
                h1.erase(old_value);
                l.del(old_list);
                //delete old_list;
            }
            h2[x] = new_list->value.insert(make_pair(x, new_value), true);
        }
        else
        {
            if (op)
                ++failed_access;
            if (total_number < max_number)
            {
                ++total_number;
                if (h1.find(y) != h1.end())
                {
                    h2[x] = h1[y]->value.insert(make_pair(x, y), true);
                }
                else
                {
                    list_node<list<pair<uint32_t, uint32_t>>> *p = new list_node<list<pair<uint32_t, uint32_t>>>(list<pair<uint32_t, uint32_t>>(), true);
                    h2[x] = p->value.insert(make_pair(x, y), true);
                    h1[y] = p;
                    p->next = l.get_head();
                    while (p->next->next != NULL && p->next->value.get_head()->value.second < y)
                        p->next = p->next->next;
                    p->pre = p->next->pre;
                    p->pre->next = p;
                    p->next->pre = p;
                }
            }
            else
            {
                list_node<list<pair<uint32_t, uint32_t>>> *p = l.get_head();
                uint32_t min_value = p->value.get_head()->value.second;
                h2.erase(p->value.get_head()->value.first);
                list_node<pair<uint32_t, uint32_t>> *tmp = p->value.get_head();
                p->value.del(tmp);
                //delete tmp;
                if (p->value.empty())
                {
                    h1.erase(min_value);
                    l.del(p);
                    //delete p;
                }
                if (h1.find(y) != h1.end())
                {
                    h2[x] = h1[y]->value.insert(make_pair(x, y), true);
                }
                else
                {
                    p = new list_node<list<pair<uint32_t, uint32_t>>>(list<pair<uint32_t, uint32_t>>(), true);
                    h2[x] = p->value.insert(make_pair(x, y), true);
                    h1[y] = p;
                    p->next = l.get_head();
                    while (p->next->next != NULL && p->next->value.get_head()->value.second < y)
                        p->next = p->next->next;
                    p->pre = p->next->pre;
                    p->pre->next = p;
                    p->next->pre = p;
                }
            }
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