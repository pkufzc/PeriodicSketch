#ifndef _list_H
#define _list_H

#include <cstdio>
using namespace std;

template <class T>
class list_node
{
public:
    bool v;
    T value;
    list_node *pre;
    list_node *next;
    list_node(T x = T(), bool y = false) : value(x), v(y)
    {
        pre = next = NULL;
    }
};

template <class T>
class list
{
private:
    list_node<T> *head;
    list_node<T> *tail;

public:
    list_node<T> *insert(T x = T(), bool y = false)
    {
        //puts("insert");
        list_node<T> *p = new list_node<T>(x, y);
        /*
        printf("%p %p %p\n", p, p->pre, p->next);
        printf("%p %p %p\n", head, head->pre, head->next);
        printf("%p %p %p\n", tail, tail->pre, tail->next);
        */
        p->pre = tail->pre;
        p->next = tail;
        p->pre->next = p;
        p->next->pre = p;
        //check();
        return p;
    }
    list(int len = 0)
    {
        head = new list_node<T>;
        tail = new list_node<T>;
        head->next = tail;
        tail->pre = head;
        for (int i = 0; i < len; ++i)
            insert();
        //check();
    }
    list(const list<T> &x)
    {
        head = new list_node<T>;
        tail = new list_node<T>;
        head->next = tail;
        tail->pre = head;
        for (list_node<T> *it = x.head->next; it != x.tail; it = it->next)
            insert(it->value, it->v);
        //check();
    }
    ~list()
    {
        for (list_node<T> *it = head; it != NULL;)
        {
            list_node<T> *tmp = it->next;
            delete it;
            it = tmp;
        }
    }
    void del(list_node<T> *p)
    {
        //puts("del");
        p->pre->next = p->next;
        p->next->pre = p->pre;
        delete p;
        //check();
    }
    void to_tail(list_node<T> *p)
    {
        p->pre->next = p->next;
        p->next->pre = p->pre;

        p->pre = tail->pre;
        p->next = tail;
        tail->pre->next = p;
        tail->pre = p;
    }
    list_node<T> *get_head()
    {
        return head->next;
    }
    bool empty()
    {
        return head->next == tail;
    }
    void check()
    {
        for (list_node<T> *it = head->next; it != tail; it = it->next)
            if (it->pre->next != it || it->next->pre != it)
            {
                puts("error1");
                exit(0);
            }
        if (head->next->pre != head)
        {
            puts("error2");
            exit(0);
        }
        if (tail->pre->next != tail)
        {
            puts("error3");
            exit(0);
        }
    }
    void operator=(const list<T> &x)
    {
        head = new list_node<T>;
        tail = new list_node<T>;
        head->next = tail;
        tail->pre = head;
        for (list_node<T> *it = x.head->next; it != x.tail; it = it->next)
            insert(it->value, it->v);
        //check();
    }
};

#endif