#ifndef UTIL_H
#define UTIL_H
#include<cstdint>
#include <unordered_map>
#include <chrono>
#include "hash.h"
#include <functional>

#define DELTA 1000

typedef uint64_t DATA_TYPE;
typedef uint64_t TIME_TYPE;
typedef int32_t COUNT_TYPE;

struct ItemPair{
    TIME_TYPE time;
    DATA_TYPE item;

    ItemPair(TIME_TYPE _time = 0, DATA_TYPE _item = 0):
        time(_time), item(_item){}
};

bool operator==(const ItemPair& a, const ItemPair& b) {
    return a.time == b.time && a.item == b.item;
}
namespace std {
    template<>
    struct hash<ItemPair> {
        size_t operator()(const ItemPair& p) const noexcept {
            auto h1 = hash<TIME_TYPE>{}(p.time);
            auto h2 = hash<DATA_TYPE>{}(p.item);
            // 使用 boost 风格的组合（广泛推荐）
            return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
        }
    };
}
typedef std::unordered_map<ItemPair, COUNT_TYPE> HashMap;

typedef std::chrono::high_resolution_clock::time_point TP;

inline TP now(){
    return std::chrono::high_resolution_clock::now();
}

inline double durationms(TP finish, TP start){
    return std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1,1000000>>>(finish - start).count();
}

#endif
