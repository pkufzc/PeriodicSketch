#ifndef UTIL_H
#define UTIL_H

#include <unordered_map>

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

bool operator == (const ItemPair& a, const ItemPair& b){
    return memcmp(&a, &b, sizeof(ItemPair)) == 0;
}

namespace std{
    template<>
    struct hash<ItemPair>{
        size_t operator()(const ItemPair& item) const noexcept
        {
            return Hash::BOBHash32((uint8_t*)&item, sizeof(ItemPair), 0);
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
