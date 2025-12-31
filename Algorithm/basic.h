#ifndef BASIC_H
#define BASIC_H
#include<iostream>
#include <fstream>
#include "Abstract.h"
#include "../common/RandomUtil.h"

#define CELLNUM 8
template <uint32_t HASH_NUM>
class Ours : public Abstract{
public:

    struct Bucket{
        struct Cell{
            DATA_TYPE item;
            uint32_t interval;
            COUNT_TYPE count;
        };

        Cell cells[CELLNUM];
        COUNT_TYPE fail;
    };

    Ours(uint32_t _MEMORY, double ratio = 0.15, std::string _name = "Ours"){
        NAME = _name;
        MEMORY = _MEMORY;
        if (_MEMORY == 0) {
        std::cerr << "ERROR: LIGHT_LENGTH is 0! Use default value 100000." << std::endl;
        MEMORY = 100000;
    }
        LIGHT_LENGTH = MEMORY * ratio / HASH_NUM / sizeof(COUNT_TYPE);
        HEAVY_LENGTH = MEMORY * (1 - ratio) / sizeof(Bucket);
        // std::cout <<NAME<<' '<< "LIGHT_LENGTH: " << LIGHT_LENGTH << std::endl;
        // std::cout <<NAME<<' '<<"HEAVY_LENGTH: " << HEAVY_LENGTH << std::endl;
        for(uint32_t i = 0;i < HASH_NUM;++i){
            counter[i] = new TIME_TYPE[LIGHT_LENGTH];
            memset(counter[i], 0, sizeof(TIME_TYPE) * LIGHT_LENGTH);
        }

        buckets = new Bucket[HEAVY_LENGTH];
        memset(buckets, 0, sizeof(Bucket) * HEAVY_LENGTH);
    }

    ~Ours(){
        for(uint32_t i = 0;i < HASH_NUM;++i)
            delete [] counter[i];
        delete [] buckets;
    }

    void Insert(const ItemPair& item){
        TIME_TYPE minTime = UINT64_MAX;

        for(uint32_t i = 0;i < HASH_NUM;++i){
            uint32_t position = ::hash(item.item, i) % LIGHT_LENGTH;
            minTime = MIN(minTime, counter[i][position]);
            counter[i][position] = item.time;
        }

        ItemPair temp(item.time - minTime, item.item);
        if(temp.time >= 1250)
            return;
        uint32_t position = ::hash(temp, HASH_NUM) % HEAVY_LENGTH, min_pos = 0;
        COUNT_TYPE min_count = INT32_MAX;

        for(uint32_t j = 0;j < CELLNUM;++j){
            if(buckets[position].cells[j].item == temp.item &&
                    buckets[position].cells[j].interval == temp.time){
                buckets[position].cells[j].count += 1;
                return;
            }

            if(buckets[position].cells[j].count < min_count){
                min_count = buckets[position].cells[j].count;
                min_pos = j;
            }
        }
        std::uniform_int_distribution<> dist(0, (2 * min_count - buckets[position].fail));
        if (dist(getRNG()) ==0) {
            buckets[position].cells[min_pos].item = temp.item;
            buckets[position].cells[min_pos].interval = temp.time;
            if(min_count == 0)
                buckets[position].cells[min_pos].count = 1;
            else
                buckets[position].cells[min_pos].count += buckets[position].fail / min_count;
            buckets[position].fail = 0;
        }
        else{
            buckets[position].fail += 1;
        }
    }

    HashMap Report1(COUNT_TYPE HIT){

        HashMap ret;
        for(uint32_t i = 0;i < HEAVY_LENGTH;++i){
            for(uint32_t j = 0;j < CELLNUM;++j){
                if(buckets[i].cells[j].count > HIT){
                    ret[ItemPair(buckets[i].cells[j].interval,
                            buckets[i].cells[j].item)]
                    = buckets[i].cells[j].count;
                }
            }
        }

        return ret;
    }

private:
    uint32_t LIGHT_LENGTH;
    uint32_t HEAVY_LENGTH;
    inline std::mt19937& getRNG() {
        thread_local std::mt19937 rng = getLocalRNG(); // 由 RandomUtil.h 中的函数初始化
        return rng;
    }
    TIME_TYPE* counter[HASH_NUM];
    Bucket* buckets;

};

#endif//OURS_TEST_H
