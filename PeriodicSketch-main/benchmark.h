#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <sys/stat.h>
#include <fstream>

#include <chrono>
#include <vector>

#include "Baseline.h"
#include "Ours.h"

#include "MMap.h"

class BenchMark{
public:

    typedef std::vector<Abstract*> AbsVector;

    BenchMark(const char* _PATH):
            PATH(_PATH){
        result = Load(PATH);

        start = (ItemPair*)result.start;
        SIZE = result.length / sizeof(ItemPair);
        TOTAL = SIZE;

        std::unordered_map<DATA_TYPE, TIME_TYPE> last_time;

        for(uint32_t i = 0;i < SIZE;++i){
            start[i].time /= DELTA;
            if(last_time.find(start[i].item) != last_time.end()){
                TIME_TYPE interval = start[i].time - last_time[start[i].item];
                mp[ItemPair(interval, start[i].item)] += 1;
            }
            last_time[start[i].item] = start[i].time;
        }
    }

    ~BenchMark(){
        UnLoad(result);
    }

    void TopKError(double alpha){
        AbsVector FPIs = {
                new Ours<2>(100000),
                new Baseline(100000),
        };

        BenchInsert(FPIs);

        for(auto FPI : FPIs){
            std::cout << "Check," << FPI->NAME << "," << FPI->MEMORY << std::endl;
            CheckError(FPI, alpha * TOTAL);
            delete FPI;
        }
    }

private:
    LoadResult result;
    ItemPair* start;

    uint64_t SIZE;
    double TOTAL;

    HashMap mp;
    const char* PATH;

    void BenchInsert(AbsVector& sketches){
        for(uint32_t i = 0;i < SIZE;++i){
            for(auto sketch : sketches)
                sketch->Insert(start[i]);
        }
    }

    void InsertThp(Abstract* sketch){
        TP initial, finish;

        initial = now();
        for(uint32_t i = 0;i < SIZE;++i){
            sketch->Insert(start[i]);
        }
        finish = now();

	    std::cout << "Insertion Thp: " << SIZE / durationms(finish, initial) << std::endl;
    }


    void CheckError(Abstract* sketch, uint32_t HIT){
        double real = 0, estimate = 0, both = 0;
        double aae = 0, are = 0, cr = 0, pr = 0, f1 = 0;

        double hot = 0;

        HashMap temp = sketch->Report(HIT);

        estimate = temp.size();

        for(auto it = mp.begin();it != mp.end();++it){
            if(it->second > HIT){
                real += 1;
                if(temp.find(it->first) != temp.end()){
                    both += 1;
                    COUNT_TYPE value = temp[it->first];
                    aae += abs(it->second - value);
                    are += abs(it->second - value) / (double)it->second;
                }
            }
        }

        if(both <= 0){
            std::cout << "Not Real" << std::endl;
        }
        else{
            aae /= both;
            are /= both;
        }

        cr = both / real;

        if(estimate <= 0){
            std::cout << "Not Find" << std::endl;
        }
        else{
            pr = both / estimate;
        }

	    std::cout << "AAE," << aae << std::endl
		    << "ARE," << are << std::endl
		    << "Recall-Rate," << cr << std::endl
		    << "Precision-Rate," << pr << std::endl;
    }

};

#endif
