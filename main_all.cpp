#include "Benchmark.h"  
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <windows.h>
#include <io.h>
#include <fcntl.h>

int main(int argc, char *argv[]) {      
    std::string set_name;
    std::vector<GroupConfig> all_groups(1);

    set_name=argv[1];
    all_groups[0].shared_memory = std::stoi(argv[2]);
    //  all_groups[1].shared_memory = std::stoi(argv[2]);
    // all_groups[2].shared_memory = std::stoi(argv[2]);
    // all_groups[3].shared_memory = std::stoi(argv[2]);
    // all_groups[4].shared_memory = std::stoi(argv[2]);

    // all_groups[0].levels = 1;
    // all_groups[0].ratios = {1};
    // all_groups[0].counters = {32};
    // all_groups[0].tower_alpha = 0.1;

    // all_groups[0].levels = 2;
    // all_groups[0].ratios = {0.5,0.5};
    // all_groups[0].counters = {32,64};
    // all_groups[0].tower_alpha = 0.1;

    // all_groups[0].levels = 3;
    // all_groups[0].ratios = {0.33, 0.33, 0.33};
    // all_groups[0].counters = {16,32,64};
    // all_groups[0].tower_alpha = 0.1;

    all_groups[0].levels = 4;
    all_groups[0].ratios = {0.25,0.25,0.25,0.25};
    all_groups[0].counters = {8,16,32,64};
    all_groups[0].tower_alpha = 0.1;

    // all_groups[0].levels = 5;
    // all_groups[0].ratios = {0.2,0.2,0.2,0.2,0.2};
    // all_groups[0].counters = {8,16,32,64,64};
    // all_groups[0].tower_alpha = 0.1;    

    for(int i=0;i<5;i++){
    BenchMark bench(set_name.c_str());
    bench.SingleRunAllGroups(
        0.0001,
        2,
        all_groups  
    );
    }
    return 0;
}