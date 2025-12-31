#include "Benchmark_thp.h"  
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <windows.h>
#include <io.h>
#include <fcntl.h>

int main(int argc, char *argv[]) {      
    std::string set_name;
    std::string mode;
    std::vector<GroupConfig> all_groups(1);

    set_name=argv[1];
    mode = argv[2];
    all_groups[0].levels = 4;
    all_groups[0].ratios = {0.25,0.25,0.25,0.25};
    all_groups[0].counters = {8,16,32,64};
    all_groups[0].tower_alpha = 0.1;
    BenchMark bench(set_name.c_str());
    bench.insert(mode);

    return 0;
}