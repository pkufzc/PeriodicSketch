#ifndef BENCHMARK_H
#define BENCHMARK_H


#include <sys/stat.h>
#include <fstream>
#include <cstring>
#include <chrono>
#include <vector>
#include <sstream>
#include <unordered_map>  

#include "Algorithm/Baseline.h"
#include "Algorithm/Ours.h"
#include "Algorithm/tower.h"
#include "common/MMap.h"
static std::ofstream MP("mp.txt");

struct GroupConfig {
    uint32_t shared_memory;
    int levels;
    double tower_alpha;
    std::vector<double> ratios;
    std::vector<uint8_t> counters;
};

class BenchMark {
public:
    typedef std::vector<Abstract*> AbsVector;
    BenchMark(const std::string& _dataset_name) : dataset_name(_dataset_name){
    if (_dataset_name=="maccdc")
        PATH=" ";
    }

    ~BenchMark() {
        
    }

    void insert(const std::string& mode) {
         std::vector<ItemPair> pairs;
        pairs.reserve(20000000);       
           const std::string csv_path = PATH;
    std::ifstream csv_file(csv_path);
    if (!csv_file.is_open()) {
        std::cerr << "无法打开 CSV 文件: " << csv_path << std::endl;
        exit(1);
    }

    uint64_t total_lines = 0;    
    std::string line;
    std::unordered_map<DATA_TYPE, TIME_TYPE> last_time;  

    if (std::getline(csv_file, line)) {
        //std::cout << "跳过表头: " << line << std::endl;
    } else {
        std::cerr << "CSV 文件为空！" << std::endl;
        csv_file.close();
        exit(1);
    }

    while (std::getline(csv_file, line)) {
        total_lines++;
        std::istringstream iss(line);
        std::string timestamp_str, src_port_str, dst_port_str;

        
        if (!std::getline(iss, timestamp_str, ',') ||
            !std::getline(iss, src_port_str, ',') ||
            !std::getline(iss, dst_port_str, ',')) {
            std::cerr << "跳过格式错误的行（列数不足）: " << line << std::endl;
            continue;
        }

       
        double raw_timestamp;  
        uint32_t src_port;            
        uint32_t dst_port;            
        try {
            raw_timestamp = std::stod(timestamp_str);
            src_port = static_cast<uint32_t>(std::stoul(src_port_str));
            dst_port = static_cast<uint32_t>(std::stoul(dst_port_str));

        } catch (const std::exception& e) {
            std::cerr << "转换失败，跳过行: " << line << "，错误: " << e.what() << std::endl;
            continue;
        }

        TIME_TYPE processed_time = static_cast<uint64_t>(raw_timestamp * 1e5);

        DATA_TYPE item = (static_cast<DATA_TYPE>(src_port) << 32) | static_cast<DATA_TYPE>(dst_port);

        ItemPair pair{processed_time, item};
        pairs.push_back(pair);

    }
            HashMap est;
            const uint64_t total_records = total_lines;
        constexpr int32_t mem_base = 0;
        constexpr int32_t mem_inc = 200000;
        constexpr int32_t mem_var = 5;
        constexpr int32_t cmp_num = 3;
        constexpr int32_t round = 5;

        Abstract* sketches[mem_var][cmp_num];
        for (int i = 0; i < mem_var; ++i) {
            int memory = (mem_base + mem_inc * (i + 1)) / 1000;
            std::cout << "Memory size: " << memory << "KB" << std::endl << std::endl;

            double thp[round][cmp_num] = {};
            double time[round][cmp_num] = {};
            double avg_thp[cmp_num] = {};
            double avg_time[cmp_num] = {};

            for (int j = 0; j < round; ++j) {
                
                sketches[i][1] = new Ours<2>((i + 1) * mem_inc);
                sketches[i][0] = new Ours_tower((i + 1) * mem_inc,0.1,{0.1,0.35,0.2,0.35},{8,16,32,64});
                sketches[i][2] = new Baseline((i + 1) * mem_inc);
                
                for (int l = 0; l < cmp_num; ++l) {
                    if (sketches[i][l] == nullptr) continue;

                    if(mode=="q"){
                        
                        for (uint64_t k = 0; k < total_records; ++k) {
                            sketches[i][l]->Insert(pairs[k]);  
                        }  
                                        
                        TP start = now();  
                        HashMap est=sketches[i][l]->Report1(0.0001*total_records);  
                        TP end = now();  
                        double qsize=static_cast<double>(est.size());
                        
                        
                        double duration = std::chrono::duration<double>(end - start).count();
                        thp[j][l] = qsize / duration;
                        time[j][l]=duration;
                        }
                    else if(mode=="i")
                    {   
                        TP start = now(); 
                        for (uint64_t k = 0; k < total_records; ++k) {                           
                            sketches[i][l]->Insert(pairs[k]);    
                        }  
                        TP end = now();
                        double duration = std::chrono::duration<double>(end - start).count();
                        thp[j][l] = total_records / duration;
                        time[j][l]=duration;
                    }
                    avg_thp[l] += thp[j][l];
                    avg_time[l] += time[j][l];
                    
                    
                    if (j != round - 1) {
                        delete sketches[i][l];
                        sketches[i][l] = nullptr;
                    }
                }
            }

            
            for (int l = 0; l < cmp_num; ++l) {
                if (sketches[i][l] != nullptr) {             
                        std::cout<<sketches[i][l]->NAME<<"平均吞吐量："<< avg_thp[l] / round<<" records/sec "<<"运行时间： "<<avg_time[l]/round<<" sec"<<std::endl;
                    delete sketches[i][l];  
                    sketches[i][l] = nullptr;
                }
            }
            std::cout << std::endl;
        }

        std::cout << "共插入 " << total_records << " 条记录" << std::endl;
    }
    
private:
    uint64_t SIZE;
    double TOTAL;
    HashMap mp; 
    std::string PATH;
    std::string dataset_name;

};

#endif // BENCHMARK_TEST2_H
