#ifndef BENCHMARK_H
#define BENCHMARK_H
#include <sys/stat.h>
#include <fstream>
#include <cstring>
#include <chrono>
#include <vector>
#include <sstream>
#include <unordered_map>  
#include <iomanip>

#include "Algorithm/Abstract.h"
#include "Algorithm/basic.h"
#include "Algorithm/Baseline.h"
#include "Algorithm/adv.h"
#include "common/MMap.h"


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
        PATH="C:\\Users\\ADOL\\Desktop\\maccdc1935.csv";
        result = Load(PATH.c_str());
        uint8_t* data = static_cast<uint8_t*>(result.start);
        uint64_t total_len = result.length;
        const uint64_t record_size = 21;
        SIZE = total_len / record_size;
        TOTAL = static_cast<double>(SIZE);
    }

    ~BenchMark() {
        UnLoad(result);  
    }

    void SingleRunAllGroups(double alpha,
                           int num_algorithms,
                           const std::vector<GroupConfig>& all_groups) {
        
        AbsVector all_sketches;  
        Abstract* ours_baseline = nullptr;  
        Abstract* baseline = nullptr;
        if(all_groups[0].shared_memory!=0)
          {  
            ours_baseline = new Ours<2>(all_groups[0].shared_memory, 0.15, "Ours");
            baseline = new Baseline((all_groups[0].shared_memory)*10);
            all_sketches.push_back(ours_baseline);  
            all_sketches.push_back(baseline);
        
            for (size_t g = 0; g < all_groups.size(); ++g) {
                const auto& group = all_groups[g];

                std::string group_name =  "our_tower_group" + std::to_string(g + 1);
                          
                Abstract* tower_sketch = new Ours_tower(
                    group.shared_memory,                
                    group.tower_alpha,            
                    group.ratios,                 
                    group.counters,               
                    group_name.c_str()            
                );
                all_sketches.push_back(tower_sketch);}
            }
    

            insert(all_sketches);


        uint32_t threshold = alpha * SIZE;

        if (baseline != nullptr) {
           CheckError(baseline, threshold, 0);  
        }
        if (ours_baseline != nullptr) {
            CheckError(ours_baseline, threshold, 0);  
        }
 
        size_t tower_start_idx = (ours_baseline != nullptr) ? 2 : 0;
        for (size_t g = 0; g < all_groups.size(); ++g) {
            size_t sketch_idx = tower_start_idx + g ;
            if (sketch_idx >= all_sketches.size()) break;
            Abstract* tower_sketch = all_sketches[sketch_idx];
            CheckError(tower_sketch, threshold, g + 1);  // 传递组号标识
        }


        for (auto sketch : all_sketches) {
            delete sketch;
        }
        all_sketches.clear();
    
    }

private:
    LoadResult result;
    uint64_t SIZE;
    double TOTAL;
    HashMap mp; 
    std::string PATH;
    std::string dataset_name;

void insert(AbsVector& sketches) {
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

        TIME_TYPE interval = processed_time - last_time[item];

        mp[{interval, item}]++;

        last_time[item] = processed_time;

        ItemPair pair{processed_time, item};
        for (auto sketch : sketches) {
            if (sketch != nullptr) {
                sketch->Insert(pair);
            }
        }
    }

    csv_file.close();
    SIZE = total_lines;
    std::cout << "数据插入完成！共插入 " << SIZE << " 条有效记录（已跳过表头）" << std::endl;
}

void CheckError(Abstract* sketch, uint32_t threshold, int group_id) {
    if (!sketch) {
        std::cerr << "ERROR: sketch is null in CheckError (group " << group_id << ")" << std::endl;
        return;
    }

    HashMap estimated1 = sketch->Report1(threshold);
    double heavy_aae = 0.0;
    double heavy_are = 0.0;
    double heavy_both = 0.0;
    double heavy_est_total = static_cast<double>(estimated1.size());
    double heavy_real_total = 0.0;
    double hit_count=0;
    double hit_real=0;

    for (const auto& [key, real_count] : mp) {
        if (real_count > threshold) {
            std::ofstream MP("mp", std::ios::app);
            MP << '(' << key.time << ',' << key.item << ',' << real_count << ')';
            MP.flush();
            hit_real+=real_count;

            heavy_real_total += 1.0;

            if (estimated1.find(key) != estimated1.end()) {
                
                heavy_both += 1.0;
                COUNT_TYPE est_count1 = estimated1[key];
                hit_count+=est_count1;
                heavy_aae += std::abs(static_cast<int64_t>(est_count1 - real_count));
                heavy_are += (real_count > 0) ? std::abs(static_cast<double>(est_count1 - real_count)) / real_count : 0.0;
            }
        }
    }

    double heavy_recall = (heavy_real_total > 0) ? (heavy_both / heavy_real_total) : 0.0;
    double heavy_precision = (heavy_est_total > 0) ? (heavy_both / heavy_est_total) : 0.0;
    double heavy_f1 = (heavy_precision + heavy_recall > 0) ? 
                      (2.0 * heavy_precision * heavy_recall) / (heavy_precision + heavy_recall) : 0.0;
    
    heavy_aae = (heavy_both > 0) ? (heavy_aae / heavy_both) : 0.0;
    heavy_are = (heavy_both > 0) ? (heavy_are / heavy_both) : 0.0;

    std::cout << "算法名称: " << sketch->NAME << std::endl;
    std::cout << "内存配置: " << sketch->MEMORY << " bytes" << std::endl;

    std::cout << "AAE: " << heavy_aae << std::endl;
    std::cout << "ARE: " << heavy_are << std::endl;
    std::cout << "Recall: " << heavy_recall << std::endl;
    std::cout << "Precision: " << heavy_precision << std::endl;
    std::cout << "F1-Score: " << heavy_f1 << std::endl;
    std::cout << hit_real<< std::endl;
    std::cout << hit_count<< "个" << std::endl;



    std::ofstream outfile("performance_results.txt", std::ios::app);

    if (outfile.is_open()) {        
        outfile << sketch->NAME << "\t"
                << sketch->MEMORY << "\t"
                << std::fixed << std::setprecision(6)
                << heavy_aae << "\t"
                << heavy_are << "\t"
                << heavy_recall << "\t"
                << heavy_precision << "\t"
                << heavy_f1 << "\n";
        outfile.close();
    } else {
        std::cerr << "ERROR: 无法打开文件 performance_results.txt" << std::endl;
    }
}

void InitPerformanceFile() {
    std::ofstream outfile("performance_results.txt", std::ios::trunc);
    if (outfile.is_open()) {
        outfile << "GroupID\tAlgorithm\tMemory(bytes)\tAAE\tARE\tRecall\tPrecision\tF1-Score\n";
        outfile.close();
    } else {
        std::cerr << "ERROR: 无法创建文件 performance_results.txt" << std::endl;
    }
}
};

#endif // BENCHMARK_H
