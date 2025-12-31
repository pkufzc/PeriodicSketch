// RandomUtil.h 修正版（核心部分）
#ifndef RANDOM_UTIL_H
#define RANDOM_UTIL_H

#include <random>
#include <thread>
#include <cstdint>

// 1. 全局统一种子
#include <chrono>  // 必须包含：用于获取高精度时间戳

// 全局种子：结合时间戳+硬件随机数（优先）+线程ID，确保每次运行不同
static const uint64_t GLOBAL_SEED = []() {
    // 1. 获取高精度时间戳（微秒级，避免秒级重复）
    auto now = std::chrono::high_resolution_clock::now().time_since_epoch();
    uint64_t time_seed = std::chrono::duration_cast<std::chrono::microseconds>(now).count();

    // 2. 尝试获取硬件随机数（增强随机性）
    std::random_device rd;
    uint64_t hardware_seed = 0;
    if (rd.entropy() > 0) {  // 检测是否支持硬件随机（如Linux /dev/urandom）
        hardware_seed = static_cast<uint64_t>(rd()) << 32 | static_cast<uint64_t>(rd());  // 64位硬件种子
    } else {
        hardware_seed = 20041003;  // 备用固定值（仅当无硬件随机时使用）
    }

    // 3. 混合时间戳和硬件种子（异或运算，确保双重随机性）
    return time_seed ^ hardware_seed;
}();

// 线程局部随机数生成器：每个线程独立初始化，避免并发冲突
inline std::mt19937& getLocalRNG() {
    thread_local std::mt19937 rng = [](){
        // 线程种子 = 全局种子 ^ 线程ID哈希（确保每个线程序列不同）
        uint64_t thread_seed = GLOBAL_SEED ^ std::hash<std::thread::id>()(std::this_thread::get_id());
        return std::mt19937(static_cast<uint32_t>(thread_seed));  // mt19937接受32位种子
    }();
    return rng;
}
#endif // RANDOM_UTIL_H