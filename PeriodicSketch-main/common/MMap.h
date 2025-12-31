// mmap.h
#ifndef MMAP_H
#define MMAP_H

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <vector>

struct LoadResult {
    void* start;
    uint64_t length;
    std::vector<uint8_t> data;
};

LoadResult Load(const char* PATH) {
    LoadResult ret;
    std::ifstream file(PATH, std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file");
    }

    ret.length = file.tellg();
    file.seekg(0, std::ios::beg);

    ret.data.resize(ret.length);
    file.read(reinterpret_cast<char*>(ret.data.data()), ret.length);
    ret.start = ret.data.data();

    if (!file) {
        throw std::runtime_error("Failed to read file");
    }

    return ret;
}

void UnLoad(LoadResult result) {
    // vector 自动释放
}

#endif