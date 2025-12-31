cmake_minimum_required(VERSION 3.14)
project(Periodic)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)
include_directories(Algorithm)
include_directories(common)

add_executable(Periodic
        Algorithm/Abstract.h
        Algorithm/Baseline.h
        Algorithm/basic.h
        common/BitMap.h
        common/CuckooMap.h
        common/hash.h
        common/Heap.h
        common/MMap.h
        common/StreamSummary.h
        common/Util.h
        benchmark.h
        main_all.cpp)
