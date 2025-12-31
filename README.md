PeriodicSketch
============

Introduction
--------
Periodicity is a crucial pattern in data streams, with periodic items referring to those arriving at fixed time intervals. "Recurrence" denotes the occurrence count of a specific interval, and mining such items is a novel yet under-explored issue in data stream mining.  

To address this, we first propose PeriodicSketch, the first sketch dedicated to mining top-K high-frequency periodic items, leveraging the Guaranteed Soft Uniform (GSU) replacement strategy to enhance the selection of high-value items.  

To overcome PeriodicSketch’s key limitation of static memory allocation, we present its advanced version PeriodicSketch+. Retaining the core "interval estimation → periodic element storage" architecture, it achieves performance leaps via three collaborative optimizations: 1) replacing the original Cover-Min sketch with customized TowerSketch to boost interval estimation accuracy and reduce errors; 2) upgrading to an elastic GSU sketch for dynamic memory management; 3) adding a collaborative memory reallocation mechanism to dynamically adjust resource distribution between modules under a fixed total memory budget, enabling robust adaptation to dynamic data streams. Experiments show PeriodicSketch+ outperforms PeriodicSketch and Baseline in precision, recall, and error metrics with stronger scene adaptability under the same memory budget.
 If you have any questions about PeriodicSketch, please contact the following email: fanzc@pku.edu.cn or yangtongemail@gmail.com

Repository structure
--------------------
*  `common/`: the hash function and bitmap data structure used by algorithms
*  `Algorithm/`: the implementation of the baseline solution, the basic version PeriodicSketch (basic.h), and the advanced version PeriodicSketch+ (adv.h).
*  `benchmark.h`: C++ header of some benchmarks about accuracy

Requirements
-------
- cmake
- g++

How to run
-------

```bash
$ git clone https://github.com/pkufzc/PeriodicSketch.git
$ cd ./PeriodicSketch/PeriodicSketch-main/
$ cmake .
$ make
$ ./Periodic {your-dataset}
```
