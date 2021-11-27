PeriodicSketch
============

Introduction
--------
Periodic is an important pattern in data streams, and periodic items refers to those items which arrive with a fixed interval. Any interval may occur many times, we use recurrence to denote the number of occurrence of an interval. Finding periodic items is a new and important issue in data streams mining, but there is no study carried out at present. To find periodic items, we propose a novel sketch, PeriodicSketch, aiming to accurately record top-𝐾 periodic items. To the best of our knowledge, this is the first work to find periodic items. To pick out periodic items with higher probability, we propose a key technique called Guaranteed Soft Uniform replacement strategy (GSU). Our theoretical proofs show that when replacement is successful, it is more likely that the new item has a higher recurrence than the current smallest recurrence; and GSU can ensure that our items in the sketch will approach the true periodic items closer and closer. And as soon as we get all the periodic items, the state would not change worse in high probability. We conduct extensive experiments on four real-world datasets. Experimental results show that the average AAE of our sketch using 1/10 memory is 737 times smaller than the baseline solution.

Repository structure
--------------------
*  `common/`: the hash function and bitmap data structure used by algorithms
*  `Algorithm/`: the implementation of our algorithm and baseline
*  `benchmark.h`: C++ header of some benchmarks about accuracy

Requirements
-------
- cmake
- g++

How to run
-------

```bash
$ cmake .
$ make
$ ./Periodic your-dataset
```
