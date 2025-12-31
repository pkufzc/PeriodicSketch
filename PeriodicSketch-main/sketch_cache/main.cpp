#include <cstdio>
#include <cstring>
#include <algorithm>
#include <ctime>
#include <winsock2.h>
#include <ws2tcpip.h>
#pragma comment(lib, "ws2_32.lib")
#include <vector>
#include "../common/Util.h"
#include "LRU.h"
#include "LFU.h"
#include "ours_LFU.h"
#include "ours_LRU.h"
#include "tower_LRU.h"
#include "tower_LFU.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <queue>
using namespace std;

struct ktp
{
    uint64_t timestamp;
    uint32_t key;
};

vector<ktp> ktpv;

bool comp(ktp x, ktp y)
{
    return x.timestamp < y.timestamp;
}

void work(uint64_t T, uint64_t advance, int ours_length)
{
    int lens[] = {10000};
    int num_lens = sizeof(lens) / sizeof(int);
    for (int i = 0; i < num_lens; ++i)
    {
       tower_lru tl(lens[i], T, advance, ours_length);
        ours_lru ol(lens[i], T, advance, ours_length);
        lru l(lens[i]);
        tower_lfu tlf(lens[i], T, advance, ours_length);
        ours_lfu olf(lens[i], T, advance, ours_length);
        lfu lf(lens[i]);
        int count = 0;
        uint64_t last_timestamp = 0;
        for (vector<ktp>::iterator it = ktpv.begin(); it != ktpv.end(); ++it)
        {
            ++count;
            //printf("%d\n", count);
            if (count == 1000000)
            {
               tl.clear_access();
               ol.clear_access();
               l.clear_access();
               tlf.clear_access();
                olf.clear_access();
               lf.clear_access();
            }
            //printf("%u %lu\n", it->key, it->timestamp);
            last_timestamp = it->timestamp;
            tl.access(it->key, it->timestamp);
           ol.access(it->key, it->timestamp);
           //l.access(it->key);
            tlf.access(it->key, it->timestamp);
            olf.access(it->key, it->timestamp);
          // lf.access(it->key);
        }
        printf("%d,%.9lf,%.9lf,%.9lf,%.9lf,%.9lf,%.9lf,%d,%d\n", lens[i], lf.get_hit_rate(), olf.get_hit_rate(),tlf.get_hit_rate(), l.get_hit_rate(), ol.get_hit_rate(),tl.get_hit_rate(), olf.periodic_count, count);
    }
}

void caida()
{
    int lens[] = {30, 40, 50, 70, 100, 150, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000, 100000};
    int num_lens = sizeof(lens) / sizeof(int);
    
    FILE *file = fopen("", "rb");
    if (file == NULL) {
        perror("fopen failed");
        return; 
    }
    
    for (int i = 0; i < num_lens; ++i)
    {   
        tower_lru tl(lens[i], 100000000, 100, 300000);
        ours_lru ol(lens[i], 100000000, 100, 300000);
        lru l(lens[i]);
        tower_lfu tlf(lens[i], 100000000, 100, 300000);
        ours_lfu olf(lens[i], 100000000, 100, 300000);
        lfu lf(lens[i]);
        
        char buf[21];
        uint32_t last_key = 0;
        uint64_t last_timestamp = 0;
        int count = 0;

        rewind(file);

        while (1)
        {
            if (!fread(buf, 21, 1, file)){
                if (feof(file)) {
                    // printf("已读取完所有数据\n");
                } else {
                    perror("fread failed");
                }
                break;
            }

            ++count;
            if (count == 1000000)
            {
                l.clear_access();
                ol.clear_access();
                tl.clear_access();
                lf.clear_access();
                olf.clear_access();
                tlf.clear_access();
            }

            uint32_t key = Hash::BOBHash32((const uint8_t *)buf, 13, 1000);
            double timestamp;
            memcpy(&timestamp, buf + 13, sizeof(double));
            tl.access(key, uint64_t(timestamp * 100000000));
            ol.access(key, uint64_t(timestamp * 100000000));
            l.access(key);
            tlf.access(key, uint64_t(timestamp * 100000000));
            olf.access(key, uint64_t(timestamp * 100000000));
            lf.access(key);
        }

        printf("%d,%.9lf,%.9lf,%.9lf,%.9lf,%.9lf,%.9lf,%d,%d\n", 
               lens[i], 
               lf.get_hit_rate(), 
               olf.get_hit_rate(),
               tlf.get_hit_rate(), 
               l.get_hit_rate(), 
               ol.get_hit_rate(), 
               tl.get_hit_rate(),
               olf.periodic_count, 
               count);
    }

    fclose(file);
}

void criteo()
{
    int lens[] = {30, 40, 50, 70, 100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000, 100000};
    int num_lens = sizeof(lens) / sizeof(int);
    for (int i = 0; i < num_lens; ++i)
    {
        lru l(lens[i]);
        ours_lru ol(lens[i], 50000, 2, 300000);
        tower_lru tl(lens[i], 50000, 2, 300000);
        lfu lf(lens[i]);
        ours_lfu olf(lens[i], 50000, 2, 300000);
        tower_lfu tlf(lens[i], 50000, 2, 300000);
        FILE *file = fopen("", "r");
        char buf[1000];
        fscanf(file, "%[^\n]", buf);
        uint32_t last_key = 0;
        uint64_t last_timestamp = 0, id;
        int count = 0;
        while (1)
        {
            if (fscanf(file, "%lu%lu%lu", &last_timestamp, &id, &id) == EOF)
                break;

            fscanf(file, "%[^\n]", buf);
            ++count;
            if (count == 1000000)
            {
                l.clear_access();
                ol.clear_access();
                tl.clear_access();
                lf.clear_access();
                olf.clear_access();
                tl.clear_access();
            }
            if (count == 10000000)
                break;
            uint32_t key = Hash::BOBHash32((const uint8_t *)&id, 8, 1000);
            l.access(key);
            ol.access(key, last_timestamp);
            tl.access(key, last_timestamp);
            lf.access(key);
            olf.access(key, last_timestamp);
            tlf.access(key, last_timestamp);
        }
        printf("%d,%.9lf,%.9lf,%.9lf,%.9lf,%.9lf,%.9lf,%d,%d\n", lens[i], lf.get_hit_rate(), olf.get_hit_rate(), tlf.get_hit_rate(),l.get_hit_rate(), ol.get_hit_rate(),tl.get_hit_rate(), olf.periodic_count, count);
        fclose(file);
    }
}

void random_data()
{
     for (int i = 5000; i <= 50000; i += 5000){
       
        uint64_t count = 0;
        ktpv.clear();
        
        for (int j = 1; j <= 2000; ++j)
        {

            for (int k = 1; k <= i; ++k){
                ktp temp;
                temp.timestamp = ++count;
                temp.key = k;
                ktpv.push_back(temp);
            }

            for (int k = i + 1; k <= 50000; ++k){
                ktp temp1;
                temp1.timestamp = ++count;
                
                temp1.key =  (((uint32_t)rand() << 15) | rand());  
                       
                ktpv.push_back(temp1);
            }
        }
        
        work(500000, 1, 1800000);
}
}



int main(int argc, char **argv)
{
    srand(time(0));
    if (argc != 2)
    {
        puts("argc must be 2");
        return 0;
    }
    if (strcmp(argv[1], "caida") == 0)
        caida();
    else if (strcmp(argv[1], "criteo") == 0)
        criteo();

    else if (strcmp(argv[1], "random") == 0)
        random_data();

    return 0;
}