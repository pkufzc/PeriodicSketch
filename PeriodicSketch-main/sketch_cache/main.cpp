#include <cstdio>
#include <cstring>
#include <algorithm>
#include <ctime>
#include <arpa/inet.h>
#include <vector>
#include "LRU.h"
#include "ours_LRU.h"
#include "LFU.h"
#include "ours_LFU.h"
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
    int lens[] = {30, 40, 50, 70, 100, 150, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000, 100000};
    int num_lens = sizeof(lens) / sizeof(int);
    for (int i = 0; i < num_lens; ++i)
    {
        lru l(lens[i]);
        ours_lru ol(lens[i], T, advance, ours_length);
        lfu lf(lens[i]);
        ours_lfu olf(lens[i], T, advance, ours_length);
        int count = 0;
        uint64_t last_timestamp = 0;
        for (vector<ktp>::iterator it = ktpv.begin(); it != ktpv.end(); ++it)
        {
            ++count;
            //printf("%d\n", count);
            if (count == 1000000)
            {
                l.clear_access();
                ol.clear_access();
                lf.clear_access();
                olf.clear_access();
            }
            //printf("%u %lu\n", it->key, it->timestamp);
            last_timestamp = it->timestamp;
            l.access(it->key);
            ol.access(it->key, it->timestamp);
            lf.access(it->key);
            olf.access(it->key, it->timestamp);
        }
        printf("%d,%.9lf,%.9lf,%.9lf,%.9lf,%d,%d\n", lens[i], lf.get_hit_rate(), olf.get_hit_rate(), l.get_hit_rate(), ol.get_hit_rate(), olf.periodic_count, count);
    }
}

void mawi()
{
    int lens[] = {30, 40, 50, 70, 100, 150, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000, 100000};
    int num_lens = sizeof(lens) / sizeof(int);
    for (int i = 0; i < num_lens; ++i)
    {
        lru l(lens[i]);
        ours_lru ol(lens[i], 20000, 2, 300000);
        lfu lf(lens[i]);
        ours_lfu olf(lens[i], 20000, 2, 300000);
        FILE *file = fopen("/share/pcap_zhangyd/final.txt", "r");
        char buf[21], sip_s[20], dip_s[20];
        uint32_t sip, dip;
        uint16_t sport, dport;
        uint8_t proto;
        uint32_t last_key = 0;
        uint64_t last_timestamp = 0;
        int count = 0;
        while (1)
        {
            if (fscanf(file, "%s%s%hu%hu%hhu", sip_s, dip_s, &sport, &dport, &proto) == EOF)
                break;

            sip = inet_addr(sip_s);
            dip = inet_addr(dip_s);

            memcpy(buf, &sip, 4);
            memcpy(buf + 4, &dip, 4);
            memcpy(buf + 8, &sport, 2);
            memcpy(buf + 10, &dport, 2);
            memcpy(buf + 12, &proto, 1);

            ++count;
            if (count == 1000000)
            {
                l.clear_access();
                ol.clear_access();
                lf.clear_access();
                olf.clear_access();
            }
            uint32_t key = Hash::BOBHash32((const uint8_t *)buf, 13, 1000);
            last_timestamp++;
            l.access(key);
            ol.access(key, last_timestamp);
            lf.access(key);
            olf.access(key, last_timestamp);
        }
        printf("%d,%.9lf,%.9lf,%.9lf,%.9lf,%d,%d\n", lens[i], lf.get_hit_rate(), olf.get_hit_rate(), l.get_hit_rate(), ol.get_hit_rate(), olf.periodic_count, count);
        fclose(file);
    }
}

void caida()
{
    int lens[] = {30, 40, 50, 70, 100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000, 100000};
    int num_lens = sizeof(lens) / sizeof(int);
    for (int i = 0; i < num_lens; ++i)
    {
        lru l(lens[i]);
        ours_lru ol(lens[i], 100000000, 100, 300000);
        lfu lf(lens[i]);
        ours_lfu olf(lens[i], 100000000, 100, 300000);
        FILE *file = fopen("/root/dataset/130000.dat", "r");
        char buf[21];
        uint32_t last_key = 0;
        uint64_t last_timestamp = 0;
        int count = 0;
        while (1)
        {
            if (!fread(buf, 21, 1, file))
                break;

            ++count;
            if (count == 1000000)
            {
                l.clear_access();
                ol.clear_access();
                lf.clear_access();
                olf.clear_access();
            }
            uint32_t key = Hash::BOBHash32((const uint8_t *)buf, 13, 1000);

            double timestamp;
            memcpy(&timestamp, buf + 13, sizeof(double));
            l.access(key);
            ol.access(key, uint64_t(timestamp * 100000000));
            lf.access(key);
            olf.access(key, uint64_t(timestamp * 100000000));
        }
        printf("%d,%.9lf,%.9lf,%.9lf,%.9lf,%d,%d\n", lens[i], lf.get_hit_rate(), olf.get_hit_rate(), l.get_hit_rate(), ol.get_hit_rate(), olf.periodic_count, count);
        fclose(file);
    }
}

void webdocs()
{
    int lens[] = {30, 40, 50, 70, 100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000, 100000};
    int num_lens = sizeof(lens) / sizeof(int);
    for (int i = 0; i < num_lens; ++i)
    {
        lru l(lens[i]);
        ours_lru ol(lens[i], 20000, 2, 300000);
        lfu lf(lens[i]);
        ours_lfu olf(lens[i], 20000, 2, 300000);
        FILE *file = fopen("/root/webdocs.dat", "r");
        char buf[21];
        uint32_t last_key = 0;
        uint64_t last_timestamp = 0, id;
        int count = 0;
        while (1)
        {
            if (fscanf(file, "%lu", &id) == EOF)
                break;

            ++count;
            if (count == 1000000)
            {
                l.clear_access();
                ol.clear_access();
                lf.clear_access();
                olf.clear_access();
            }
            if (count == 10000000)
                break;
            uint32_t key = Hash::BOBHash32((const uint8_t *)&id, 8, 1000);
            ++last_timestamp;
            l.access(key);
            ol.access(key, last_timestamp);
            lf.access(key);
            olf.access(key, last_timestamp);
        }
        printf("%d,%.9lf,%.9lf,%.9lf,%.9lf,%d,%d\n", lens[i], lf.get_hit_rate(), olf.get_hit_rate(), l.get_hit_rate(), ol.get_hit_rate(), olf.periodic_count, count);
        fclose(file);
    }
}

void criteo()
{
    int lens[] = {30, 40, 50, 70, 100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000, 100000};
    int num_lens = sizeof(lens) / sizeof(int);
    for (int i = 0; i < num_lens; ++i)
    {
        lru l(lens[i]);
        ours_lru ol(lens[i], 50000, 2, 300000);
        lfu lf(lens[i]);
        ours_lfu olf(lens[i], 50000, 2, 300000);
        FILE *file = fopen("/root/dataset_criteo/criteo_attribution_dataset.tsv", "r");
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
                lf.clear_access();
                olf.clear_access();
            }
            if (count == 10000000)
                break;
            uint32_t key = Hash::BOBHash32((const uint8_t *)&id, 8, 1000);
            l.access(key);
            ol.access(key, last_timestamp);
            lf.access(key);
            olf.access(key, last_timestamp);
        }
        printf("%d,%.9lf,%.9lf,%.9lf,%.9lf,%d,%d\n", lens[i], lf.get_hit_rate(), olf.get_hit_rate(), l.get_hit_rate(), ol.get_hit_rate(), olf.periodic_count, count);
        fclose(file);
    }
}

void dc()
{
    int lens[] = {30, 40, 50, 70, 100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000, 100000};
    int num_lens = sizeof(lens) / sizeof(int);
    for (int i = 0; i < num_lens; ++i)
    {
        lru l(lens[i]);
        ours_lru ol(lens[i], 5, 1, 300000);
        lfu lf(lens[i]);
        ours_lfu olf(lens[i], 5, 1, 300000);
        FILE *file = fopen("/share/IMC_DC_IP_trace/univ1/univ1_pt1", "r");
        fseek(file, 24, SEEK_CUR);
        char buf[3000];
        uint32_t last_key = 0;
        uint64_t last_timestamp = 0;
        int count = 0;
        while (1)
        {
            unsigned header_buf[4];
            if (!fread(header_buf, 4, 4, file))
                break;

            time_t ut = header_buf[0];
            unsigned cap_len = header_buf[2];

            if (cap_len > 2048 || cap_len <= 20)
            {
                printf("error\n");
                return;
            }
            fread(buf, cap_len, 1, file);

            ++count;
            if (count == 1000000)
            {
                l.clear_access();
                ol.clear_access();
                lf.clear_access();
                olf.clear_access();
            }

            uint32_t key = Hash::BOBHash32((const uint8_t *)(buf + 12), 13, 1000);
            l.access(key);
            ol.access(key, ut);
            lf.access(key);
            olf.access(key, ut);
        }
        printf("%d,%.9lf,%.9lf,%.9lf,%.9lf,%d,%d\n", lens[i], lf.get_hit_rate(), olf.get_hit_rate(), l.get_hit_rate(), ol.get_hit_rate(), olf.periodic_count, count);
        fclose(file);
    }
}

void caida16()
{
    FILE *file = fopen("/share/datasets/CAIDA2016/formatted00.dat", "r");
    char buf[1000];
    while (1)
    {
        if (!fread(buf, 16, 1, file))
            break;

        uint32_t key = Hash::BOBHash32((const uint8_t *)(buf + 8), 8, 1000);
        uint64_t timestamp;
        memcpy(&timestamp, buf, sizeof(uint64_t));
        ktpv.push_back((ktp){timestamp, key});
        if (ktpv.size() > 10000000)
            break;
    }
    fclose(file);
    sort(ktpv.begin(), ktpv.end(), comp);
    work(20000, 2, 300000);
}

void random_data()
{
    for (int i = 5000; i <= 50000; i += 5000)
    {
        uint64_t count = 0;
        ktpv.clear();
        for (int j = 1; j <= 2000; ++j)
        {
            for (int k = 1; k <= i; ++k)
                ktpv.push_back((ktp){++count, k});
            for (int k = i + 1; k <= 50000; ++k)
                ktpv.push_back((ktp){++count, rand()});
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
    if (strcmp(argv[1], "mawi") == 0)
        mawi();
    else if (strcmp(argv[1], "caida") == 0)
        caida();
    else if (strcmp(argv[1], "webdocs") == 0)
        webdocs();
    else if (strcmp(argv[1], "criteo") == 0)
        criteo();
    else if (strcmp(argv[1], "dc") == 0)
        dc();
    else if (strcmp(argv[1], "caida16") == 0)
        caida16();
    else if (strcmp(argv[1], "random") == 0)
        random_data();
    return 0;
}