#ifndef TOWER_SIMD_H
#define TOWER_SIMD_H

#include <immintrin.h>
#include <numeric>
#include <iostream>
#include <fstream>
#include "Abstract.h"
#include "../common/Util.h"
#include "../common/hash.h"
#include <queue>
#include <vector>
#include <mutex>
#include <atomic>
#include <thread>
#include <chrono>
#include <cstring>
#include <cassert>
#include <random>
#include <unordered_map>
#include <algorithm>
#include "../common/RandomUtil.h"

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

class Matrix {
public:
    uint8_t counter_len, counter_per_int;
    uint64_t w, h, mask;
    uint64_t** counters;  
    uint64_t counter_max;

    uint64_t old_w = 0;
    uint64_t** old_counters = nullptr;
    uint64_t old_actual_w = 0;
    uint32_t transition_cnt = 0;
    const uint32_t transition_window = 2000;

    Matrix() {}
    Matrix(uint64_t in_w, uint64_t in_h, uint8_t in_counter_len) {
        init(in_w, in_h, in_counter_len);
    }
    ~Matrix() { clear(); }

    void init(uint64_t in_w, uint64_t in_h, uint8_t in_counter_len) {
        w = in_w, h = in_h, counter_len = in_counter_len;
        counter_per_int = 64 / counter_len;
        mask = (counter_len == 64) ? 0xFFFFFFFFFFFFFFFFULL : (1ULL << counter_len) - 1;
        counter_max = (counter_len == 64) ? mask : mask+1;

        uint32_t actual_w = (w + counter_per_int - 1) / counter_per_int;
        counters = new uint64_t*[h];
        *counters = new uint64_t[h * actual_w];
        for (uint64_t i = 1; i < h; ++i) {
            counters[i] = *counters + i * actual_w;
        }
    }

    void clear() {
        if (counters) {
            delete[] *counters;
            delete[] counters;
            counters = nullptr;
        }
        if (old_counters) {
            delete[] *old_counters;
            delete[] old_counters;
            old_counters = nullptr;
        }
    }

    bool adjustW(uint64_t target_w) {
        if (target_w == 0 || target_w == w) return false;

        old_w = w;
        old_actual_w = (w + counter_per_int - 1) / counter_per_int;
        old_counters = new uint64_t*[h];
        *old_counters = new uint64_t[h * old_actual_w];
        memcpy(*old_counters, *counters, h * old_actual_w * sizeof(uint64_t));
        for (uint64_t i = 1; i < h; ++i) {
            old_counters[i] = *old_counters + i * old_actual_w;
        }

        uint32_t new_actual_w = (target_w + counter_per_int - 1) / counter_per_int;
        delete[] *counters;
        delete[] counters;
        counters = new uint64_t*[h];
        *counters = new uint64_t[h * new_actual_w]();
        for (uint64_t i = 1; i < h; ++i) {
            counters[i] = *counters + i * new_actual_w;
        }

        w = target_w;
        transition_cnt = transition_window;
        return true;
    }

    bool adjustWByStep(int32_t step) {
        uint64_t target_w = w + step;
        if (target_w < 1) target_w = 1;
        return adjustW(target_w);
    }

    inline TIME_TYPE query(DATA_TYPE key, uint32_t row_id, uint8_t i) {
        uint32_t index = ::hash(key, i) % w;
        uint64_t buf = counters[row_id][index / counter_per_int];
        uint8_t shift = (index % counter_per_int) * counter_len;
        TIME_TYPE val = (buf >> shift) & mask;

        if (old_w > 0 && transition_cnt > 0 && val == 0) {
            uint32_t old_index = ::hash(key, i) % old_w;
            uint64_t old_buf = old_counters[row_id][old_index / counter_per_int];
            uint8_t old_shift = (old_index % counter_per_int) * counter_len;
            val = (old_buf >> old_shift) & mask;
        }

        return val;
    }

    inline void insert(DATA_TYPE key, uint32_t row_id, uint8_t i, TIME_TYPE f) {
        uint32_t index = ::hash(key, i) % w;
        uint64_t buf = counters[row_id][index / counter_per_int];
        uint8_t shift = (index % counter_per_int) * counter_len;
        TIME_TYPE remain = f % counter_max;
        counters[row_id][index / counter_per_int] = (buf & ~(mask << shift)) | (remain << shift);

        if (old_w > 0 && transition_cnt > 0) {
            uint32_t old_index = ::hash(key, i) % old_w;
            uint64_t old_buf = old_counters[row_id][old_index / counter_per_int];
            uint8_t old_shift = (old_index % counter_per_int) * counter_len;
            old_counters[row_id][old_index / counter_per_int] = (old_buf & ~(mask << old_shift)) | (remain << old_shift);
        }

        if (old_w > 0 && --transition_cnt == 0) {
            delete[] *old_counters;
            delete[] old_counters;
            old_counters = nullptr;
            old_w = 0;
        }
    }

    uint32_t memoryUsage() const {
        uint32_t actual_w = (w + counter_per_int - 1) / counter_per_int;
        return h * actual_w * sizeof(uint64_t);
    }
};

class Ours_tower : public Abstract {
public:        
struct Bucket {
    struct Cell {
            DATA_TYPE item = 0;
            uint32_t interval = 0;
            COUNT_TYPE count = 0;
        };
    DATA_TYPE* items = nullptr;
    uint32_t* intervals = nullptr;
    COUNT_TYPE* counts = nullptr;
    uint8_t cell_num = 0;
    COUNT_TYPE fail = 0;

    Bucket() : cell_num(8) {
        items = new DATA_TYPE[cell_num]();
        intervals = new uint32_t[cell_num]();
        counts = new COUNT_TYPE[cell_num]();
    }

    Bucket(uint8_t init_cell_num) : cell_num(init_cell_num) {
        assert(init_cell_num >= MIN_CELL_NUM && init_cell_num <= MAX_CELL_NUM);
        items = new DATA_TYPE[cell_num]();
        intervals = new uint32_t[cell_num]();
        counts = new COUNT_TYPE[cell_num]();
    }

    ~Bucket() {
        delete[] items;
        delete[] intervals;
        delete[] counts;
        items = nullptr;
        intervals = nullptr;
        counts = nullptr;
    }

    Bucket(Bucket&& other) noexcept 
        : items(other.items), intervals(other.intervals), counts(other.counts),
          cell_num(other.cell_num), fail(other.fail) {
        other.items = nullptr;
        other.intervals = nullptr;
        other.counts = nullptr;
        other.cell_num = 0;
    }

    Bucket& operator=(Bucket&& other) noexcept {
        if (this != &other) {
            delete[] items;
            delete[] intervals;
            delete[] counts;
            items = other.items;
            intervals = other.intervals;
            counts = other.counts;
            cell_num = other.cell_num;
            fail = other.fail;
            other.items = nullptr;
            other.intervals = nullptr;
            other.counts = nullptr;
            other.cell_num = 0;
        }
        return *this;
    }

    Bucket(const Bucket&) = delete;
    Bucket& operator=(const Bucket&) = delete;
};

    static constexpr uint32_t REFLOW_INTERVAL = 200000;
    static constexpr uint32_t MIN_CELL_NUM = 4;
    static constexpr uint32_t MAX_CELL_NUM = 32;

    uint8_t level = 0;
    Matrix* mat = nullptr;
    std::vector<double> layer_ratios;

    Bucket* buckets = nullptr;
    uint32_t HEAVY_LENGTH = 0;

    std::atomic<uint32_t> total_inserts{0};
    std::atomic<bool> is_reflowing{false};

    struct alignas(4) BucketStats {
        std::atomic<uint16_t> attempts{0};
        std::atomic<uint16_t> hits{0};

        void reset() {
            attempts.store(0, std::memory_order_relaxed);
            hits.store(0, std::memory_order_relaxed);
        }

        double hitRate() const {
            uint16_t a = attempts.load(std::memory_order_relaxed);
            uint16_t h = hits.load(std::memory_order_relaxed);
            return (a > 0) ? static_cast<double>(h) / a : 0.0;
        }
    };
    BucketStats* bucket_stats;

    Ours_tower(uint32_t _MEMORY , double _ratio = 0.1,
               const std::vector<double>& tower_ratios = {},
               const std::vector<uint8_t>& counter_len = {},
               std::string _name = "Ours_tower")  
    {
        assert(tower_ratios.size() == counter_len.size() && "分层比例与计数器长度数量需一致");
        level = tower_ratios.size();
        layer_ratios = tower_ratios;

        MEMORY = _MEMORY;
        NAME = _name;
        uint32_t light_memory = static_cast<uint32_t>(MEMORY * _ratio);
        uint32_t heavy_memory = MEMORY - light_memory;

        mat = new Matrix[level];
        for (uint8_t j = 0; j < level; ++j) {
            uint32_t layer_light = static_cast<uint32_t>(light_memory * layer_ratios[j]);
            uint64_t w = (layer_light * 8) / counter_len[j];
            mat[j].init(w, 1, counter_len[j]);
        }

        uint32_t avg_cell_per_bucket = 8;
        uint32_t cell_size = sizeof(Bucket::Cell);
        HEAVY_LENGTH = heavy_memory / (avg_cell_per_bucket * cell_size+4+16);
        if (HEAVY_LENGTH == 0) HEAVY_LENGTH = 1;
        bucket_stats = new BucketStats[HEAVY_LENGTH]();
        buckets = reinterpret_cast<Bucket*>(new char[HEAVY_LENGTH * sizeof(Bucket)]);
        for (uint32_t i = 0; i < HEAVY_LENGTH; ++i) {
            new (&buckets[i]) Bucket(8);
        }
    }

~Ours_tower() {
    if (mat) {
        for (uint8_t j = 0; j < level; ++j) mat[j].clear();
        delete[] mat;
        mat = nullptr;
    }
    if (buckets) {
        for (uint32_t i = 0; i < HEAVY_LENGTH; ++i) {
            buckets[i].~Bucket();
        }
        delete[] reinterpret_cast<char*>(buckets);
        buckets = nullptr;
    }
    if (bucket_stats) {
        delete[] bucket_stats;
        bucket_stats = nullptr;
    }
}

    void Insert(const ItemPair& item) {
        TIME_TYPE maxinterval = 0;
        for (uint8_t i = 0; i < level; ++i) {
            TIME_TYPE last_time = mat[i].query(item.item, 0, i);
            TIME_TYPE interval = calculateInterval(last_time, item.time, mat[i].counter_max);
            maxinterval = std::max(maxinterval, interval);
            mat[i].insert(item.item, 0, i, item.time);
        }

        if (maxinterval > 500) return;

        ItemPair temp(maxinterval, item.item);
        uint32_t pos = ::hash(temp, 2) % HEAVY_LENGTH;

        bool is_hit = processBucket(buckets[pos], temp);
        
        updateBucketStats(pos, is_hit);

        maybeTriggerReflow();
    }

 HashMap Report1(COUNT_TYPE threshold) {
    HashMap ret;
    __m256i v_threshold = _mm256_set1_epi32(threshold);

    for (uint32_t i = 0; i < HEAVY_LENGTH; ++i) {
        Bucket& bucket = buckets[i];
        uint32_t cell_num = bucket.cell_num;
        
        COUNT_TYPE* counts = bucket.counts;
        DATA_TYPE* items = bucket.items;
        uint32_t* intervals = bucket.intervals;

        uint32_t j = 0;
        for (; j <= cell_num - 8; j += 8) {
            __m256i v_counts = _mm256_load_si256(reinterpret_cast<const __m256i*>(&counts[j]));
            __m256i v_mask = _mm256_cmpgt_epi32(v_counts, v_threshold);
            if (_mm256_testz_si256(v_mask, v_mask)) continue;

            int mask = _mm256_movemask_epi8(v_mask);
            for (int k = 0; k < 8; ++k) {
                if (mask & (0xF << (k * 4))) {
                    ret[ItemPair(intervals[j + k], items[j + k])] = counts[j + k];
                }
            }
        }

        for (; j < cell_num; ++j) {
            if (counts[j] > threshold) {
                ret[ItemPair(intervals[j], items[j])] = counts[j];
            }
        }
    }
    return ret;
}

private:
    TIME_TYPE calculateInterval(TIME_TYPE last, TIME_TYPE current, TIME_TYPE mod) {
        if (current % mod >= last) return current % mod - last;
        else return current % mod + mod - last;
    }

bool processBucket(Bucket& bucket, const ItemPair& item) {
    const uint32_t cell_num = bucket.cell_num;
    
    #if defined(__AVX2__)
    if (cell_num >= 8) {
        __m256i target_item = _mm256_set1_epi64x(item.item);
        
        uint32_t simd_end = (cell_num / 4) * 4;
        if (simd_end > cell_num - 4) {
            simd_end = (cell_num >= 4) ? (cell_num - 4) : 0;
        }
        
        for (uint32_t j = 0; j < simd_end; j += 4) {
            if (j + 3 >= cell_num) break;
            
            __m256i items_vec = _mm256_loadu_si256((__m256i*)&bucket.items[j]);
            __m256i items_cmp = _mm256_cmpeq_epi64(items_vec, target_item);
            
            __m128i intervals_128 = _mm_loadu_si128((__m128i*)&bucket.intervals[j]);
            __m256i intervals_vec = _mm256_cvtepu32_epi64(intervals_128);
            
            __m128i target_interval_128 = _mm_set1_epi32(item.time);
            __m256i target_interval_vec = _mm256_cvtepu32_epi64(target_interval_128);
            __m256i intervals_cmp = _mm256_cmpeq_epi64(intervals_vec, target_interval_vec);
            
            __m256i match_mask = _mm256_and_si256(items_cmp, intervals_cmp);
            int mask = _mm256_movemask_epi8(match_mask);
            
            if (mask != 0) {
                int bit_pos = __builtin_ctz(mask);
                uint32_t match_idx = j + (bit_pos / 8);
                
                if (match_idx < cell_num) {
                    bucket.counts[match_idx]++;
                    return true;
                }
            }
        }
        
        for (uint32_t j = simd_end; j < cell_num; ++j) {
            if (bucket.items[j] == item.item && bucket.intervals[j] == item.time) {
                bucket.counts[j]++;
                return true;
            }
        }
        
        uint32_t min_pos = 0;
        COUNT_TYPE min_count = INT32_MAX;
        
        uint32_t simd_end_8 = (cell_num / 8) * 8;
        
        if (simd_end_8 > 0 && cell_num >= 8) {
            __m256i min_counts_vec = _mm256_set1_epi32(INT32_MAX);
            __m256i min_indices_vec = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
            __m256i current_min_idx = min_indices_vec;
            
            for (uint32_t j = 0; j < simd_end_8; j += 8) {
                if (j + 7 >= cell_num) break;
                
                __m256i counts_vec = _mm256_loadu_si256((__m256i*)&bucket.counts[j]);
                __m256i indices_vec = _mm256_add_epi32(min_indices_vec, _mm256_set1_epi32(j));
                
                __m256i cmp_mask = _mm256_cmpgt_epi32(min_counts_vec, counts_vec);
                min_counts_vec = _mm256_blendv_epi8(min_counts_vec, counts_vec, cmp_mask);
                current_min_idx = _mm256_blendv_epi8(current_min_idx, indices_vec, cmp_mask);
            }
            
            alignas(32) COUNT_TYPE counts_arr[8];
            alignas(32) uint32_t indices_arr[8];
            _mm256_store_si256((__m256i*)counts_arr, min_counts_vec);
            _mm256_store_si256((__m256i*)indices_arr, current_min_idx);
            
            for (int i = 0; i < 8; ++i) {
                if (counts_arr[i] < min_count) {
                    min_count = counts_arr[i];
                    min_pos = indices_arr[i];
                }
            }
            
            for (uint32_t j = simd_end_8; j < cell_num; ++j) {
                if (bucket.counts[j] < min_count) {
                    min_count = bucket.counts[j];
                    min_pos = j;
                }
            }
        } else {
            for (uint32_t j = 0; j < cell_num; ++j) {
                if (bucket.counts[j] < min_count) {
                    min_count = bucket.counts[j];
                    min_pos = j;
                }
            }
        }
        
        if (min_count == 0) {
            bucket.items[min_pos] = item.item;
            bucket.intervals[min_pos] = item.time;
            bucket.counts[min_pos] = 1;
            return true;
        }
        
        COUNT_TYPE f_min = min_count;
        double threshold = 1.0 / (2 * f_min - bucket.fail + 1);
        std::uniform_real_distribution<> dist(0.0, 1.0);
        double rand_val = dist(getRNG());
        
        if (rand_val <= threshold) {
            bucket.items[min_pos] = item.item;
            bucket.intervals[min_pos] = item.time;
            bucket.counts[min_pos] = f_min + bucket.fail / f_min;
            bucket.fail = 0;
        } else {
            bucket.fail++;
        }
        
        return false;
        
    } else {
    #endif
        
        bool is_hit = false;
        uint32_t min_pos = 0;
        COUNT_TYPE min_count = INT32_MAX;

        for (uint32_t j = 0; j < cell_num; ++j) {
            if (bucket.items[j] == item.item && bucket.intervals[j] == item.time) {
                bucket.counts[j]++;
                is_hit = true;
                return is_hit;
            }
            if (bucket.counts[j] < min_count) {
                min_count = bucket.counts[j];
                if (!min_count) {
                    bucket.items[j] = item.item;
                    bucket.intervals[j] = item.time;
                    bucket.counts[j] = 1;
                    is_hit = true;
                    return is_hit;
                }
                min_pos = j;
            }
        }

        if (!is_hit) {
            COUNT_TYPE f_min = min_count;
            double threshold = 1.0 / (2 * f_min - bucket.fail + 1);
            std::uniform_real_distribution<> dist(0.0, 1.0);
            double rand_val = dist(getRNG());
            if (rand_val <= threshold) {
                bucket.items[min_pos] = item.item;
                bucket.intervals[min_pos] = item.time;
                bucket.counts[min_pos] = f_min + bucket.fail / f_min;
                bucket.fail = 0;
            } else {
                bucket.fail++;
            }
        }
        return is_hit;
        
    #if defined(__AVX2__)
    }
    #endif
}
   
    void updateBucketStats(uint32_t bucket_idx, bool is_hit) {
        bucket_stats[bucket_idx].attempts.fetch_add(1, std::memory_order_relaxed);
        if (is_hit) {
            bucket_stats[bucket_idx].hits.fetch_add(1, std::memory_order_relaxed);
        }
    }

    void maybeTriggerReflow() {
        uint32_t inserts = total_inserts.fetch_add(1, std::memory_order_acq_rel) + 1;
        if (inserts % REFLOW_INTERVAL == 0) {
            triggerBatchReflow();
        }
    }

std::pair<int, int> computeZTestDeltas(
    std::vector<int32_t>& deltas,
    uint32_t current_total_cells,
    double z_threshold_expand = 0,
    double z_threshold_shrink = 0
) {
    z_threshold_expand = (MEMORY < 100000) ? 1.45: 
                        ((MEMORY < 150000) ? 1.811 : 2.5);
    z_threshold_shrink = (MEMORY < 100000) ? -1.35: 
                        ((MEMORY < 150000) ? -1.609 : -5);
    
    deltas.assign(HEAVY_LENGTH, 0);
    
    #if defined(__AVX2__) && defined(__FMA__)
    alignas(32) double densities[HEAVY_LENGTH];
    double sum_density = 0.0;
    int valid_count = 0;
    
    __m256d sum_vec = _mm256_setzero_pd();
    const __m256d epsilon = _mm256_set1_pd(1e-9);
    
    uint32_t simd_end = (HEAVY_LENGTH / 4) * 4;
    
    for (uint32_t i = 0; i < simd_end; i += 4) {
        __m128i cell_nums_i32 = _mm_set_epi32(
            buckets[i+3].cell_num,
            buckets[i+2].cell_num,
            buckets[i+1].cell_num,
            buckets[i].cell_num
        );
        __m256d cell_nums = _mm256_cvtepi32_pd(cell_nums_i32);
        
        __m128i attempts_i16 = _mm_set_epi16(
            0, 0, 0, 0,
            bucket_stats[i+3].attempts.load(std::memory_order_relaxed),
            bucket_stats[i+2].attempts.load(std::memory_order_relaxed),
            bucket_stats[i+1].attempts.load(std::memory_order_relaxed),
            bucket_stats[i].attempts.load(std::memory_order_relaxed)
        );
        __m128i hits_i16 = _mm_set_epi16(
            0, 0, 0, 0,
            bucket_stats[i+3].hits.load(std::memory_order_relaxed),
            bucket_stats[i+2].hits.load(std::memory_order_relaxed),
            bucket_stats[i+1].hits.load(std::memory_order_relaxed),
            bucket_stats[i].hits.load(std::memory_order_relaxed)
        );
        
        __m128i conflicts_i16 = _mm_sub_epi16(attempts_i16, hits_i16);
        __m128i conflicts_i32 = _mm_cvtepu16_epi32(conflicts_i16);
        __m256d conflicts = _mm256_cvtepi32_pd(conflicts_i32);
        
        __m256d is_valid = _mm256_cmp_pd(cell_nums, epsilon, _CMP_GT_OQ);
        __m256d density = _mm256_div_pd(conflicts, _mm256_max_pd(cell_nums, epsilon));
        
        __m256d masked_density = _mm256_and_pd(density, is_valid);
        sum_vec = _mm256_add_pd(sum_vec, masked_density);
        
        _mm256_store_pd(&densities[i], density);
        
        int mask = _mm256_movemask_pd(is_valid);
        valid_count += __builtin_popcount(mask);
    }
    
    for (uint32_t i = simd_end; i < HEAVY_LENGTH; ++i) {
        uint32_t cell_num = buckets[i].cell_num;
        if (cell_num == 0) {
            densities[i] = 0.0;
            continue;
        }
        uint32_t conflicts = bucket_stats[i].attempts.load(std::memory_order_relaxed) 
                           - bucket_stats[i].hits.load(std::memory_order_relaxed);
        densities[i] = static_cast<double>(conflicts) / cell_num;
        sum_density += densities[i];
        valid_count++;
    }
    
    alignas(32) double sum_arr[4];
    _mm256_store_pd(sum_arr, sum_vec);
    sum_density += sum_arr[0] + sum_arr[1] + sum_arr[2] + sum_arr[3];
    
    if (valid_count == 0) {
        return {0, 0};
    }
    
    double avg_density = sum_density / valid_count;
    __m256d avg_vec = _mm256_set1_pd(avg_density);
    __m256d sum_sq_vec = _mm256_setzero_pd();
    
    for (uint32_t i = 0; i < simd_end; i += 4) {
        __m256d density_vec = _mm256_load_pd(&densities[i]);
        __m256d diff = _mm256_sub_pd(density_vec, avg_vec);
        __m256d sq_diff = _mm256_mul_pd(diff, diff);
        
        __m128i cell_nums_i32 = _mm_set_epi32(
            buckets[i+3].cell_num,
            buckets[i+2].cell_num,
            buckets[i+1].cell_num,
            buckets[i].cell_num
        );
        __m256d cell_nums = _mm256_cvtepi32_pd(cell_nums_i32);
        __m256d is_valid = _mm256_cmp_pd(cell_nums, epsilon, _CMP_GT_OQ);
        __m256d masked_sq_diff = _mm256_and_pd(sq_diff, is_valid);
        
        sum_sq_vec = _mm256_add_pd(sum_sq_vec, masked_sq_diff);
    }
    
    double sum_sq_diff = 0.0;
    for (uint32_t i = simd_end; i < HEAVY_LENGTH; ++i) {
        if (buckets[i].cell_num == 0) continue;
        double diff = densities[i] - avg_density;
        sum_sq_diff += diff * diff;
    }
    
    alignas(32) double sq_arr[4];
    _mm256_store_pd(sq_arr, sum_sq_vec);
    sum_sq_diff += sq_arr[0] + sq_arr[1] + sq_arr[2] + sq_arr[3];
    
    double std_density = std::sqrt(sum_sq_diff / valid_count);
    if (std_density < 1e-9) {
        std_density = 1.0;
    }
    
    int expand_count = 0;
    int shrink_count = 0;
    
    __m256d std_vec = _mm256_set1_pd(std_density);
    __m256d z_expand_vec = _mm256_set1_pd(z_threshold_expand);
    __m256d z_shrink_vec = _mm256_set1_pd(z_threshold_shrink);
    
    for (uint32_t i = 0; i < simd_end; i += 4) {
        __m256d density_vec = _mm256_load_pd(&densities[i]);
        __m256d diff = _mm256_sub_pd(density_vec, avg_vec);
        __m256d z_vec = _mm256_div_pd(diff, std_vec);
        
        __m128i cell_nums_i32 = _mm_set_epi32(
            buckets[i+3].cell_num,
            buckets[i+2].cell_num,
            buckets[i+1].cell_num,
            buckets[i].cell_num
        );
        
        alignas(32) double z_arr[4];
        _mm256_store_pd(z_arr, z_vec);
        
        for (int j = 0; j < 4; ++j) {
            uint32_t idx = i + j;
            uint32_t cell_num = buckets[idx].cell_num;
            if (cell_num == 0) continue;
            
            double z = z_arr[j];
            int32_t delta = 0;
            int32_t max_delta = std::max(1, static_cast<int32_t>(cell_num * 0.25));
            max_delta = MIN(max_delta, 4);
            
            if (z > z_threshold_expand) {
                delta = max_delta;
                expand_count++;
            } else if (z < z_threshold_shrink) {
                delta = -max_delta;
                shrink_count++;
            }
            
            deltas[idx] = delta;
        }
    }
    
    for (uint32_t i = simd_end; i < HEAVY_LENGTH; ++i) {
        uint32_t cell_num = buckets[i].cell_num;
        if (cell_num == 0) continue;
        
        double z = (densities[i] - avg_density) / std_density;
        int32_t delta = 0;
        int32_t max_delta = std::max(1, static_cast<int32_t>(cell_num * 0.25));
        max_delta = MIN(max_delta, 4);
        
        if (z > z_threshold_expand) {
            delta = max_delta;
            expand_count++;
        } else if (z < z_threshold_shrink) {
            delta = -max_delta;
            shrink_count++;
        }
        
        deltas[i] = delta;
    }
    
    return {expand_count, shrink_count};
    
    #else
    std::vector<double> densities;
    double sum_density = 0.0;
    int valid_count = 0;

    for (uint32_t i = 0; i < HEAVY_LENGTH; ++i) {
        uint32_t cell_num = buckets[i].cell_num;
        if (cell_num == 0) continue;

        uint32_t conflicts = bucket_stats[i].attempts.load(std::memory_order_relaxed)
                           - bucket_stats[i].hits.load(std::memory_order_relaxed);
        double density = static_cast<double>(conflicts) / cell_num;

        densities.push_back(density);
        sum_density += density;
        valid_count++;
    }

    if (valid_count == 0) {
        return {0, 0};
    }

    double avg_density = sum_density / valid_count;
    double sum_sq_diff = 0.0;
    for (double d : densities) {
        double diff = d - avg_density;
        sum_sq_diff += diff * diff;
    }
    double std_density = std::sqrt(sum_sq_diff / valid_count);

    if (std_density < 1e-9) {
        std_density = 1.0;
    }

    int expand_count = 0;
    int shrink_count = 0;

    size_t idx = 0;
    for (uint32_t i = 0; i < HEAVY_LENGTH; ++i) {
        uint32_t cell_num = buckets[i].cell_num;
        if (cell_num == 0) continue;

        double density = densities[idx++];
        double z = (density - avg_density) / std_density;

        int32_t delta = 0;
        int32_t max_delta = std::max(1, static_cast<int32_t>(cell_num * 0.25));
        max_delta = MIN(max_delta, 4);

        if (z > z_threshold_expand) {
            delta = max_delta;
            expand_count++;
        } else if (z < z_threshold_shrink) {
            delta = -max_delta;
            shrink_count++;
        }

        deltas[i] = delta;
    }

    return {expand_count, shrink_count};
    #endif
}
void triggerBatchReflow() {
    if (is_reflowing) {
        return;
    }
    is_reflowing = true;

    auto start = std::chrono::high_resolution_clock::now();

    uint32_t current_total_cells = 0;
    std::vector<double> hit_rates;
    for (uint32_t i = 0; i < HEAVY_LENGTH; ++i) {
        current_total_cells += buckets[i].cell_num;
    }
    
    std::vector<int32_t> deltas;
    auto [expand_count, shrink_count] = computeZTestDeltas(deltas, current_total_cells);

    const double USAGE_RATIO = 0.95;
    uint64_t max_allowed_bytes = static_cast<uint64_t>(MEMORY * USAGE_RATIO);
    uint32_t max_allowed_cells = static_cast<uint32_t>((max_allowed_bytes-20*HEAVY_LENGTH) / sizeof(Bucket::Cell));
    if (max_allowed_cells < current_total_cells) {
        max_allowed_cells = current_total_cells;
    }
    int32_t available_increase = std::max(static_cast<int32_t>(max_allowed_cells - current_total_cells), 0);

    int32_t total_cell_change = 0;
    for (uint32_t i = 0; i < HEAVY_LENGTH; ++i) {
        int32_t delta = deltas[i];
        uint32_t current_num = buckets[i].cell_num;

        if (delta == 0) continue;

        int32_t new_num = static_cast<int32_t>(current_num) + delta;
        new_num = (new_num < static_cast<int32_t>(MIN_CELL_NUM)) ? static_cast<int32_t>(MIN_CELL_NUM) : new_num;
        new_num = (new_num > static_cast<int32_t>(MAX_CELL_NUM)) ? static_cast<int32_t>(MAX_CELL_NUM) : new_num;
        int32_t net_change = new_num - static_cast<int32_t>(current_num);

        if (net_change > 0 && total_cell_change + net_change > available_increase) {
            continue;
        }

        if (buckets[i].cell_num != current_num) {
            continue;
        }
        resizeBucketCells(buckets[i], static_cast<uint32_t>(new_num));
        total_cell_change += net_change;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    int32_t delta_memory = total_cell_change * sizeof(Bucket::Cell);

    for (int j = 0; j < level; ++j) {
        double ratio = layer_ratios[j];
        double counter_bytes = mat[j].counter_len / 8.0;
        int32_t memory_for_layer = static_cast<int32_t>(delta_memory * ratio);
        int32_t step = static_cast<int32_t>(memory_for_layer / counter_bytes);      
        mat[j].adjustWByStep(-step);      
    }

        for (uint32_t i = 0; i < HEAVY_LENGTH; ++i) {
            bucket_stats[i].reset();
        }
        is_reflowing = false;
    
}

void resizeBucketCells(Bucket& bucket, uint32_t target_num) {
    uint32_t old_num = bucket.cell_num;
    if (target_num == old_num) return;
    if (target_num < MIN_CELL_NUM || target_num > MAX_CELL_NUM) return;

    DATA_TYPE* old_items = bucket.items;
    uint32_t* old_intervals = bucket.intervals;
    COUNT_TYPE* old_counts = bucket.counts;

    DATA_TYPE* new_items = new DATA_TYPE[target_num]();
    uint32_t* new_intervals = new uint32_t[target_num]();
    COUNT_TYPE* new_counts = new COUNT_TYPE[target_num]();

    if (target_num > old_num) {
        std::copy(old_items, old_items + old_num, new_items);
        std::copy(old_intervals, old_intervals + old_num, new_intervals);
        std::copy(old_counts, old_counts + old_num, new_counts);
    } else {
        std::vector<uint32_t> indices(old_num);
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](uint32_t a, uint32_t b) {
            return old_counts[a] > old_counts[b];
        });
        for (uint32_t i = 0; i < target_num; ++i) {
            uint32_t idx = indices[i];
            new_items[i] = old_items[idx];
            new_intervals[i] = old_intervals[idx];
            new_counts[i] = old_counts[idx];
        }
    }

    delete[] old_items;
    delete[] old_intervals;
    delete[] old_counts;
    bucket.items = new_items;
    bucket.intervals = new_intervals;
    bucket.counts = new_counts;
    bucket.cell_num = target_num;
    bucket.fail = 0;
}

    inline std::mt19937& getRNG() {
        thread_local std::mt19937 rng = getLocalRNG();
        return rng;
    }
};

#endif // TOWER_H