
#include "test_suite.h"

#include <iostream>
#include <chrono>
#include <random>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sys/time.h>

// #include "../BTreeOLC/BTreeOLC.h"
#include "../BTreeOLC/BTreeOLC_child_layout.h"

#include<cilk/cilk.h>
#include<thread>

static long get_usecs() {
  struct timeval st;
  gettimeofday(&st, NULL);
  return st.tv_sec * 1000000 + st.tv_usec;
}

using namespace std;

using Key=uint64_t;

using TID=uint64_t;

// index types
enum {
    TYPE_BTREEOLC,
    TYPE_ART,
    TYPE_HOT,
    TYPE_BWTREE,
    TYPE_MASSTREE,
    TYPE_CLHT,
    TYPE_FASTFAIR,
    TYPE_LEVELHASH,
    TYPE_CCEH,
    TYPE_WOART,
};

enum {
    OP_INSERT,
    OP_UPDATE,
    OP_READ,
    OP_SCAN,
    OP_SCAN_END,
    OP_DELETE,
};

enum {
    WORKLOAD_A,
    WORKLOAD_B,
    WORKLOAD_C,
    WORKLOAD_D,
    WORKLOAD_E,
    WORKLOAD_X,
    WORKLOAD_Y,
};

enum {
    RANDINT_KEY,
    STRING_KEY,
};

enum {
    UNIFORM,
    ZIPFIAN,
};

namespace Dummy {
    inline void mfence() {asm volatile("mfence":::"memory");}

    inline void clflush(char *data, int len, bool front, bool back)
    {
        if (front)
            mfence();
        volatile char *ptr = (char *)((unsigned long)data & ~(64 - 1));
        for (; ptr < data+len; ptr += 64){
#ifdef CLFLUSH
            asm volatile("clflush %0" : "+m" (*(volatile char *)ptr));
#elif CLFLUSH_OPT
            asm volatile(".byte 0x66; clflush %0" : "+m" (*(volatile char *)(ptr)));
#elif CLWB
            asm volatile(".byte 0x66; xsaveopt %0" : "+m" (*(volatile char *)(ptr)));
#endif
        }
        if (back)
            mfence();
    }
}


static uint64_t LOAD_SIZE = 100000000;
static uint64_t RUN_SIZE = 100000000;

void loadKey(TID tid, Key &key) {
    return ;
}

void ycsb_load_run_string(int index_type, int wl, int kt, int ap, int num_thread,
        std::vector<Key *> &init_keys,
        std::vector<Key *> &keys,
        std::vector<int> &ranges,
        std::vector<int> &ops)
{
    
}

template <typename F> inline void parallel_for(size_t start, size_t end, F f) {
  cilk_for(size_t i = start; i < end; i++) f(i);
}

double findMedian(vector<double>& vec) {
    size_t size = vec.size();
    if (size == 0) {
        return 0;
    }
    sort(vec.begin(), vec.end());
    if (size % 2 == 0) {
        return (vec[size / 2 - 1] + vec[size / 2]) / 2;
    } else {
        return vec[size / 2];
    }
}

template <class T>
std::vector<T> create_random_data(size_t n, size_t max_val,
                                  std::seed_seq &seed) {

  std::mt19937_64 eng(seed); // a source of random data

  std::uniform_int_distribution<T> dist(0, max_val);
  std::vector<T> v(n);

  generate(begin(v), end(v), bind(dist, eng));
  return v;
}

void ycsb_load_run_randint(int index_type, int wl, int kt, int ap, int num_thread,
        std::vector<uint64_t> &init_keys,
        std::vector<uint64_t> &keys,
        std::vector<uint64_t> &range_end,
        std::vector<int> &ranges,
        std::vector<int> &ops)
{
    std::string init_file;
    std::string txn_file;

    if (ap == UNIFORM) {
        if (kt == RANDINT_KEY && wl == WORKLOAD_A) {
            init_file = "../ycsb/index-microbench/workloads/uniform/loada_unif_int.dat";
            txn_file = "../ycsb/index-microbench/workloads/uniform/txnsa_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_B) {
            init_file = "../ycsb/index-microbench/workloads/uniform/loadb_unif_int.dat";
            txn_file = "../ycsb/index-microbench/workloads/uniform/txnsb_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_C) {
            init_file = "../ycsb/index-microbench/workloads/uniform/loadc_unif_int.dat";
            txn_file = "../ycsb/index-microbench/workloads/uniform/txnsc_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_D) {
            init_file = "../ycsb/index-microbench/workloads/uniform/loadd_unif_int.dat";
            txn_file = "../ycsb/index-microbench/workloads/uniform/txnsd_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_E) {
            init_file = "../ycsb/index-microbench/workloads/uniform/loade_unif_int.dat";
            txn_file = "../ycsb/index-microbench/workloads/uniform/txnse_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_X) {
            init_file = "../ycsb/index-microbench/workloads/uniform/loadx_unif_int.dat";
            txn_file = "../ycsb/index-microbench/workloads/uniform/txnsx_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_Y) {
            init_file = "../ycsb/index-microbench/workloads/uniform/loady_unif_int.dat";
            txn_file = "../ycsb/index-microbench/workloads/uniform/txnsy_unif_int.dat";
        }
    } else {
        if (kt == RANDINT_KEY && wl == WORKLOAD_A) {
            init_file = "../ycsb/index-microbench/workloads/zipfian/loada_unif_int.dat";
            txn_file = "../ycsb/index-microbench/workloads/zipfian/txnsa_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_B) {
            init_file = "../ycsb/index-microbench/workloads/zipfian/loadb_unif_int.dat";
            txn_file = "../ycsb/index-microbench/workloads/zipfian/txnsb_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_C) {
            init_file = "../ycsb/index-microbench/workloads/zipfian/loadc_unif_int.dat";
            txn_file = "../ycsb/index-microbench/workloads/zipfian/txnsc_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_D) {
            init_file = "../ycsb/index-microbench/workloads/zipfian/loadd_unif_int.dat";
            txn_file = "../ycsb/index-microbench/workloads/zipfian/txnsd_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_E) {
            init_file = "../ycsb/index-microbench/workloads/zipfian/loade_unif_int.dat";
            txn_file = "../ycsb/index-microbench/workloads/zipfian/txnse_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_X) {
            init_file = "../ycsb/index-microbench/workloads/zipfian/loadx_unif_int.dat";
            txn_file = "../ycsb/index-microbench/workloads/zipfian/txnsx_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_Y) {
            init_file = "../ycsb/index-microbench/workloads/zipfian/loadx_unif_int.dat";
            txn_file = "../ycsb/index-microbench/workloads/zipfian/txnsy_unif_int.dat";
        }
    }

    std::ifstream infile_load(init_file);

    std::string op;
    uint64_t key;
    uint64_t rend;
    int range;

    std::string insert("INSERT");
    std::string update("UPDATE");
    std::string read("READ");
    std::string scan("SCAN");
    std::string scanend("SCANEND");

    int count = 0;
    while ((count < LOAD_SIZE) && infile_load.good()) {
        infile_load >> op >> key;
        if (op.compare(insert) != 0) {
            std::cout << "READING LOAD FILE FAIL!\n";
            return ;
        }
        init_keys.push_back(key);
        count++;
    }

    fprintf(stderr, "Loaded %d keys\n", count);

    std::ifstream infile_txn(txn_file);

    count = 0;
    while ((count < RUN_SIZE) && infile_txn.good()) {
        infile_txn >> op >> key;
        if (op.compare(insert) == 0) {
            ops.push_back(OP_INSERT);
            keys.push_back(key);
            ranges.push_back(1);
            range_end.push_back(1);
        } else if (op.compare(update) == 0) {
            ops.push_back(OP_UPDATE);
            keys.push_back(key);
            ranges.push_back(1);
            range_end.push_back(1);
        } else if (op.compare(read) == 0) {
            ops.push_back(OP_READ);
            keys.push_back(key);
            ranges.push_back(1);
            range_end.push_back(1);
        } else if (op.compare(scan) == 0) {
            infile_txn >> range;
            ops.push_back(OP_SCAN);
            keys.push_back(key);
            ranges.push_back(range);
            range_end.push_back(1);
        } else if (op.compare(scanend) == 0) {
            infile_txn >> rend;
            ops.push_back(OP_SCAN_END);
            keys.push_back(key);
            range_end.push_back(rend);
            ranges.push_back(1);
        } else {
            std::cout << "UNRECOGNIZED CMD!\n";
            return;
        }
        count++;
    }

    std::atomic<int> range_complete, range_incomplete;
    range_complete.store(0);
    range_incomplete.store(0);

    fprintf(stderr, "Loaded %d more keys\n", count);

    std::this_thread::sleep_for(std::chrono::nanoseconds(3000000000));

    fprintf(stderr, "Slept\n");
    
    // /*
    vector<double> load_tpts;
    vector<double> run_tpts;

    for (int trial = 0; trial <= 5; trial++) {

        if (index_type == TYPE_BWTREE) {
            TreeType *t = GetEmptyTree();
            t->UpdateThreadLocal(1);
            t->AssignGCID(0);
            std::atomic<int> next_thread_id;

            std::vector<uint64_t> query_results_keys(RUN_SIZE);
            std::vector<uint64_t> query_results_vals(RUN_SIZE);

            {
                // Load
                auto starttime = std::chrono::system_clock::now();
                next_thread_id.store(0);
                t->UpdateThreadLocal(num_thread);
                auto func = [&]() {
                    int thread_id = next_thread_id.fetch_add(1);
                    uint64_t start_key = LOAD_SIZE / num_thread * (uint64_t)thread_id;
                    uint64_t end_key = start_key + LOAD_SIZE / num_thread;
                    if(thread_id == num_thread - 1){
                        end_key = LOAD_SIZE;
                    }
                    t->AssignGCID(thread_id);
                    for (uint64_t i = start_key; i < end_key; i++) {
                        t->Insert(init_keys[i], init_keys[i]);
                    }
                    t->UnregisterThread(thread_id);
                };

                std::vector<std::thread> thread_group;

                for (int i = 0; i < num_thread; i++)
                    thread_group.push_back(std::thread{func});

                for (int i = 0; i < num_thread; i++)
                    thread_group[i].join();
                t->UpdateThreadLocal(1);
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::system_clock::now() - starttime);
                printf("Throughput: load, %f ,ops/us\n", (LOAD_SIZE * 1.0) / duration.count());
                if(trial != 0) load_tpts.push_back((LOAD_SIZE * 1.0) / duration.count());
            }

            {
                // Run
                uint64_t err_count = 0;
                auto starttime = std::chrono::system_clock::now();
                next_thread_id.store(0);
                t->UpdateThreadLocal(num_thread);
                auto func = [&]() {
                    std::vector<unsigned long> v{};
                    v.reserve(1);
                    int thread_id = next_thread_id.fetch_add(1);
                    uint64_t start_key = RUN_SIZE / num_thread * (uint64_t)thread_id;
                    uint64_t end_key = start_key + RUN_SIZE / num_thread;
                    if(thread_id == num_thread - 1){
                        end_key = RUN_SIZE;
                    }


                    t->AssignGCID(thread_id);
                    for (uint64_t i = start_key; i < end_key; i++) {
                        if (ops[i] == OP_INSERT) {
                            t->Insert(keys[i], keys[i]);
                        } else if (ops[i] == OP_READ) {
                            v.clear();
                            t->GetValue(keys[i], v);
                            if (v[0] != keys[i]) {
                                err_count++;
                                // std::cout << "[BWTREE] wrong key read: " << v[0] << " expected:" << keys[i] << std::endl;
                            }
                        } else if (ops[i] == OP_SCAN) {
                            auto it = t->Begin(keys[i]);
                            uint64_t key_sum = 0, val_sum = 0;

                            int resultsFound = 0;
                            while (it.IsEnd() != true && resultsFound != ranges[i]) {
                                key_sum += it->first;
                                val_sum += it->second;
                                resultsFound++;
                                it++;
                            }
                            query_results_keys[i] = key_sum;
                            query_results_vals[i] = val_sum;
                        } else if (ops[i] == OP_SCAN_END) {
                            auto it = t->Begin(keys[i]);
                            uint64_t key_sum = 0, val_sum = 0;
                            while (it.IsEnd() != true && it->first != range_end[i]) {
                                key_sum += it->first;
                                val_sum += it->second;
                                it++;
                            }
                            query_results_keys[i] = key_sum;
                            query_results_vals[i] = val_sum;
                        } else if (ops[i] == OP_UPDATE) {
                            std::cout << "NOT SUPPORTED CMD!\n";
                            exit(0);
                        }
                    }
                    t->UnregisterThread(thread_id);
                };

                std::vector<std::thread> thread_group;

                for (int i = 0; i < num_thread; i++)
                    thread_group.push_back(std::thread{func});

                for (int i = 0; i < num_thread; i++)
                    thread_group[i].join();
                t->UpdateThreadLocal(1);
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::system_clock::now() - starttime);
                printf("Throughput: run, %f ,ops/us\n", (RUN_SIZE * 1.0) / duration.count());
                printf("\n errors count = %lu \n",err_count);
                if(trial != 0) run_tpts.push_back((RUN_SIZE * 1.0) / duration.count());

            }
            uint64_t key_sum = 0;
            uint64_t val_sum = 0;
            for(int i = 0; i < RUN_SIZE; i++) {
                key_sum += query_results_keys[i];
                val_sum += query_results_vals[i];
            }
            printf("\ttotal key sum = %lu, total val sum = %lu\n\n", key_sum, val_sum);
        }


        if (index_type == TYPE_BTREEOLC) {
            btreeolc::BTree<uint64_t, uint64_t> t;
            
            std::vector<uint64_t> query_results_keys(RUN_SIZE);
            std::vector<uint64_t> query_results_vals(RUN_SIZE);

            {
                // Load
                auto starttime = std::chrono::system_clock::now();
                parallel_for(0, LOAD_SIZE, [&](const uint64_t &i) {
                    t.insert(init_keys[i], init_keys[i]);
                });
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::system_clock::now() - starttime);
                printf("Throughput: load, %f ,ops/us\n", (LOAD_SIZE * 1.0) / duration.count());
            }

            {
                // Run
                auto starttime = std::chrono::system_clock::now();
                parallel_for(0, RUN_SIZE, [&](const uint64_t &i) {
                    if (ops[i] == OP_INSERT) {
                        t.insert(keys[i], keys[i]);
                    } else if (ops[i] == OP_READ) {
                        uint64_t result;
                        t.lookup(keys[i],result);
                        if( result && result != keys[i]) {
                            std::cout << "[BTREE] wrong key read: " << result << " expected:" << keys[i] << std::endl;
                            exit(0);
                        }
                    } else if (ops[i] == OP_SCAN) {
                        uint64_t start= keys[i];
			            uint64_t key_sum = 0, val_sum = 0;

                        uint64_t results[range];
                        uint64_t count = t.scan(keys[i], ranges[i], results);
                        for (size_t j = 0; j < count; j++) {
                            key_sum += results[j];
                            val_sum += results[j];
                        }

                        query_results_keys[i] = key_sum;
                        query_results_vals[i] = val_sum;
                    } else if (ops[i] == OP_SCAN_END) {
			            std::cout << "NOT SUPPORTED CMD!\n";
                        exit(0);
                    } else if (ops[i] == OP_UPDATE) {
                        std::cout << "NOT SUPPORTED CMD!\n";
                        exit(0);
                    }
                });
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::system_clock::now() - starttime);
                printf("Throughput: run, %f ,ops/us\n", (RUN_SIZE * 1.0) / duration.count());
            }
        }
    }

    printf("Median load : , %f ,ops/us\n", findMedian(load_tpts));
    printf("Median run : , %f ,ops/us\n", findMedian(run_tpts));

    return;
}


/*
 * GetThreadNum() - Returns the number of threads used for multithreaded testing
 *
 * By default 40 threads are used
 */
uint64_t GetThreadNum() {    
    uint64_t thread_num = 40;
    bool ret = Envp::GetValueAsUL("THREAD_NUM", &thread_num);
    if(ret == false) {
      throw "THREAD_NUM must be an unsigned ineteger!"; 
    } else {
      printf("Using thread_num = %lu\n", thread_num); 
    }
    
    return thread_num;
}

int main(int argc, char **argv) {
  bool run_benchmark_all = false;
  bool run_test = false;
  bool run_benchmark_bwtree = false;
  bool run_benchmark_bwtree_full = false;
  bool run_benchmark_btree_full = false;
  bool run_benchmark_art_full = false;
  bool run_stress = false;
  bool run_epoch_test = false;
  bool run_infinite_insert_test = false;
  bool run_email_test = false;
  bool run_mixed_test = false;
  bool run_ycsb_test = false;

  int opt_index = 1;
  // while(opt_index < argc) {
    char *opt_p = argv[opt_index];

    if(strcmp(opt_p, "--benchmark-all") == 0) {
      run_benchmark_all = true;
    } else if(strcmp(opt_p, "--test") == 0) {
      run_test = true;
    } else if(strcmp(opt_p, "--benchmark-bwtree") == 0) {
      run_benchmark_bwtree = true;
    } else if(strcmp(opt_p, "--benchmark-bwtree-full") == 0) {
      run_benchmark_bwtree_full = true;
    } else if(strcmp(opt_p, "--benchmark-btree-full") == 0) {
      run_benchmark_btree_full = true;
    } else if(strcmp(opt_p, "--benchmark-art-full") == 0) {
      run_benchmark_art_full = true;
    } else if(strcmp(opt_p, "--stress-test") == 0) {
      run_stress = true;
    } else if(strcmp(opt_p, "--epoch-test") == 0) {
      run_epoch_test = true;
    } else if(strcmp(opt_p, "--infinite-insert-test") == 0) {
      run_infinite_insert_test = true;
    } else if(strcmp(opt_p, "--email-test") == 0) {
      run_email_test = true;
    } else if(strcmp(opt_p, "--mixed-test") == 0) {
      run_mixed_test = true;
    } else if(strcmp(opt_p, "--ycsb-test") == 0) {
      run_ycsb_test = true;
    } else {
      printf("ERROR: Unknown option: %s\n", opt_p);

      return 0;
    }

  //   opt_index++;
  // }

  // bwt_printf("RUN_BENCHMARK_ALL = %d\n", run_benchmark_all);
  // bwt_printf("RUN_BENCHMARK_BWTREE_FULL = %d\n", run_benchmark_bwtree_full);
  // bwt_printf("RUN_BENCHMARK_BWTREE = %d\n", run_benchmark_bwtree);
  // bwt_printf("RUN_BENCHMARK_ART_FULL = %d\n", run_benchmark_art_full);
  // bwt_printf("RUN_TEST = %d\n", run_test);
  // bwt_printf("RUN_STRESS = %d\n", run_stress);
  // bwt_printf("RUN_EPOCH_TEST = %d\n", run_epoch_test);
  // bwt_printf("RUN_INFINITE_INSERT_TEST = %d\n", run_infinite_insert_test);
  // bwt_printf("RUN_EMAIL_TEST = %d\n", run_email_test);
  // bwt_printf("RUN_MIXED_TEST = %d\n", run_mixed_test);
  bwt_printf("RUN_YCSB_TEST = %d\n", run_ycsb_test);
  bwt_printf("======================================\n");

  //////////////////////////////////////////////////////
  // Next start running test cases
  //////////////////////////////////////////////////////

  TreeType *t1 = nullptr;

  if(run_ycsb_test == true) {
    t1 = GetEmptyTree();

    printf("\n\nStarting ycsb testing...\n");
    if (argc != 7) {
        std::cout << "Usage: ./ycsb [index type] [ycsb workload type] [key distribution] [access pattern] [number of threads]\n";
        std::cout << "1. index type: art hot bwtree masstree clht\n";
        std::cout << "               fastfair levelhash cceh woart\n";
        std::cout << "2. ycsb workload type: a, b, c, e\n";
        std::cout << "3. key distribution: randint, string\n";
        std::cout << "4. access pattern: uniform, zipfian\n";
        std::cout << "5. number of threads (integer)\n";
        return 1;
    }
     printf("%s, workload%s, %s, %s, threads %s\n", argv[2], argv[3], argv[4], argv[5], argv[6]);
    // LaunchParallelTestID(t1, mixed_thread_num, MixedTest1, t1);

    int index_type;
    if (strcmp(argv[2], "art") == 0)
        index_type = TYPE_ART;

    else if (strcmp(argv[2], "btree") == 0)
        index_type = TYPE_BTREEOLC;
    else if (strcmp(argv[2], "hot") == 0) {
#ifdef HOT
        index_type = TYPE_HOT;
#else
        return 1;
#endif
    } else if (strcmp(argv[2], "bwtree") == 0)
        index_type = TYPE_BWTREE;
    else if (strcmp(argv[2], "masstree") == 0)
        index_type = TYPE_MASSTREE;
    else if (strcmp(argv[2], "clht") == 0)
        index_type = TYPE_CLHT;
    else if (strcmp(argv[2], "fastfair") == 0)
        index_type = TYPE_FASTFAIR;
    else if (strcmp(argv[2], "levelhash") == 0)
        index_type = TYPE_LEVELHASH;
    else if (strcmp(argv[2], "cceh") == 0)
        index_type = TYPE_CCEH;
    else if (strcmp(argv[2], "woart") == 0)
        index_type = TYPE_WOART;
    else {
        fprintf(stderr, "Unknown index type: %s\n", argv[2]);
        exit(1);
    }

    int wl;
    if (strcmp(argv[3], "a") == 0) {
        wl = WORKLOAD_A;
    } else if (strcmp(argv[3], "b") == 0) {
        wl = WORKLOAD_B;
    } else if (strcmp(argv[3], "c") == 0) {
        wl = WORKLOAD_C;
    } else if (strcmp(argv[3], "d") == 0) {
        wl = WORKLOAD_D;
    } else if (strcmp(argv[3], "e") == 0) {
        wl = WORKLOAD_E;
    } else if (strcmp(argv[3], "x") == 0) {
        wl = WORKLOAD_X;
    } else if (strcmp(argv[3], "y") == 0) {
        wl = WORKLOAD_Y;
    } else {
        fprintf(stderr, "Unknown workload: %s\n", argv[3]);
        exit(1);
    }

    int kt;
    if (strcmp(argv[4], "randint") == 0) {
        kt = RANDINT_KEY;
    } else if (strcmp(argv[4], "string") == 0) {
        kt = STRING_KEY;
    } else {
        fprintf(stderr, "Unknown key type: %s\n", argv[4]);
        exit(1);
    }

    int ap;
    if (strcmp(argv[5], "uniform") == 0) {
        ap = UNIFORM;
    } else if (strcmp(argv[5], "zipfian") == 0) {
        ap = ZIPFIAN;
    } else {
        fprintf(stderr, "Unknown access pattern: %s\n", argv[5]);
        exit(1);
    }

    int num_thread = atoi(argv[6]);
    // tbb::task_scheduler_init init(num_thread);

    if (kt != STRING_KEY) {
        std::vector<uint64_t> init_keys;
        std::vector<uint64_t> keys;
        std::vector<uint64_t> ranges_end;
        std::vector<int> ranges;
        std::vector<int> ops;

        init_keys.reserve(LOAD_SIZE);
        keys.reserve(RUN_SIZE);
        ranges_end.reserve(RUN_SIZE);
        ranges.reserve(RUN_SIZE);
        ops.reserve(RUN_SIZE);

        memset(&init_keys[0], 0x00, LOAD_SIZE * sizeof(uint64_t));
        memset(&keys[0], 0x00, RUN_SIZE * sizeof(uint64_t));
        memset(&ranges_end[0], 0x00, RUN_SIZE * sizeof(uint64_t));
        memset(&ranges[0], 0x00, RUN_SIZE * sizeof(int));
        memset(&ops[0], 0x00, RUN_SIZE * sizeof(int));

        ycsb_load_run_randint(index_type, wl, kt, ap, num_thread, init_keys, keys,ranges_end, ranges, ops);
    } 

    printf("Finished ycsb testing\n\n");

    // PrintStat(t1);

    // MixedGetValueTest(t1);

    // DestroyTree(t1);
  }
}