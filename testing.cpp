#include <functional>
#include "mpcops.hpp"
#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include "rdpf.hpp"
#include "shapes.hpp"
#include "testing.hpp"

// #define SEGTREE_VERBOSE
// To enable timing/stat instrumentation define SEGTREE_VERBOSE2
#define SEGTREE_VERBOSE2

#ifdef SEGTREE_VERBOSE2
#define STATS_PRE() do { tio.sync_lamport(); mpcio.reset_stats(); tio.reset_lamport(); } while(0)
#define STATS_POST(MSG) do { tio.sync_lamport(); std::cout << MSG << std::endl; mpcio.dump_stats(std::cout); mpcio.reset_stats(); tio.reset_lamport(); } while(0)
#else
#define STATS_PRE() do {} while(0)
#define STATS_POST(MSG) do {} while(0)
#endif

int count = 10; // Number of times to access in test functions to get a better average

// size_t arrayIndexToLevelPos(size_t idx) {
//     // Compute level as the floor of log2(idx)
//     size_t level = static_cast<size_t>(std::log2(idx));
//     // For a complete binary tree with root at index 1,
//     // leftmost index in this level is 2^level.
//     size_t pos = idx - (1ULL << level);

//     return pos;
// }

// void SegmentTree::init(MPCTIO &tio, yield_t & yield) {
//     auto SegTreeArray = oram.flat(tio, yield);
//     size_t num_leaves = 1ULL << (depth - 1);
//     size_t leaf_start = num_leaves;

//     // Create and initialize segment tree array
//     std::vector<uint64_t> segTree(num_items, 0);

//     // Initialize leaf nodes with values 0, 100, 200, ..., (num_leaves-1)*100
//     for (size_t i = 0; i < num_leaves; ++i) {
//         segTree[leaf_start + i] = i * 100;
//     }

//     // Build internal nodes bottom-up by computing range sums
//     for (int i = num_leaves - 1; i >= 1; --i) {
//         segTree[i] = segTree[2*i] + segTree[2*i + 1];
//     }

//     SegTreeArray.init([this, &segTree] (size_t i) -> size_t {
//         if (i >= 1 && i < num_items) {
//             return segTree[i];
//         } else {
//             return size_t(0x7fffffffffffffff);
//         }
//     });

//    auto isEvenArray = isEven.flat(tio, yield);
//    isEvenArray.init([this] (size_t i) -> size_t {
//      if (i >= 1 && i < num_items) {
//           return (i % 2 == 0) ? 1 : 0;
//      } else {
//           return size_t(1); // 0 is even, just for consistency
//      }
//     });

//     auto siblingArray = sibling.flat(tio, yield);
//     siblingArray.init([this] (size_t i) -> size_t {
//         if (i >= 1 && i < num_items) {
//             return arrayIndexToLevelPos((i % 2 == 0) ? i + 1 : i - 1); // we took 1s sibling as 0 and 0s as 1 which should not cause any issue.
//         } else {
//             return size_t(0);
//         }
//     });

//     auto parentArray = parent.flat(tio, yield);
//     parentArray.init([this] (size_t i) -> size_t {
//         if (i >= 1 && i < num_items) {
//             return arrayIndexToLevelPos(i/2);
//         } else {
//             return size_t(0);
//         }
//     });
// }

// Original Update function with detailed timing per operation

// Test function to isolate Flat object creation overhead
void SegmentTreeTest::Test_FlatCreation(MPCTIO &tio, MPCIO &mpcio, yield_t & yield) {
    auto SegTreeArray = oram.flat(tio, yield);
    auto parentArray = parent.flat(tio, yield);
    
    // Just create Flat objects for each level without doing any operations
    for(size_t i=1; i<=depth; i++) {
        size_t level = depth - i;
        typename Duoram < RegAS > ::Flat parentLevel(parentArray, tio, yield, (1ULL << level), (1ULL << level)+1);
        typename Duoram < RegAS > ::Flat segTreeLevel(SegTreeArray, tio, yield, (1ULL << level), (1ULL << level)+1);
        // No operations - just creation overhead
    }
}

// Test 1: Fresh local ORAM with depth number of RegAS accesses
void SegmentTreeTest::Test_LocalORAM(MPCTIO &tio, MPCIO &mpcio, yield_t & yield) {
    // Create a fresh local ORAM of same size
    Duoram<RegAS> localarray(tio.player(), num_items);
    auto localArray = localarray.flat(tio, yield);
    
    // Access using RegAS shares depth number of times
    RegAS index;  
    RegAS one;
    one.set(tio.player() == 0 ? 1 : 0);  

    for(size_t i = 0; i < depth; i++) {
        auto startTime = std::chrono::high_resolution_clock::now();

        // access it count times to get a better average
        for (size_t j = 0; j < count; j++) {
            RegAS value = localArray[index];
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::micro> duration = endTime - startTime;
        std::cout << "Access " << i << " time: " << duration.count() / 1000 << " ms" << std::endl;
    }
}

// Test 2: Use existing SegTreeArray depth number of times
void SegmentTreeTest::Test_SegTreeArray(MPCTIO &tio, MPCIO &mpcio, yield_t & yield) {
    auto SegTreeArray = oram.flat(tio, yield);
    
    RegAS index; // 0th index
    RegAS one;
    one.set(tio.player() == 0 ? 1 : 0);  

    // Access SegTreeArray depth number of times
    for(size_t i = 0; i < depth; i++) {
        auto startTime = std::chrono::high_resolution_clock::now();
        
        // access it count times to get a better average
        for (size_t j = 0; j < count; j++) {
            RegAS value = SegTreeArray[index];
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::micro> duration = endTime - startTime;
        std::cout << "Access " << i << " time: " << duration.count() / 1000 << " ms" << std::endl;
    }
}

// Test 3: Level-wise Flat objects with accesses (like segment tree pattern)
void SegmentTreeTest::Test_LevelWiseAccess(MPCTIO &tio, MPCIO &mpcio, yield_t & yield) {
    auto SegTreeArray = oram.flat(tio, yield);
    
    RegAS index;
    RegAS one;
    one.set(tio.player() == 0 ? 1 : 0);  

    // For each level, create Flat object and access
    for(size_t i = 1; i <= depth; i++) {
        size_t level = depth - i;
        auto timebeforeCreatingLevelArray = std::chrono::high_resolution_clock::now();
        typename Duoram<RegAS>::Flat segTreeLevel(SegTreeArray, tio, yield, (1ULL << level), (1ULL << level)+1);
        auto timeafterCreatingLevelArray = std::chrono::high_resolution_clock::now();

        // Access the level count times
        for (size_t j = 0; j < count; j++) {
            RegAS value = segTreeLevel[index];
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::micro> constructionDuration = timeafterCreatingLevelArray - timebeforeCreatingLevelArray;
        std::chrono::duration<double, std::micro> accessDuration = endTime - timeafterCreatingLevelArray;
        std::cout << "Level " << level << " Flat construction time: " << constructionDuration.count() / 1000 << " ms" << " access time: " << accessDuration.count() / 1000 << " ms" << std::endl;
    }
}

void SegTreeTest(MPCIO &mpcio, const PRACOptions &opts, char **args) {
    // Parse command line arguments
    int nargs = 0;
    while (args[nargs] != nullptr) {
        ++nargs;
    }

    nbits_t depth = 5;
    size_t n_updates = 1;
    size_t n_queries = 1;

    for (int i = 0; i < nargs; i += 2) {
        std::string option = args[i];
        if (option == "-d" && i + 1 < nargs) {
            depth = std::atoi(args[i + 1]);
        } else if (option == "-u" && i + 1 < nargs) {
            n_updates = std::atoi(args[i + 1]);
        } else if (option == "-q" && i + 1 < nargs) {
            n_queries = std::atoi(args[i + 1]);
        }
        else if (option == "-c" && i + 1 < nargs) {
            count = std::atoi(args[i + 1]);
        }
    }

    address_t len = (1<<depth);

    MPCTIO tio(mpcio, 0, opts.num_cpu_threads);

    run_coroutines(tio, [&tio, &mpcio, len, depth, n_updates, n_queries] (yield_t &yield) {
        
        std::cout << "===== SegmentTree Initialization Test =====" << std::endl;        
        // Create SegmentTree object
        STATS_PRE();
        SegmentTreeTest segTree(tio.player(), len, depth);
        STATS_POST("[SEGTREE] Object creation");
        
        // Initialize the segment tree
        // STATS_PRE();
        // segTree.init(tio, yield);
        // STATS_POST("[SEGTREE] Initialization");
        
        std::cout << "===== Initialization Complete =====" << std::endl;
        
        // Test both update methods for comparison
        RegAS index, value;
        index.set(tio.player() == 0 ? 5 : 0);  // Update index 5
        value.set(tio.player() == 0 ? 1000 : 0);  // New value 1000
        
        std::cout << "\n===== Test 1: Local ORAM Access (" << depth << " accesses) =====" << std::endl;
        STATS_PRE();
        segTree.Test_LocalORAM(tio, mpcio, yield);
        STATS_POST("[SEGTREE] Local ORAM Access Test");
        
        std::cout << "\n===== Test 2: SegTreeArray Direct Access (" << depth << " accesses) =====" << std::endl;
        STATS_PRE();
        segTree.Test_SegTreeArray(tio, mpcio, yield);
        STATS_POST("[SEGTREE] SegTreeArray Direct Access Test");
        
        std::cout << "\n===== Test 3: Level-wise Flat Access Pattern =====" << std::endl;
        STATS_PRE();
        segTree.Test_LevelWiseAccess(tio, mpcio, yield);
        STATS_POST("[SEGTREE] Level-wise Flat Access Test");
        
        std::cout << "\n===== Performance Comparison Complete =====" << std::endl;
    });
}
