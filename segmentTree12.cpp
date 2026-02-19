#include "segmentTree12.hpp"
#include "cell.hpp"
#include "coroutine.hpp"
#include "duoram.hpp"
#include "logger.hpp"
#include "mpcops.hpp"
#include "rdpf.hpp"
#include "shapes.hpp"
#include "types.hpp"
#include <chrono>
#include <functional>

// uncomment to enable intermediate reconstructions and logging
// #define SEGTREE_VERBOSE

/*
The segment tree data structure in SegmentTree12 uses an ARRAY layout (heap-style indexing).

ARRAY STRUCTURE:
- Complete binary tree stored in an array with root at index 1
- Each level k has 2^k nodes
- Total nodes in the tree array: 2^d (indices [1 .. 2^d - 1] are used)
- Level k starts at array index: 2^k

LEVEL-WISE LAYOUT EXAMPLE:
For depth d=3, array [100, 200, 300, 400] (4 elements):

Level 0: [600]           (Root level: root_sum)
Level 1: [300][700]      (Level 1: left_sum + right_sum)
Level 2: [100][200][300][400]  (Leaf level: leaf_values)

Array storage:
Index:  [1][2][3][4][5][6][7]
Value:  [600][300][700][100][200][300][400]
Level:   0   1   1   2   2   2  2

MATHEMATICAL FORMULAS:
- Level k start index: 2^k
- Level k length: 2^k
- Total nodes (including index 0 unused): 2^d
*/

/*
RANGE SUM QUERY OPERATION (Level-wise Processing)

The range sum query computes the sum of elements in a given range [left, right]
using a level-wise segment tree structure. The algorithm processes levels using
level-specific Flat objects for performance in MPC.

Indexing notes:
- Parent of node at array index i is at i/2 (i >> 1). In this implementation,
  indices are maintained as level-relative indices while traversing levels, so
  shifting right still yields the parent index at the next higher level.
- Sibling indices are obtained via a separate Sibling ORAM (SiblingOram) that
  maps the current node's level-relative index to its sibling's level-relative
  index. This is the key distinction of SegmentTree12 versus variants that
  compute sibling indices via MPC directly.

Algorithm Overview:
1. Start with leaf-level indices (left, right)
2. Precompute the path indices for all levels going up the tree
3. For each level from leaves to root:
   - Create level-specific Flat objects for the Tree ORAM and the Sibling ORAM
   - Compute Additive Shares of XOR shares of left and right indices
   - Use additive-share comparison to compute isNotDone
   - Use the least significant bit of XOR shares of left/right indices as a
     boolean-share isOdd
   - Read left/right sibling indices from the Sibling ORAM; then read the sibling
     values from the Tree ORAM
   - Add left sibling value if left is a left child: (!isOdd && isNotDone)
   - Add right sibling value if right is a right child: (isOdd && isNotDone)
*/

/*
UPDATE OPERATION (Level-wise Propagation)

The update operation modifies a single array element and propagates the change
up through all levels using level-specific addressing and parent relationships.

Algorithm Overview:
1. Calculate difference between new and current values
2. For each level from leaves to root:
   - Create level-specific Flat objects
   - Update current node with difference
   - Move to parent using level-relative parent index
*/

struct LevelStats {
    double start_ms;
    double indices_ms;
    double is_not_done_ms;
    double siblings_start_ms;
    double left_oram_access_ms;
    double right_oram_access_ms;
    double siblings_ms;
    double end_ms;
};

// NOTE:
// SegmentTree12 stores sibling pointers directly as absolute tree indices
// (i.e., the sibling of node i is i^1 / i±1 in heap indexing), so we no longer
// need arrayIndexToLevelIndex12().

void SegmentTree12::init(MPCTIO &tio, yield_t &yield) {

    // Create a flat reference to the main segment tree ORAM
    auto SegTreeArray = TreeOram.flat(tio, yield);
    size_t num_leaves = 1ULL << (depth - 1);
    size_t leaf_start = num_leaves;

    // Create and initialize a local segment tree array
    std::vector<uint64_t> segTree(num_items, 0);

    // Initialize leaf nodes with values 0, 100, 200, ..., (num_leaves-1)*100
    for (size_t i = 0; i < num_leaves; ++i) {
        segTree[leaf_start + i] = i * 100;
    }

    // Build internal nodes bottom-up by computing range sums
    for (int i = num_leaves - 1; i >= 1; --i) {
        segTree[i] = segTree[2 * i] + segTree[2 * i + 1];
    }

    // Initialize the ORAM segment tree with computed values
    SegTreeArray.init([this, &segTree](size_t i) -> size_t {
        if (i >= 1 && i < num_items) {
            return segTree[i];
        } else {
            return size_t(0x7fffffffffffffff);
        }
    });

    // Initialize sibling ORAM:
    // Store absolute sibling tree index (heap layout).
    // For i in [1 .. num_items-1], sibling is (i ^ 1).
    // (This assumes root is at 1 and index 0 is unused.)
    auto SiblingArray = SiblingOram.flat(tio, yield);
    SiblingArray.init([this](size_t i) -> size_t {
        if (i >= 1 && i < num_items) {
            return (i ^ 1ULL);
        } else {
            return size_t(0);
        }
    });
}

// helper function to reconstruct and print SegTreeArray
void SegmentTree12::printSegmentTree(MPCTIO &tio, yield_t &yield) {
    auto SegTreeArray = TreeOram.flat(tio, yield);
    auto SegTreeRecons = SegTreeArray.reconstruct();
    for (size_t i = 1; i < num_items; i++) {
        std::cout << "SegTreeArray[" << i << "] = " << SegTreeRecons[i].share() << std::endl;
    }
}

// Main function to compute range sum directly (optimized - no intermediate
// bitVec)
RegAS SegmentTree12::computeRangeSum(MPCTIO &tio, MPCIO &mpcio, yield_t &yield, RegXS leftIndex,
                                    RegXS rightIndex, PerformanceLogger *logger, size_t operation_id) {

    auto SegTreeArray = TreeOram.flat(tio, yield);
    auto SiblingArray = SiblingOram.flat(tio, yield);

    // Phase 1: Pre-compute the path indices for all levels
    auto start_time = std::chrono::high_resolution_clock::now();
    // Store absolute tree indices for each level (root-to-leaf layout in a heap array)
    std::vector<RegXS> leftPathIndex(depth);
    std::vector<RegXS> rightPathIndex(depth);

    // Compute parent indices for all levels going up the tree.
    // Inputs are absolute indices in the flattened tree (heap indexing), so parent is (idx >> 1).
    for (uint32_t i = 1; i <= depth; i++) {
        size_t level = depth - i; // Current level (leaf to root)

        leftPathIndex[level] = leftIndex;
        rightPathIndex[level] = rightIndex;

        leftIndex >>= 1;
        rightIndex >>= 1;
    }

    auto step1_end_time = std::chrono::high_resolution_clock::now();
    auto step1_duration = std::chrono::duration_cast<std::chrono::milliseconds>(step1_end_time - start_time).count();
    std::cout << "Step 1 (Path Computation) Time: " << step1_duration << " ms" << std::endl;

    // Phase 2: Compute sum directly while identifying nodes to include
    std::vector<RegAS> levelSums(depth);
    for (size_t level = 0; level < depth; level++) {
        levelSums[level].set(0);
    }

    std::vector<LevelStats> stats(depth);
    auto parallel_start_time = std::chrono::high_resolution_clock::now();

    std::vector<coro_t> sum_coroutines;
    for (uint32_t i = 1; i <= depth; i++) {
        size_t level = depth - i;

        sum_coroutines.emplace_back([&tio, level, i, this, &leftPathIndex, &rightPathIndex, &SegTreeArray,
                                     &SiblingArray, &levelSums, &stats, parallel_start_time](yield_t &sub_yield) {
            auto level_start = std::chrono::high_resolution_clock::now();
            stats[level].start_ms =
                std::chrono::duration<double, std::milli>(level_start - parallel_start_time).count();

            // size_t levelStart = 1ULL << level;
            // size_t levelLength = (1ULL << level) + 1;

            // Create one, zero which are required in mpc operations further
            RegAS one;
            one.set(tio.player());

            RegAS zero;
            zero.set(0);

            // NOTE:
            // We intentionally do NOT create level-wise Flat sub-ORAMs here.
            // All accesses use the main ORAM flats directly (SegTreeArray, SiblingArray).
            // (void)levelStart;
            // (void)levelLength;

            // Absolute tree indices for this level
            RegXS leftIndex = leftPathIndex[level];
            RegXS rightIndex = rightPathIndex[level];

            RegAS leftIndexAS;
            RegAS rightIndexAS;
            run_coroutines(
                sub_yield,
                [&tio, leftIndex, &leftIndexAS](yield_t &yield) { mpc_xs_to_as(tio, yield, leftIndexAS, leftIndex); },
                [&tio, rightIndex, &rightIndexAS](yield_t &yield) {
                    mpc_xs_to_as(tio, yield, rightIndexAS, rightIndex);
                });

            stats[level].indices_ms = std::chrono::duration<double, std::milli>(
                                          std::chrono::high_resolution_clock::now() - parallel_start_time)
                                          .count();

            // if left and right are adjacent or same then the marking of nodes is
            // already completed so we set isDone here
            RegBS isNotDone;
            CDPF cdpf1 = tio.cdpf(sub_yield);
            RegAS diff1 = rightIndexAS - (leftIndexAS + one); // diff1 = right - (left + 1)
            auto [lt_c1, eq_c1, gt_c1] = cdpf1.compare(tio, sub_yield, diff1, tio.aes_ops());

            stats[level].is_not_done_ms = std::chrono::duration<double, std::milli>(
                                              std::chrono::high_resolution_clock::now() - parallel_start_time)
                                              .count();

            isNotDone = gt_c1; // isNotDone = (right > left + 1)

            stats[level].siblings_start_ms = std::chrono::duration<double, std::milli>(
                                                 std::chrono::high_resolution_clock::now() - parallel_start_time)
                                                 .count();

            // include left and right independently based on isValid and isDifferent
            // and if left is left child and right is right child if xor of last bit
            // is 1 then its odd -> right child else left child if leftLastBit is 1
            // then left is right child else left child
            RegBS leftLastBit = leftIndex.bitat(0);
            RegBS rightLastBit = rightIndex.bitat(0);

            RegAS leftSum, rightSum;
            run_coroutines(
                sub_yield,
                [&tio, &SegTreeArray, &SiblingArray, isNotDone, leftIndex, &leftSum, leftLastBit, &stats, level,
                 parallel_start_time](yield_t &yield) {
                    auto SiblingArrayCoro = SiblingArray.context(yield);
                    auto SegTreeArrayCoro = SegTreeArray.context(yield);

                    // Indices are already absolute in the flattened tree
                    RegAS leftSibIndex = SiblingArrayCoro[leftIndex];
                    RegAS leftSibValue = SegTreeArrayCoro[leftSibIndex];

                    stats[level].left_oram_access_ms =
                        std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() -
                                                                  parallel_start_time)
                            .count();

                    RegBS isLeftIncluded;
                    mpc_and(tio, yield, isLeftIncluded, isNotDone, leftLastBit ^ tio.player());
                    mpc_flagmult(tio, yield, leftSum, isLeftIncluded, leftSibValue);
                },
                [&tio, &SegTreeArray, &SiblingArray, isNotDone, rightIndex, &rightSum, rightLastBit, &stats,
                 level, parallel_start_time](yield_t &yield) {
                    auto SiblingArrayCoro = SiblingArray.context(yield);
                    auto SegTreeArrayCoro = SegTreeArray.context(yield);

                    // Indices are already absolute in the flattened tree
                    RegAS rightSibIndex = SiblingArrayCoro[rightIndex];
                    RegAS rightSibValue = SegTreeArrayCoro[rightSibIndex];

                    stats[level].right_oram_access_ms =
                        std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() -
                                                                  parallel_start_time)
                            .count();

                    RegBS isRightIncluded;
                    mpc_and(tio, yield, isRightIncluded, isNotDone, rightLastBit);
                    mpc_flagmult(tio, yield, rightSum, isRightIncluded, rightSibValue);
                });

            stats[level].siblings_ms = std::chrono::duration<double, std::milli>(
                                           std::chrono::high_resolution_clock::now() - parallel_start_time)
                                           .count();

            levelSums[level] += leftSum + rightSum;

            // if its the leaf level then we also have to include the leftIndex and
            // rightIndex themselves after checking validation and duplication
            if (i == 1) {
                CDPF cdpf2 = tio.cdpf(sub_yield);
                RegAS diff2 = rightIndexAS - leftIndexAS; // diff1 = right - left
                auto [lt_c2, eq_c2, gt_c2] = cdpf2.compare(tio, sub_yield, diff2, tio.aes_ops());

                RegAS leftLeafContribution;
                RegAS rightLeafContribution;
                run_coroutines(
                    sub_yield,
                    [&tio, &lt_c2, &eq_c2, &gt_c2, &zero, &one, &leftIndex, &SegTreeArray,
                     &leftLeafContribution](yield_t &yield) {
                        // if right >= left add left
                        // that is equivalent to !(left > right) => lt_c2 ^ tio.player()
                        auto SegTreeArrayCoro = SegTreeArray.context(yield);
                        RegAS leftLeafValue = SegTreeArrayCoro[leftIndex];
                        mpc_select(tio, yield, leftLeafContribution, lt_c2 ^ tio.player(), zero, leftLeafValue);
                    },
                    [&tio, &lt_c2, &eq_c2, &gt_c2, &zero, &one, &rightIndex, &SegTreeArray,
                     &rightLeafContribution](yield_t &yield) {
                        // Only add right leaf if it's different from left (to avoid
                        // double counting) that is if right > left
                        auto SegTreeArrayCoro = SegTreeArray.context(yield);
                        RegAS rightLeafValue = SegTreeArrayCoro[rightIndex];
                        mpc_select(tio, yield, rightLeafContribution, gt_c2, zero, rightLeafValue);
                    });

                levelSums[level] += leftLeafContribution + rightLeafContribution;
            }

            auto level_end = std::chrono::high_resolution_clock::now();
            stats[level].end_ms = std::chrono::duration<double, std::milli>(level_end - parallel_start_time).count();
        });
    }

    run_coroutines(tio, sum_coroutines);

    // Combine all level sums
    RegAS totalSum;
    totalSum.set(0);
    for (size_t level = 0; level < depth; level++) {
        totalSum += levelSums[level];
    }

    auto step2_end_time = std::chrono::high_resolution_clock::now();
    auto step2_duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(step2_end_time - step1_end_time).count();

    std::cout << "Step 2 (Direct Sum Computation) Time: " << step2_duration << " ms" << std::endl;

    // Log detailed level-wise stats
    // if (logger) {
    //     std::ostringstream oss;
    //     oss << "\nLevel-wise Timing in ms (relative to parallel start) "
    //            "[T=Timestamp, "
    //            "D=Duration]:\n";
    //     oss << "Level | Start(T)   | Indices(T) | IsDone(T)  | SibStart(T) | LeftORAM(T) | RightORAM(T) | Siblings(T) "
    //            "| "
    //            "End(T)     | Total(D)\n";
    //     oss << "------+------------+------------+------------+-------------+-------------+--------------+-------------+"
    //            "-------"
    //            "-----+----------\n";
    //     for (size_t level = 0; level < depth; level++) {
    //         if (stats[level].end_ms > 0) {
    //             double t_start = stats[level].start_ms;
    //             double t_idx = stats[level].indices_ms;
    //             double t_done = stats[level].is_not_done_ms;
    //             double t_sib_start = stats[level].siblings_start_ms;
    //             double t_left_oram = stats[level].left_oram_access_ms;
    //             double t_right_oram = stats[level].right_oram_access_ms;
    //             double t_sib = stats[level].siblings_ms;
    //             double t_end = stats[level].end_ms;

    //             oss << std::setw(5) << level << " | " << std::setw(10) << std::fixed << std::setprecision(3) << t_start
    //                 << " | " << std::setw(10) << (t_idx) << " | " << std::setw(10) << (t_done) << " | " << std::setw(11)
    //                 << t_sib_start << " | " << std::setw(11) << t_left_oram << " | " << std::setw(12) << t_right_oram
    //                 << " | " << std::setw(11) << (t_sib) << " | " << std::setw(10) << t_end << " | " << std::setw(8)
    //                 << (t_end - t_start) << "\n";

    //             // oss << std::setw(5) << level << " | " << std::setw(10) << std::fixed << std::setprecision(3) <<
    //             // t_start
    //             //     << " | " << std::setw(10) << t_idx << " | " << std::setw(10) << t_done << " | " << std::setw(11)
    //             //     << t_sib << " | " << std::setw(10) << t_end << " | " << std::setw(8) << (t_end - t_start) <<
    //             //     "\n";
    //         }
    //     }
    //     logger->log_output(oss.str());
    // }

    // Log metrics if logger is provided
    if (logger) {
        logger->log_metric("query", operation_id, "path_computation_time", step1_duration);
        logger->log_metric("query", operation_id, "direct_sum_time", step2_duration);
    }

    return totalSum;
}

// Main RangeSum function
void SegmentTree12::RangeSum(MPCTIO &tio, MPCIO &mpcio, yield_t &yield, RegXS left, RegXS right,
                            PerformanceLogger *logger, size_t operation_id) {
    RegAS sum = computeRangeSum(tio, mpcio, yield, left, right, logger, operation_id);

#ifdef SEGTREE_VERBOSE
    value_t answer = mpc_reconstruct(tio, yield, sum);
    std::cout << "Sum = " << answer << std::endl;
#endif
}

// Main Update function
void SegmentTree12::Update(MPCTIO &tio, MPCIO &mpcio, yield_t &yield, RegXS index, RegAS value,
                          PerformanceLogger *logger, size_t operation_id) {

    auto SegTreeArray = TreeOram.flat(tio, yield);

    // Index is already an absolute tree index in the flattened tree
    RegAS currVal = SegTreeArray[index];
    RegAS diff = value - currVal;

// Debugging and intermediate reconstructions
#ifdef SEGTREE_VERBOSE
    auto recons_index = mpc_reconstruct(tio, yield, index);
    auto recons_currVal = mpc_reconstruct(tio, yield, currVal);
    auto recons_newvalue = mpc_reconstruct(tio, yield, value);
    auto recons_diff = mpc_reconstruct(tio, yield, diff);
    std::cout << "Absolute Index to be updated in the leaf level = " << recons_index << std::endl;
    std::cout << "Current Value at index = " << recons_currVal << std::endl;
    std::cout << "New Value to be updated = " << recons_newvalue << std::endl;
    std::cout << "Diff = " << recons_diff << std::endl;
#endif

    auto update_start = std::chrono::high_resolution_clock::now();
    // Phase 1: Pre-compute the path indices for all levels (similar to
    // getBitVector)
    std::vector<RegXS> updatePathIndex(depth);

    // Compute parent indices for all levels going up the tree
    for (size_t i = 1; i <= depth; i++) {
        size_t level = depth - i; // Current level (leaf to root)

        updatePathIndex[level] = index;

        index >>= 1;
    }

    auto phase1_end = std::chrono::high_resolution_clock::now();
    auto phase1_duration = std::chrono::duration_cast<std::chrono::milliseconds>(phase1_end - update_start).count();
    std::cout << "Update Phase 1 (Path Computation) Time: " << phase1_duration << " ms" << std::endl;

    // Phase 2: Parallelize the updates across all levels
    std::vector<coro_t> update_coroutines;
    for (size_t i = 1; i <= depth; i++) {
        size_t level = depth - i;

        update_coroutines.emplace_back([&tio, level, this, &diff, &updatePathIndex, &SegTreeArray](yield_t &sub_yield) {
            // IMPORTANT:
            // When performing ORAM accesses inside a coroutine, always use the coroutine's yield
            // to create the context, and keep all accesses within that same context.
            auto SegTreeArrayCoro = SegTreeArray.context(sub_yield);

            RegXS absIndex = updatePathIndex[level];
            SegTreeArrayCoro[absIndex] += diff;

#ifdef SEGTREE_VERBOSE
            auto recons_Index = mpc_reconstruct(tio, sub_yield, absIndex);
            RegAS updated_val = SegTreeArrayCoro[absIndex];
            auto recons_updated = mpc_reconstruct(tio, sub_yield, updated_val);
            std::cout << "Updated Index = " << (recons_Index) << " with value = " << recons_updated << std::endl;
#endif
        });
    }

    run_coroutines(tio, update_coroutines);

    auto phase2_end = std::chrono::high_resolution_clock::now();
    auto phase2_duration = std::chrono::duration_cast<std::chrono::milliseconds>(phase2_end - phase1_end).count();

    std::cout << "Update Phase 2 (Parallel Updates) Time: " << phase2_duration << " ms" << std::endl;

    // Log metrics if logger is provided
    if (logger) {
        logger->log_metric("update", operation_id, "path_computation_time", phase1_duration);
        logger->log_metric("update", operation_id, "parallel_updates_time", phase2_duration);
    }
}

// Main function to run Segment Tree 8 operations
void SegTree12(MPCIO &mpcio, const PRACOptions &opts, char **args) {
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
    }

    // Calculate new array size: 2^d
    address_t len = (1 << depth);

    MPCTIO tio(mpcio, 0, opts.num_cpu_threads);

    // Create experiment ID based on timestamp and parameters
    auto now = std::time(nullptr);
    auto tm = *std::localtime(&now);
    std::ostringstream exp_id_stream;
    exp_id_stream << "st12_d" << (int)depth << "_u" << n_updates << "_q" << n_queries << "_"
                  << std::put_time(&tm, "%Y%m%d_%H%M%S");
    std::string experiment_id = exp_id_stream.str();

    // Initialize logger (only player 0 will actually create log files)
    PerformanceLogger logger(experiment_id, depth, n_updates, n_queries, tio.player(), tio.getMode());

    run_coroutines(tio, [&tio, &mpcio, len, depth, n_updates, n_queries, &logger](yield_t &yield) {
        // Time initialization
        auto init_start = std::chrono::high_resolution_clock::now();

        SegmentTree12 segTree(tio.player(), len, depth);
        segTree.init(tio, yield);

        auto init_end = std::chrono::high_resolution_clock::now();
        auto init_duration = std::chrono::duration_cast<std::chrono::milliseconds>(init_end - init_start);

        logger.log_section("Segment Tree Init Stats");
        std::ostringstream oss;
        oss << "Updates: " << n_updates << ", Queries: " << n_queries << "\n";
        oss << "Initialization Time: " << init_duration.count() / 1000.0 << " seconds\n";
        logger.log_output(oss.str());

        logger.log_metric("init", 0, "initialization_time", init_duration.count());

        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        mpcio.reset_stats();
        tio.reset_lamport();

// Print initial segment tree
#ifdef SEGTREE_VERBOSE
        logger.log_section("Initial Segment Tree");
        segTree.printSegmentTree(tio, yield);
#endif

        // Time updates
        auto update_start = std::chrono::high_resolution_clock::now();

        // Perform updates
        for (size_t u = 0; u < n_updates; ++u) {
            logger.log_section("Update " + std::to_string(u + 1) + " begins");

            auto single_update_start = std::chrono::high_resolution_clock::now();

            RegXS index;
            size_t leafStart = (1ULL << (depth - 1));
            size_t idx_val = leafStart + (u % (1 << (depth - 1)));
            index.set(tio.player() == 0 ? idx_val : 0);

            RegAS value;
            size_t val_to_set = (u + 1) * 50;
            value.set(tio.player() == 0 ? val_to_set : 0);

            segTree.Update(tio, mpcio, yield, index, value, &logger, u + 1);

            auto single_update_end = std::chrono::high_resolution_clock::now();
            auto single_update_duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(single_update_end - single_update_start).count();

            logger.log_metric("update", u + 1, "total_update_time", single_update_duration);
            logger.log_output("Update " + std::to_string(u + 1) + " ends\n");
        }

        auto update_end = std::chrono::high_resolution_clock::now();
        auto update_duration = std::chrono::duration_cast<std::chrono::milliseconds>(update_end - update_start);

        logger.log_section("Updates Stats");
        oss.str("");
        oss << "Total Updates Time: " << update_duration.count() / 1000.0 << " seconds\n";
        logger.log_output(oss.str());

        logger.log_metric("update_summary", 0, "total_all_updates_time", update_duration.count());

        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        mpcio.reset_stats();
        tio.reset_lamport();

// Print updated segment tree
#ifdef SEGTREE_VERBOSE
        logger.log_section("Updated Segment Tree");
        segTree.printSegmentTree(tio, yield);
#endif

        // Time range queries
        auto query_start = std::chrono::high_resolution_clock::now();

        // Perform range sum queries
        for (size_t q = 0; q < n_queries; ++q) {
            logger.log_section("Range Sum Query " + std::to_string(q + 1) + " begins");

            auto single_query_start = std::chrono::high_resolution_clock::now();

            RegXS left_index, right_index;

            size_t leafStart = (1ULL << (depth - 1));
            size_t max_leaf_idx = (1 << (depth - 1)) - 1;
            size_t left_val = (q % (1 << (depth - 1)));
            size_t right_val = left_val + (q % (max_leaf_idx + 1 - left_val));

            left_index.set(tio.player() == 0 ? (leafStart + left_val) : 0);
            right_index.set(tio.player() == 0 ? (leafStart + right_val) : 0);

#ifdef SEGTREE_VERBOSE
            auto recons_left = mpc_reconstruct(tio, yield, left_index);
            auto recons_right = mpc_reconstruct(tio, yield, right_index);
            oss.str("");
            oss << "Range Sum Query [" << recons_left << ", " << recons_right << "]\n";
            logger.log_output(oss.str());
#endif

            segTree.RangeSum(tio, mpcio, yield, left_index, right_index, &logger, q + 1);

            auto single_query_end = std::chrono::high_resolution_clock::now();
            auto single_query_duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(single_query_end - single_query_start).count();

            logger.log_metric("query", q + 1, "total_query_time", single_query_duration);
            logger.log_output("Range Sum Query " + std::to_string(q + 1) + " ends\n");
        }

        auto query_end = std::chrono::high_resolution_clock::now();
        auto query_duration = std::chrono::duration_cast<std::chrono::milliseconds>(query_end - query_start);

        logger.log_section("Range Sum Stats");
        oss.str("");
        oss << "Total Range Queries Time: " << query_duration.count() / 1000.0 << " seconds\n";
        logger.log_output(oss.str());

        logger.log_metric("query_summary", 0, "total_all_queries_time", query_duration.count());
    });
}
