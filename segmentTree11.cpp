#include "segmentTree11.hpp"
#include "coroutine.hpp"
#include "duoram.hpp"
#include "mpcops.hpp"
#include "rdpf.hpp"

/*
The segment tree data structure in SegmentTree11 uses an ARRAY layout (heap-style indexing).

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
using the level-wise segment tree structure. The algorithm processes levels using
level-specific Flat objects for performance in MPC.

Indexing notes:
- Parent of node at array index i is at i/2 (i >> 1). In this implementation,
  indices are maintained as level-relative indices while traversing levels, so
  shifting right still yields the parent index at the next higher level.
- Unlike SegmentTree9, SegmentTree11 does NOT maintain a separate Sibling ORAM.
  Instead, the sibling index is derived via MPC from the current index (using the
  index parity / last bit), and then used to read the sibling value from the Tree
  ORAM.

Algorithm Overview:
1. Start with leaf-level indices (left, right)
2. Precompute the path indices for all levels going up the tree
3. For each level from leaves to root:
   - Create a level-specific Flat object for the Tree ORAM
   - Compute Additive Shares of XOR shares of left and right indices
   - Use additive-share comparison to compute isNotDone
   - Use the least significant bit of XOR shares of left/right indices as a
     boolean-share isOdd
   - Compute left/right sibling indices via MPC (no sibling ORAM), then read the
     sibling values from the Tree ORAM
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

void SegmentTree11::init(MPCTIO &tio, yield_t &yield) {

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
}

// helper function to reconstruct and print SegTreeArray
void SegmentTree11::printSegmentTree(MPCTIO &tio, yield_t &yield) {
    auto SegTreeArray = TreeOram.flat(tio, yield);
    auto SegTreeRecons = SegTreeArray.reconstruct();
    for (size_t i = 1; i < num_items; i++) {
        std::cout << "SegTreeArray[" << i << "] = " << SegTreeRecons[i].share() << std::endl;
    }
}

// Range sum over absolute indices in the tree.
// NOTE: left and right are absolute tree indices (NOT level-relative leaf indices).
RegAS SegmentTree11::RangeSum(MPCTIO &tio, MPCIO &mpcio, yield_t &yield, RegXS leftLevelIndex, RegXS rightLevelIndex) {

    auto SegTreeArray = TreeOram.flat(tio, yield);

    // Pre-compute the path indices for all levels
    std::vector<RegXS> leftPathIndex(depth);
    std::vector<RegXS> rightPathIndex(depth);

    // Compute parent indices for all levels going up the tree
    for (uint32_t i = 1; i <= depth; i++) {
        size_t level = depth - i; // Current level (leaf to root)

        leftPathIndex[level] = leftLevelIndex;
        rightPathIndex[level] = rightLevelIndex;

        leftLevelIndex >>= 1;
        rightLevelIndex >>= 1;
    }

    // Compute sum directly while identifying nodes to include
    std::vector<RegAS> levelSums(depth);
    for (size_t level = 0; level < depth; level++) {
        levelSums[level].set(0);
    }

    std::vector<coro_t> sum_coroutines;
    for (uint32_t i = 1; i <= depth; i++) {
        size_t level = depth - i;

        sum_coroutines.emplace_back([&tio, level, i, this, &leftPathIndex, &rightPathIndex, &SegTreeArray,
                                     &levelSums](yield_t &sub_yield) {
            size_t levelStart = 1ULL << level;
            size_t levelLength = (1ULL << level) + 1;

            // Create one, zero which are required in mpc operations further
            RegAS one;
            one.set(tio.player());

            RegAS zero;
            zero.set(0);

            typename Duoram<RegAS>::Flat levelTreeArray(SegTreeArray, tio, sub_yield, levelStart, levelLength);
            // typename Duoram<RegAS>::Flat levelSiblingArray(SiblingArray, tio, sub_yield, levelStart, levelLength);

            // level index which are already stored
            RegXS leftIndex = leftPathIndex[level];
            RegXS rightIndex = rightPathIndex[level];

            RegAS leftIndexAS;
            RegAS rightIndexAS;
            run_coroutines(
                sub_yield,
                [&tio, leftIndex, &leftIndexAS](yield_t &yield) {
                    mpc_xs_to_as(tio, yield, leftIndexAS, leftIndex);
                },
                [&tio, rightIndex, &rightIndexAS](yield_t &yield) {
                    mpc_xs_to_as(tio, yield, rightIndexAS, rightIndex);
                });

            // if left and right are adjacent or same then the marking of nodes is
            // already completed so we set isDone here
            RegBS isNotDone;
            CDPF cdpf1 = tio.cdpf(sub_yield);
            RegAS diff1 = rightIndexAS - (leftIndexAS + one); // diff1 = right - (left + 1)
            auto [lt_c1, eq_c1, gt_c1] = cdpf1.compare(tio, sub_yield, diff1, tio.aes_ops());

            isNotDone = gt_c1; // isNotDone = (right > left + 1)

            // include left and right independently based on isValid and isDifferent
            // and if left is left child and right is right child if xor of last bit
            // is 1 then its odd -> right child else left child if leftLastBit is 1
            // then left is right child else left child
            RegBS leftLastBit = leftIndex.bitat(0);
            RegBS rightLastBit = rightIndex.bitat(0);

            RegAS leftSum, rightSum;
            run_coroutines(
                sub_yield,
                [&tio, &levelTreeArray, one, isNotDone, leftIndex, leftIndexAS, &leftSum, leftLastBit](yield_t &yield) {

                    RegAS plusOne = leftIndexAS + one;
                    RegAS minusOne = leftIndexAS - one;
                    // leftLastBit=0 (left child)  → sibling = i+1 (plusOne)
                    // leftLastBit=1 (right child) → sibling = i-1 (minusOne)
                    RegAS leftSibIndex;
                    mpc_select(tio, yield, leftSibIndex, leftLastBit, plusOne, minusOne);

                    auto levelTreeArrayCoro = levelTreeArray.context(yield);
                    RegAS leftSibValue = levelTreeArrayCoro[leftSibIndex];

                    RegBS isLeftIncluded;
                    mpc_and(tio, yield, isLeftIncluded, isNotDone, leftLastBit ^ tio.player());
                    mpc_flagmult(tio, yield, leftSum, isLeftIncluded, leftSibValue);
                },
                [&tio, &levelTreeArray, one, isNotDone, rightIndex, rightIndexAS, &rightSum, rightLastBit](yield_t &yield) {

                    RegAS plusOne = rightIndexAS + one;
                    RegAS minusOne = rightIndexAS - one;
                    // rightLastBit=0 (left child)  → sibling = i+1 (plusOne)
                    // rightLastBit=1 (right child) → sibling = i-1 (minusOne)
                    RegAS rightSibIndex;
                    mpc_select(tio, yield, rightSibIndex, rightLastBit, plusOne, minusOne);

                    auto levelTreeArrayCoro = levelTreeArray.context(yield);
                    RegAS rightSibValue = levelTreeArrayCoro[rightSibIndex];

                    RegBS isRightIncluded;
                    mpc_and(tio, yield, isRightIncluded, isNotDone, rightLastBit);
                    mpc_flagmult(tio, yield, rightSum, isRightIncluded, rightSibValue);
                });

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
                    [&tio, &lt_c2, &eq_c2, &gt_c2, &zero, &one, &leftIndex, &levelTreeArray,
                     &leftLeafContribution](yield_t &yield) {
                        // if right >= left add left
                        // that is equivalent to !(left > right) => lt_c2 ^ tio.player()
                        auto levelTreeArrayCoro = levelTreeArray.context(yield);
                        RegAS leftLeafValue = levelTreeArrayCoro[leftIndex];
                        mpc_select(tio, yield, leftLeafContribution, lt_c2 ^ tio.player(), zero, leftLeafValue);
                    },
                    [&tio, &lt_c2, &eq_c2, &gt_c2, &zero, &one, &rightIndex, &levelTreeArray,
                     &rightLeafContribution](yield_t &yield) {
                        // Only add right leaf if it's different from left (to avoid
                        // double counting) that is if right > left
                        auto levelTreeArrayCoro = levelTreeArray.context(yield);
                        RegAS rightLeafValue = levelTreeArrayCoro[rightIndex];
                        mpc_select(tio, yield, rightLeafContribution, gt_c2, zero, rightLeafValue);
                    });

                levelSums[level] += leftLeafContribution + rightLeafContribution;
            }
        });
    }

    run_coroutines(tio, sum_coroutines);

    // Combine all level sums
    RegAS totalSum;
    totalSum.set(0);
    for (size_t level = 0; level < depth; level++) {
        totalSum += levelSums[level];
    }

    return totalSum;
}

// Main Update function
void SegmentTree11::Update(MPCTIO &tio, MPCIO &mpcio, yield_t &yield, RegXS index, RegAS value) {

    auto SegTreeArray = TreeOram.flat(tio, yield);

    size_t leafLevel = depth - 1;
    size_t leafStart = 1ULL << leafLevel;

    // Access last level Flat object and get current value to calculate difference
    size_t leafLength = (1ULL << leafLevel) + 1;
    typename Duoram<RegAS>::Flat leafLevelArray(SegTreeArray, tio, yield, leafStart, leafLength);
    RegAS currVal = leafLevelArray[index];
    RegAS diff = value - currVal;

    // Debugging and intermediate reconstructions

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

    // Phase 2: Parallelize the updates across all levels
    std::vector<coro_t> update_coroutines;
    for (size_t i = 1; i <= depth; i++) {
        size_t level = depth - i;

        update_coroutines.emplace_back([&tio, level, this, &diff, &updatePathIndex, &SegTreeArray](yield_t &sub_yield) {
            size_t levelStart = 1ULL << level;
            size_t levelLength = (1ULL << level) + 1;

            typename Duoram<RegAS>::Flat levelTreeArray(SegTreeArray, tio, sub_yield, levelStart, levelLength);

            // Update the value field at this level's index
            levelTreeArray[updatePathIndex[level]] += diff;
        });
    }

    run_coroutines(tio, update_coroutines);
}

// Main function to run Segment Tree 8 operations
void SegTree11(MPCIO &mpcio, const PRACOptions &opts, char **args) {
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
    exp_id_stream << "st11_d" << (int)depth << "_u" << n_updates << "_q" << n_queries << "_"
                  << std::put_time(&tm, "%Y%m%d_%H%M%S");
    std::string experiment_id = exp_id_stream.str();

    run_coroutines(tio, [&tio, &mpcio, len, depth, n_updates, n_queries](yield_t &yield) {
        SegmentTree11 segTree(tio.player(), len, depth);
        segTree.init(tio, yield);

        // Updates
        for (size_t u = 0; u < n_updates; ++u) {
            RegXS index;
            size_t idx_val = (u % (1 << (depth - 1)));
            index.set(tio.player() == 0 ? idx_val : 0);

            RegAS value;
            size_t val_to_set = (u + 1) * 50;
            value.set(tio.player() == 0 ? val_to_set : 0);

            segTree.Update(tio, mpcio, yield, index, value);
        }

        // Range Queries
        for (size_t q = 0; q < n_queries; ++q) {
            RegXS left_index, right_index;

            size_t max_leaf_idx = (1 << (depth - 1)) - 1;
            size_t left_val = (q % (1 << (depth - 1)));
            size_t right_val = left_val + (q % (max_leaf_idx + 1 - left_val));

            left_index.set(tio.player() == 0 ? left_val : 0);
            right_index.set(tio.player() == 0 ? right_val : 0);

            RegAS sum = segTree.RangeSum(tio, mpcio, yield, left_index, right_index);
        }
    });
}
