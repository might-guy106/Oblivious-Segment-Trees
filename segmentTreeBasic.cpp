#include "segmentTreeBasic.hpp"
#include "coroutine.hpp"
#include "duoram.hpp"
#include "mpcops.hpp"
#include "rdpf.hpp"
#include "types.hpp"

/*
The segment tree data structure in SegmentTreeBasic uses an ARRAY layout (heap-style indexing).

ARRAY STRUCTURE:
- Complete binary tree stored in an array with root at index 1
- Each level k has 2^k nodes
- Total nodes in the tree: 2^d (indices [1 .. 2^d - 1] are used; index 0 is unused)
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
RANGE SUM QUERY OPERATION (Absolute-index Processing)

The range sum query computes the sum of elements in a given range [left, right]
where left and right are ABSOLUTE indices in the flattened/heap segment tree
(not level-relative leaf indices).

Indexing notes:
- Parent of node at array index i is at i/2 (i >> 1).
- Sibling pointers are stored in a separate Sibling ORAM (SiblingOram) as ABSOLUTE
  tree indices: sibling(i) = (i ^ 1) for i >= 1.
- This implementation does NOT create per-level Flat sub-ORAMs. All accesses use
  the main ORAM flats directly.

Algorithm Overview:
1. Start with absolute indices (left, right)
2. Precompute the path indices for all levels going up the tree by repeated >> 1
3. For each level from leaves to root:
   - Use additive-share comparison to compute isNotDone
   - Use the least significant bit of XOR shares of left/right indices as a
     boolean-share isOdd
   - Read left/right sibling indices from the Sibling ORAM (absolute indices),
     then read sibling values from the Tree ORAM
   - Add left sibling value if left is a left child: (!isOdd && isNotDone)
   - Add right sibling value if right is a right child: (isOdd && isNotDone)
*/

/*
UPDATE OPERATION (Absolute-index Propagation)

The update operation modifies a single node at an ABSOLUTE tree index and
propagates the change up through all levels by repeatedly shifting to the parent.

Algorithm Overview:
1. Calculate difference between new and current values at the absolute index
2. For each level from leaf to root:
   - Add diff to the current absolute index in the Tree ORAM
   - Move to parent using index >>= 1
*/

void SegmentTreeBasic::init(MPCTIO &tio, yield_t &yield) {

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
    // Store absolute sibling tree index (heap layout.
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
void SegmentTreeBasic::printSegmentTree(MPCTIO &tio, yield_t &yield) {
    auto SegTreeArray = TreeOram.flat(tio, yield);
    auto SegTreeRecons = SegTreeArray.reconstruct();
    for (size_t i = 1; i < num_items; i++) {
        std::cout << "SegTreeArray[" << i << "] = " << SegTreeRecons[i].share() << std::endl;
    }
}

// Main function to compute range sum
RegAS SegmentTreeBasic::RangeSum(MPCTIO &tio, MPCIO &mpcio, yield_t &yield, RegXS leftIndex, RegXS rightIndex) {

    auto SegTreeArray = TreeOram.flat(tio, yield);
    auto SiblingArray = SiblingOram.flat(tio, yield);

    // Pre-compute the path indices for all levels
    std::vector<RegXS> leftPathIndex(depth);
    std::vector<RegXS> rightPathIndex(depth);

    // Compute parent indices for all levels going up the tree.
    for (uint32_t i = 1; i <= depth; i++) {
        size_t level = depth - i; // Current level (leaf to root)

        leftPathIndex[level] = leftIndex;
        rightPathIndex[level] = rightIndex;

        leftIndex >>= 1;
        rightIndex >>= 1;
    }

    // For storing sum of nodes in a level that are included in the range sum computation
    std::vector<RegAS> levelSums(depth);
    for (size_t level = 0; level < depth; level++) {
        levelSums[level].set(0);
    }

    std::vector<coro_t> sum_coroutines;
    for (uint32_t i = 1; i <= depth; i++) {
        size_t level = depth - i;

        sum_coroutines.emplace_back([&tio, level, i, this, &leftPathIndex, &rightPathIndex, &SegTreeArray,
                                     &SiblingArray, &levelSums](yield_t &sub_yield) {
            RegAS one;
            one.set(tio.player());

            RegAS zero;
            zero.set(0);

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

            RegBS isNotDone;
            CDPF cdpf1 = tio.cdpf(sub_yield);
            RegAS diff1 = rightIndexAS - (leftIndexAS + one); // diff1 = right - (left + 1)
            auto [lt_c1, eq_c1, gt_c1] = cdpf1.compare(tio, sub_yield, diff1, tio.aes_ops());

            isNotDone = gt_c1; // isNotDone = (right > left + 1)

            // include left and right independently based on isValid and isDifferent
            // and if left is left child and right is right child
            // if xor of last bit is 1 then its odd -> right child else left child
            // if leftLastBit is 1 then left is right child else left child
            RegBS leftLastBit = leftIndex.bitat(0);
            RegBS rightLastBit = rightIndex.bitat(0);

            RegAS leftSum, rightSum;
            run_coroutines(
                sub_yield,
                [&tio, &SegTreeArray, &SiblingArray, isNotDone, leftIndex, &leftSum, leftLastBit](yield_t &yield) {
                    auto SiblingArrayCoro = SiblingArray.context(yield);
                    auto SegTreeArrayCoro = SegTreeArray.context(yield);

                    // Indices are already absolute in the flattened tree
                    RegAS leftSibIndex = SiblingArrayCoro[leftIndex];
                    RegAS leftSibValue = SegTreeArrayCoro[leftSibIndex];

                    RegBS isLeftIncluded;
                    mpc_and(tio, yield, isLeftIncluded, isNotDone, leftLastBit ^ tio.player());
                    mpc_flagmult(tio, yield, leftSum, isLeftIncluded, leftSibValue);
                },
                [&tio, &SegTreeArray, &SiblingArray, isNotDone, rightIndex, &rightSum, rightLastBit](yield_t &yield) {
                    auto SiblingArrayCoro = SiblingArray.context(yield);
                    auto SegTreeArrayCoro = SegTreeArray.context(yield);

                    // Indices are already absolute in the flattened tree
                    RegAS rightSibIndex = SiblingArrayCoro[rightIndex];
                    RegAS rightSibValue = SegTreeArrayCoro[rightSibIndex];

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
void SegmentTreeBasic::Update(MPCTIO &tio, MPCIO &mpcio, yield_t &yield, RegXS index, RegAS value) {

    auto SegTreeArray = TreeOram.flat(tio, yield);

    RegAS currVal = SegTreeArray[index];
    RegAS diff = value - currVal;

    std::vector<RegXS> updatePathIndex(depth);

    for (size_t i = 1; i <= depth; i++) {
        size_t level = depth - i; // Current level (leaf to root)
        updatePathIndex[level] = index;
        index >>= 1;
    }

    std::vector<coro_t> update_coroutines;
    for (size_t i = 1; i <= depth; i++) {
        size_t level = depth - i;

        update_coroutines.emplace_back([&tio, level, this, &diff, &updatePathIndex, &SegTreeArray](yield_t &sub_yield) {
            auto SegTreeArrayCoro = SegTreeArray.context(sub_yield);

            RegXS absIndex = updatePathIndex[level];
            SegTreeArrayCoro[absIndex] += diff;
        });
    }

    run_coroutines(tio, update_coroutines);

    (void)mpcio;
}

// Main function to run Segment Tree 8 operations
void SegTreeBasic(MPCIO &mpcio, const PRACOptions &opts, char **args) {
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
    exp_id_stream << "stBasic_d" << (int)depth << "_u" << n_updates << "_q" << n_queries << "_"
                  << std::put_time(&tm, "%Y%m%d_%H%M%S");
    std::string experiment_id = exp_id_stream.str();

    run_coroutines(tio, [&tio, &mpcio, len, depth, n_updates, n_queries](yield_t &yield) {
        SegmentTreeBasic segTree(tio.player(), len, depth);
        segTree.init(tio, yield);

        // Updates
        for (size_t u = 0; u < n_updates; ++u) {
            RegXS index;
            size_t leafStart = (1ULL << (depth - 1));
            size_t idx_val = leafStart + (u % (1 << (depth - 1)));
            index.set(tio.player() == 0 ? idx_val : 0);

            RegAS value;
            size_t val_to_set = (u + 1) * 50;
            value.set(tio.player() == 0 ? val_to_set : 0);

            segTree.Update(tio, mpcio, yield, index, value);
        }

        // Range Queries
        for (size_t q = 0; q < n_queries; ++q) {
            RegXS left_index, right_index;

            size_t leafStart = (1ULL << (depth - 1));
            size_t max_leaf_idx = (1 << (depth - 1)) - 1;
            size_t left_val = (q % (1 << (depth - 1)));
            size_t right_val = left_val + (q % (max_leaf_idx + 1 - left_val));

            left_index.set(tio.player() == 0 ? (leafStart + left_val) : 0);
            right_index.set(tio.player() == 0 ? (leafStart + right_val) : 0);

            RegAS sum = segTree.RangeSum(tio, mpcio, yield, left_index, right_index);
        }
    });
}
