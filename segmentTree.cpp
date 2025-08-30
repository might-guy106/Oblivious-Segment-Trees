#include <functional>
#include "mpcops.hpp"
#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include "rdpf.hpp"
#include "shapes.hpp"
#include "segmentTree.hpp"

#define SEGTREE_VERBOSE

size_t arrayIndexToLevelPos(size_t idx) {
    // Compute level as the floor of log2(idx)
    size_t level = static_cast<size_t>(std::log2(idx));
    // For a complete binary tree with root at index 1,
    // leftmost index in this level is 2^level.
    size_t pos = idx - (1ULL << level);

    return pos;
}

void SegmentTree::init(MPCTIO &tio, yield_t & yield) {
    auto SegTreeArray = oram.flat(tio, yield);
    std::cout << "depth in init is: " << depth << std::endl;
    size_t num_leaves = 1ULL << (depth - 1);
    size_t leaf_start = num_leaves;

    // Create and initialize segment tree array
    std::vector<uint64_t> segTree(num_items, 0);

    // Initialize leaf nodes with values 0, 100, 200, ..., (num_leaves-1)*100
    for (size_t i = 0; i < num_leaves; ++i) {
        segTree[leaf_start + i] = i * 100;
    }

    // Build internal nodes bottom-up by computing range sums
    for (int i = num_leaves - 1; i >= 1; --i) {
        segTree[i] = segTree[2*i] + segTree[2*i + 1];
    }

    SegTreeArray.init([this, &segTree] (size_t i) -> size_t {
        if (i >= 1 && i < num_items) {
            return segTree[i];
        } else {
            return size_t(0x7fffffffffffffff);
        }
    });

   auto isEvenArray = isEven.flat(tio, yield);
   isEvenArray.init([this] (size_t i) -> size_t {
     if (i >= 1 && i < num_items) {
          return (i % 2 == 0) ? 1 : 0;
     } else {
          return size_t(1); // 0 is even, just for consistency
     }
    });

    auto nextLArray = nextL.flat(tio, yield);
    nextLArray.init([this] (size_t i) -> size_t {
        if (i >= 1 && i < num_items) {
            if (i > 1 && ((i + 1) & i) == 0) // i is one less than a power of two.
                return ((arrayIndexToLevelPos(i) + 1)/2 - 1);
            else
                return ((arrayIndexToLevelPos(i) + 1)/2);
        } else {
            return size_t(0);
        }
    });

    auto nextRArray = nextR.flat(tio, yield);
    nextRArray.init([this] (size_t i) -> size_t {
        if (i >= 1 && i < num_items) {
            if ((i & (i - 1)) == 0) // i is a power of two.
                return size_t(0);
            else
                return ((arrayIndexToLevelPos(i) - 1)/2);
        } else {
            return size_t(0);
        }
    });

    auto parentArray = parent.flat(tio, yield);
    parentArray.init([this] (size_t i) -> size_t {
        if (i >= 1 && i < num_items) {
            return i / 2;
        } else {
            return size_t(0);
        }
    });
}

void SegmentTree::printSegmentTree(MPCTIO &tio, yield_t & yield) {
    auto SegTreeArray = oram.flat(tio, yield);
    auto SegTreeRecons = SegTreeArray.reconstruct();
    for(size_t i=1; i<num_items; i++) {
        std::cout << "SegTreeArray[" << i << "] = " << SegTreeRecons[i].share() << std::endl;
    }
}

void SegmentTree::getBitVector(MPCTIO &tio, yield_t & yield, Duoram < RegXS > &bitVec, RegAS left, RegAS right) {

    auto bitVecArray = bitVec.flat(tio, yield);
    auto isEvenArray = isEven.flat(tio, yield);
    auto nextLArray = nextL.flat(tio, yield);
    auto nextRArray = nextR.flat(tio, yield);

    uint32_t height = static_cast<uint32_t>(std::log2(num_items));

    RegXS incl;
    incl.set(tio.player()==0 ? 1 : 0);
    RegXS excl;
    excl.set(0);
    RegBS isDone;
    isDone.set(0);
    RegBS zero;
    zero.set(0);
    RegBS one;
    one.set(tio.player()==0 ? 1 : 0);

    for(uint32_t i=1; i<=height; i++)
    {
        uint32_t level = height - i;
        typename Duoram < RegXS > ::Flat bitVecLevel(bitVecArray, tio, yield, (1ULL << level), (1ULL << level)+1);
        typename Duoram < RegXS > ::Flat isEvenLevel(isEvenArray, tio, yield, (1ULL << level), (1ULL << level)+1);
        typename Duoram < RegAS > ::Flat nextLLevel(nextLArray, tio, yield, (1ULL << level), (1ULL << level)+1);
        typename Duoram < RegAS > ::Flat nextRLevel(nextRArray, tio, yield, (1ULL << level), (1ULL << level)+1);

        CDPF cdpf = tio.cdpf(yield);
        RegAS diff = right - left;
        auto[lt_c, eq_c, gt_c] = cdpf.compare(tio, yield, diff, tio.aes_ops());

        RegXS leftIncluded = bitVecLevel[left];
        RegXS rightIncluded = bitVecLevel[right];

        RegBS valid;
        mpc_or(tio, yield, valid, eq_c, gt_c);
        mpc_or(tio, yield, isDone, valid ^ one, isDone);



        RegXS isEvenL = isEvenLevel[left];
        RegBS inclusionL = isEvenL.bitat(0);
        RegBS f;
        mpc_and(tio, yield, f, valid, (inclusionL ^ one));
        mpc_select(tio, yield, leftIncluded, f, excl, incl);

        RegXS isEvenR = isEvenLevel[right];
        RegBS inclusionR = isEvenR.bitat(0);
        RegBS g;
        mpc_and(tio, yield, g, valid, inclusionR);
        mpc_select(tio, yield, rightIncluded, g, excl, incl);

        // checks if already done . if done then both left and right are 0
        mpc_select(tio, yield, leftIncluded, isDone, leftIncluded, excl);
        mpc_select(tio, yield, rightIncluded, isDone, rightIncluded, excl);

        // if not done and if eq then leftinclude = righinclude = 1 and isDone = 1
        RegBS temp;
        mpc_and(tio, yield, temp, isDone ^ one, eq_c);
        mpc_select(tio, yield, leftIncluded, temp, leftIncluded, incl);
        mpc_select(tio, yield, rightIncluded, temp, rightIncluded, incl);
        mpc_or(tio, yield, isDone, eq_c, isDone);

        bitVecLevel[left] = leftIncluded;
        bitVecLevel[right] = rightIncluded;

        #ifdef SEGTREE_VERBOSE
        value_t recons = mpc_reconstruct(tio, yield, valid);
        auto leftIndRecons = mpc_reconstruct(tio, yield, left);
        auto rightIndRecons = mpc_reconstruct(tio, yield, right);
        auto lIsIncluded = mpc_reconstruct(tio, yield, leftIncluded);
        auto rIsIncluded = mpc_reconstruct(tio, yield, rightIncluded);
        std::cout << "Level: " << level << " Left Node [" << leftIndRecons << "] (bitVec[ "<< ((1ULL << level) + leftIndRecons) << "]) isincluded: "  << lIsIncluded << " Right Node [" << rightIndRecons << "] (bitVec[" << ((1ULL << level) + rightIndRecons) << "]) isincluded: "  << rIsIncluded << " valid: " << recons << std::endl;
        #endif

        left = nextLLevel[left];
        right = nextRLevel[right];
    }
}

void SegmentTree::RangeSum(MPCTIO &tio, yield_t & yield, RegAS left, RegAS right) {
    Duoram < RegXS > bitVec(tio.player(), num_items);
    getBitVector(tio, yield, bitVec, left, right);


    auto bitVecArray = bitVec.flat(tio, yield);
    auto SegTreeArray = oram.flat(tio, yield);

    RegAS sum;
    sum.set(0);

    for(size_t i=1; i<num_items; i++) {
        RegXS element = bitVecArray[i];
        RegBS incl = element.bitat(0);
        RegAS val = SegTreeArray[i];
        RegAS zero;
        zero.set(0);

        RegAS sum1;
        mpc_select(tio, yield, sum1, incl, zero, val);

        sum.ashare += sum1.ashare;
    }

    #ifdef SEGTREE_VERBOSE
    value_t answer = mpc_reconstruct(tio, yield, sum);
    std::cout << "Sum = " << answer << std::endl;
    #endif
}

void SegmentTree::Update(MPCTIO &tio, yield_t & yield, RegAS index, RegAS value) {
    auto SegTreeArray = oram.flat(tio, yield);
    auto parentArray = parent.flat(tio, yield);

    RegAS disp;
    disp.set(tio.player() == 0 ? (1 << (depth - 1)) : 0);

    RegAS index1 = index + disp;


    RegAS currVal = SegTreeArray[index1];
    RegAS diff = value - currVal;

    #ifdef SEGTREE_VERBOSE
    auto recons_index = mpc_reconstruct(tio, yield, index);
    auto recons_index1 = mpc_reconstruct(tio, yield, index1);
    auto recons_currVal = mpc_reconstruct(tio, yield, currVal);
    auto recons_newvalue = mpc_reconstruct(tio, yield, value);
    auto recons_diff = mpc_reconstruct(tio, yield, diff);
    std::cout << "Index to be updated in the original array = " << recons_index << std::endl;
    std::cout << "Index to be updated in the segment tree array = " << recons_index1 << std::endl;
    std::cout << "Current Value at index = " << recons_currVal << std::endl;
    std::cout << "New Value to be updated = " << recons_newvalue << std::endl;
    std::cout << "Diff = " << recons_diff << std::endl;
    #endif

    SegTreeArray[index1] = value;
    for(size_t i=1; i<=(depth-1); i++) {
        RegAS parentIndex = parentArray[index1];
        auto recons_parentIndex = mpc_reconstruct(tio, yield, parentIndex);
        SegTreeArray[parentIndex] += diff;
        #ifdef SEGTREE_VERBOSE
        auto recons_updatedParent = mpc_reconstruct(tio, yield, SegTreeArray[parentIndex]);
        std::cout << "Updated Parent Index = " << recons_parentIndex << " with value = " << recons_updatedParent << std::endl;
        #endif
        index1 = parentIndex;
    }
}


void SegTree(MPCIO &mpcio, const PRACOptions &opts, char **args) {
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

    address_t len = (1<<depth);

    MPCTIO tio(mpcio, 0, opts.num_cpu_threads);

    run_coroutines(tio, [&tio, len, depth, n_updates, n_queries] (yield_t &yield) {
        SegmentTree segTree(tio.player(), len, depth);
        segTree.init(tio, yield);

        std::cout << "===== Segment Tree Initialized =====" << std::endl;
        std::cout << "Depth: " << depth << ", Size: " << len << std::endl;
        std::cout << "Updates: " << n_updates << ", Queries: " << n_queries << std::endl;

        #ifdef SEGTREE_VERBOSE
        // Print initial segment tree
        std::cout << "\n===== Initial Segment Tree =====" << std::endl;
        segTree.printSegmentTree(tio, yield);
        #endif

        // Perform updates
        for (size_t u = 0; u < n_updates; ++u) {
            std::cout << "\n===== Update " << (u + 1) << " begins =====" << std::endl;
            RegAS index;
            size_t idx_val = u % (1 << (depth - 1));
            index.set(tio.player() == 0 ? idx_val : 0); // Cycle through leaf indices


            RegAS value;
            size_t val_to_set = (u + 1) * 50;
            value.set(tio.player() == 0 ? val_to_set : 0); // Use different values for each update

            segTree.Update(tio, yield, index, value);
            std::cout << "Update " << (u + 1) << " ends" << std::endl;
        }

        #ifdef SEGTREE_VERBOSE
        // Print updated segment tree
        std::cout << "\n===== Updated Segment Tree =====" << std::endl;
        segTree.printSegmentTree(tio, yield);
        #endif

        // Perform range sum queries
        for (size_t q = 0; q < n_queries; ++q) {
            std::cout << "\n===== Range Sum Query " << (q + 1) << " begins =====" << std::endl;
            RegAS left_index, right_index;

            // Ensure left < right for valid range queries
            size_t max_leaf_idx = (1 << (depth - 1)) - 1;
            size_t left_val = q % (1 << (depth - 1));  // Start anywhere in valid range
            size_t right_val = left_val + (q % (max_leaf_idx - left_val)); // Ensuring right > left

            left_index.set(tio.player() == 0 ? left_val : 0); // Start of range
            right_index.set(tio.player() == 0 ? right_val : 0); // End of range

            #ifdef SEGTREE_VERBOSE
            auto recons_left = mpc_reconstruct(tio, yield, left_index);
            auto recons_right = mpc_reconstruct(tio, yield, right_index);
            std::cout << "Range Sum Query [" << recons_left << ", " << recons_right << "]" << std::endl;
            #endif

            segTree.RangeSum(tio, yield, left_index, right_index);
            std::cout << "Range Sum Query " << (q + 1) << " ends" << std::endl;
        }
    });
}
