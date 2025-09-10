#include <functional>
#include "mpcops.hpp"
#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include "rdpf.hpp"
#include "shapes.hpp"
#include "segmentTree6.hpp"

#define SEGTREE_VERBOSE
// To enable timing/stat instrumentation define SEGTREE_VERBOSE2
// #define SEGTREE_VERBOSE2

#ifdef SEGTREE_VERBOSE2
#define STATS_PRE() do { tio.sync_lamport(); mpcio.reset_stats(); tio.reset_lamport(); } while(0)
#define STATS_POST(MSG) do { tio.sync_lamport(); std::cout << MSG << std::endl; mpcio.dump_stats(std::cout); mpcio.reset_stats(); tio.reset_lamport(); } while(0)
#else
#define STATS_PRE() do {} while(0)
#define STATS_POST(MSG) do {} while(0)
#endif

// Helper functions for the new level-wise array structure    
size_t findLevel6(size_t index, size_t depth) {
    for (size_t level = 0; level < depth; level++) {
        size_t levelStart = (1ULL << level) + level - 1;
        size_t levelEnd = (1ULL << (level + 1)) + level - 1;
        if (index >= levelStart && index <= levelEnd) {
            return level;
        }
    }
    return SIZE_MAX; // Invalid
}

size_t getLevelStart6(size_t level) {
    return (1ULL << level) + level - 1;  // 2^level + level - 1
}

size_t getLevelLength6(size_t level) {
    return (1ULL << level) + 1;  // 2^level + 1
}

size_t parentIndex6(size_t index, size_t depth) {
    size_t currentLevel = findLevel6(index, depth);
    
    if (currentLevel <= 0) {
        return index;  // Root level
    }
    
    size_t levelStart = getLevelStart6(currentLevel);
    
    // If it's the extra node, return itself
    if (index == levelStart) {
        return 0;
    }
    
    // For regular nodes
    size_t positionInLevel = index - levelStart;
    size_t parentPosition = (positionInLevel - 1) / 2 + 1;
    
    return parentPosition;
}

size_t getLeftSiblingIndex6(size_t index, size_t depth) {
    size_t currentLevel = findLevel6(index, depth);
    
    if (currentLevel <= 0) {
        return index;  // Root level, return itself
    }
    
    size_t levelStart = getLevelStart6(currentLevel);
    
    // If it's the extra node, return itself
    if (index == levelStart) {
        return 0;
    }
    
    // For regular nodes
    size_t positionInLevel = index - levelStart;
    
    if(positionInLevel % 2 == 1) { // odd index -> left child
        return positionInLevel + 1; // return sibling
    } else {
        return 0;
    }
}

size_t getRightSiblingIndex6(size_t index, size_t depth) {
    size_t currentLevel = findLevel6(index, depth);
    
    if (currentLevel <= 0) {
        return index;  // Root level, return itself
    }
    
    size_t levelStart = getLevelStart6(currentLevel);
    
    // If it's the extra node, return itself
    if (index == levelStart) {
        return 0;
    }
    
    // For regular nodes
    size_t positionInLevel = index - levelStart;

    if(positionInLevel % 2 == 0) { // even index -> right child
        return positionInLevel - 1; // return sibling
    } else {
        return 0;
    }
}

void SegmentTree6::init(MPCTIO &tio, yield_t & yield) {
    auto SegTreeArray = oram.flat(tio, yield);
    
    // Create and initialize new level-wise segment tree array
    std::vector<uint64_t> segTree(num_items, 0);
    
    // Initialize the original binary tree values first (for mapping to new structure)
    size_t original_num_leaves = 1ULL << (depth - 1);
    size_t original_size = original_num_leaves * 2;
    std::vector<uint64_t> originalTree(original_size, 0);
    
    // Initialize leaf nodes with values 0, 100, 200, ..., (num_leaves-1)*100
    for (size_t i = 0; i < original_num_leaves; ++i) {
        originalTree[original_num_leaves + i] = i * 100;
    }
    
    // Build internal nodes bottom-up by computing range sums
    for (int i = original_num_leaves - 1; i >= 1; --i) {
        originalTree[i] = originalTree[2*i] + originalTree[2*i + 1];
    }
    
    // Map original tree to new level-wise structure
    for (size_t level = 0; level < depth; level++) {
        size_t levelStart = getLevelStart6(level);
        size_t levelLength = getLevelLength6(level);
        
        // Set extra node (position 0) to 0
        segTree[levelStart] = 0;
        
        // Map regular nodes
        for (size_t pos = 1; pos < levelLength; pos++) {
            size_t newIndex = levelStart + pos;
            
            if (level == depth - 1) {
                // Leaf level: map directly from original leaf values
                segTree[newIndex] = originalTree[original_num_leaves + pos - 1];
            } else {
                // Internal level: map from original internal nodes
                size_t originalIndex = (1ULL << level) + pos - 1;
                if (originalIndex < original_size) {
                    segTree[newIndex] = originalTree[originalIndex];
                }
            }
        }
    }

    SegTreeArray.init([this, &segTree] (size_t i) -> size_t {
        if (i < num_items) {
            return segTree[i];
        } else {
            return size_t(0x7fffffffffffffff);
        }
    });

    auto leftChildSiblingArray = leftChildSibling.flat(tio, yield);
    leftChildSiblingArray.init([this] (size_t i) -> size_t {
        if (i < num_items) {
            return getLeftSiblingIndex6(i, depth); // if odd (left) keep sibling else 0
        } else {
            return size_t(0);
        }
    });

    auto rightChildSiblingArray = rightChildSibling.flat(tio, yield);
    rightChildSiblingArray.init([this] (size_t i) -> size_t {
        if (i < num_items) {
            return getRightSiblingIndex6(i, depth); // if even (right) keep sibling else 0
        } else {
            return size_t(0);
        }
    });

    auto parentArray = parent.flat(tio, yield);
    parentArray.init([this] (size_t i) -> size_t {
        if (i < num_items) {
            return parentIndex6(i, depth);
        } else {
            return size_t(0);
        }
    });
}

void SegmentTree6::printSegmentTree(MPCTIO &tio, yield_t & yield) {
    auto SegTreeArray = oram.flat(tio, yield);
    auto SegTreeRecons = SegTreeArray.reconstruct();
    for(size_t i=0; i<num_items; i++) {
        std::cout << "SegTreeArray[" << i << "] = " << SegTreeRecons[i].share() << std::endl;
    }
}

void SegmentTree6::getBitVector(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, Duoram < RegXS > &bitVec, RegAS leftLevelIndex, RegAS rightLevelIndex) {

    auto bitVecArray = bitVec.flat(tio, yield);
    auto leftChildSiblingArray = leftChildSibling.flat(tio, yield);
    auto rightChildSiblingArray = rightChildSibling.flat(tio, yield);
    auto parentArray = parent.flat(tio, yield);


    RegXS incl;
    incl.set(tio.player()==0 ? 1 : 0);
    RegXS excl;
    excl.set(0);
    RegBS isDone;
    isDone.set(0);
    RegBS one;
    one.set(tio.player()==0 ? 1 : 0);
    RegBS zero;
    zero.set(0);
    
    for(uint32_t i=1; i<=depth; i++) {
        // level variable retained only for logging clarity
        size_t level = depth - i;
        size_t levelStart = getLevelStart6(level);
        size_t levelLength = getLevelLength6(level);
        
        typename Duoram < RegXS > ::Flat bitVecLevel(bitVecArray, tio, yield, levelStart, levelLength);
        typename Duoram < RegAS > ::Flat leftChildSiblingLevel(leftChildSiblingArray, tio, yield, levelStart, levelLength);
        typename Duoram < RegAS > ::Flat rightChildSiblingLevel(rightChildSiblingArray, tio, yield, levelStart, levelLength);
        typename Duoram < RegAS > ::Flat parentLevel(parentArray, tio, yield, levelStart, levelLength);

        // --- ORAM base reads: parents + conditional siblings (4 accesses) ---
        STATS_PRE();
        RegAS leftParent = parentLevel[leftLevelIndex];
        RegAS rightParent = parentLevel[rightLevelIndex];
        RegAS leftSibling = leftChildSiblingLevel[leftLevelIndex];   // 0 if left not a left child
        RegAS rightSibling = rightChildSiblingLevel[rightLevelIndex]; // 0 if right not a right child
        STATS_POST("[SEGTREE][BITVEC] ORAM Reads (parents+condSiblings) Stats (level=" + std::to_string(level) + ")");

        // --- CDPF compare for siblings ---
        STATS_PRE();
        
        CDPF cdpf2 = tio.cdpf(yield);
        RegAS diff2 = leftParent - rightParent;
        auto[lt_c2, eq_c2, gt_c2] = cdpf2.compare(tio, yield, diff2, tio.aes_ops());
        
        STATS_POST("[SEGTREE][BITVEC] cdpf2.compare (siblings) Stats (level=" + std::to_string(level) + ")");

        // --- mpc_or for isDone ---
        STATS_PRE();
       
        mpc_or(tio, yield, isDone, eq_c2, isDone);
        
        STATS_POST("[SEGTREE][BITVEC] mpc_or (isDone) Stats (level=" + std::to_string(level) + ")");

        // --- CDPF compare for range ---
        STATS_PRE();
        
        CDPF cdpf = tio.cdpf(yield);
        RegAS diff = rightLevelIndex - leftLevelIndex;
        auto[lt_c, eq_c, gt_c] = cdpf.compare(tio, yield, diff, tio.aes_ops());
        
        STATS_POST("[SEGTREE][BITVEC] cdpf.compare (range) Stats (level=" + std::to_string(level) + ")");

        // --- mpc_or for valid ---
        STATS_PRE();
        
        RegBS valid;
        mpc_or(tio, yield, valid, eq_c, gt_c);
        
        STATS_POST("[SEGTREE][BITVEC] mpc_or (valid) Stats (level=" + std::to_string(level) + ")");

        // --- Prepare inclusion without extra isEven test (already zeroed sibling indices) ---
        STATS_PRE();
        RegXS leftSiblingIncluded; // will become 'incl' if Check passes later
        leftSiblingIncluded.set(0);
        RegXS rightSiblingIncluded;
        rightSiblingIncluded.set(0);
        STATS_POST("[SEGTREE][BITVEC] Init siblingIncluded placeholders Stats (level=" + std::to_string(level) + ")");

        // --- mpc_and for Check ---
        STATS_PRE();
        
        RegBS Check;
        mpc_and(tio, yield, Check, one ^ isDone, valid);
        
        STATS_POST("[SEGTREE][BITVEC] mpc_and (Check) Stats (level=" + std::to_string(level) + ")");

        // --- mpc_select for leftSiblingIncluded (Check) ---
        STATS_PRE();
        
        mpc_select(tio, yield, leftSiblingIncluded, Check, leftSiblingIncluded, incl);
        
        STATS_POST("[SEGTREE][BITVEC] mpc_select (leftSiblingIncluded, Check) Stats (level=" + std::to_string(level) + ")");

        // --- mpc_select for rightSiblingIncluded (Check) ---
        STATS_PRE();
        
        mpc_select(tio, yield, rightSiblingIncluded, Check, rightSiblingIncluded, incl);
        
        STATS_POST("[SEGTREE][BITVEC] mpc_select (rightSiblingIncluded, Check) Stats (level=" + std::to_string(level) + ")");

        // --- Set bitVec entries ---
        STATS_PRE();

        bitVecLevel[leftSibling] = leftSiblingIncluded;

        STATS_POST("[SEGTREE][BITVEC] Set bitVecLevel[leftSibling] entries Stats (level=" + std::to_string(level) + ")");

        // --- Set bitVec entries ---
        STATS_PRE();

        bitVecLevel[rightSibling] = rightSiblingIncluded;

        STATS_POST("[SEGTREE][BITVEC] Set bitVecLevel[rightSibling] entries Stats (level=" + std::to_string(level) + ")");

        if(i == 1)
        {   
            STATS_PRE();

            bitVecLevel[leftLevelIndex] = incl;
            bitVecLevel[rightLevelIndex] = incl;

            STATS_POST("[SEGTREE][BITVEC] Set bitVecLevel[left/right] entries Stats (level=" + std::to_string(level) + ")");
        }
        
        #ifdef SEGTREE_VERBOSE
        // value_t recons = mpc_reconstruct(tio, yield, valid);
        auto leftIndRecons = mpc_reconstruct(tio, yield, leftLevelIndex);
        auto rightIndRecons = mpc_reconstruct(tio, yield, rightLevelIndex);
        auto leftSiblingIndRecons = mpc_reconstruct(tio, yield, leftSibling);
        auto rightSiblingIndRecons = mpc_reconstruct(tio, yield, rightSibling);
        auto lIsIncluded = mpc_reconstruct(tio, yield, bitVecLevel[leftLevelIndex]);
        auto rIsIncluded = mpc_reconstruct(tio, yield, bitVecLevel[rightLevelIndex]);
        auto lSiblingIncluded = mpc_reconstruct(tio, yield, bitVecLevel[leftSibling]);
        auto rSiblingIncluded = mpc_reconstruct(tio, yield, bitVecLevel[rightSibling]);
        auto leftParentReccons = mpc_reconstruct(tio, yield, leftParent);
        auto rightParentReccons = mpc_reconstruct(tio, yield, rightParent);
        // auto isDoneRecon = mpc_reconstruct(tio, yield, isDone);
        // auto isLleftchild = mpc_reconstruct(tio, yield, isL_leftchild);
        // auto isRrightchild = mpc_reconstruct(tio, yield, isR_rightchild);
        // std::cout << "Level: " << level << " Left Node [" << leftIndRecons << "] (bitVec[ "<< ((1ULL << level) + leftIndRecons) << "]) isincluded: "  << lIsIncluded << " Right Node [" << rightIndRecons << "] (bitVec[" << ((1ULL << level) + rightIndRecons) << "]) isincluded: "  << rIsIncluded << " valid: " << recons << std::endl;
        std::cout << " Level: " << level << " [" << leftIndRecons << "," << rightIndRecons <<  "]" << std::endl;
        std::cout << " leftParent: " << leftParentReccons << " rightParent: " << rightParentReccons << std::endl;
        std::cout << " leftSibling: " << leftSiblingIndRecons << " rightSibling: " << rightSiblingIndRecons << std::endl;
        // std::cout << " isDone: " << isDoneRecon << " valid: " << recons << std::endl;
        // std::cout << " isLleftchild: " << isLleftchild << " isRrightchild: " << isRrightchild << std::endl;
        std::cout << " bitVec["<< (levelStart + leftIndRecons) << "] isincluded: "  << lIsIncluded << " bitVec[" << (levelStart + rightIndRecons) << "] isincluded: "  << rIsIncluded << std::endl;
        std::cout << " bitVec["<< (levelStart + leftSiblingIndRecons) << "] isincluded: "  << lSiblingIncluded << " bitVec[" << (levelStart + rightSiblingIndRecons) << "] isincluded: "  << rSiblingIncluded << std::endl;
        // std::cout << " bitVec[2]: " << mpc_reconstruct(tio, yield, bitVecArray[2]) << std::endl;
        #endif

        leftLevelIndex = leftParent;
        rightLevelIndex = rightParent;
    }
}

void SegmentTree6::RangeSum(MPCTIO &tio,  MPCIO &mpcio, yield_t & yield, RegAS left, RegAS right) {
    Duoram < RegXS > bitVec(tio.player(), num_items);
    getBitVector(tio, mpcio, yield, bitVec, left, right);

    auto bitVecArray = bitVec.flat(tio, yield);
    auto SegTreeArray = oram.flat(tio, yield);

    #ifdef SEGTREE_VERBOSE
    std::cout << "Bit Vector:" << std::endl;
    auto bitVecRecons = bitVecArray.reconstruct();
    for(size_t i=0; i<num_items; i++) {
        std::cout << "bitVecArray[" << i << "] = " << bitVecRecons[i].share() << std::endl;
    }
    #endif

    RegAS sum;
    sum.set(0);

    // --- Measure: RangeSum accumulation loop (level-wise) ---
    STATS_PRE();
    for(size_t level = 0; level < depth; level++) {
        size_t levelStart = getLevelStart6(level);
        size_t levelLength = getLevelLength6(level);

        typename Duoram<RegXS>::Flat bitVecLevel(bitVecArray, tio, yield, levelStart, levelLength);
        typename Duoram<RegAS>::Flat segTreeLevel(SegTreeArray, tio, yield, levelStart, levelLength);
        
        // Iterate from 1 to avoid the extra node (position 0)
        for(size_t pos = 1; pos < levelLength; pos++) {
            
            RegXS element = bitVecLevel[pos];
            RegBS incl = element.bitat(0);
            RegAS val = segTreeLevel[pos];
            RegAS zero;
            zero.set(0);

            RegAS sum1;
            mpc_select(tio, yield, sum1, incl, zero, val);

            sum.ashare += sum1.ashare;
        }
    }
    STATS_POST("[SEGTREE][RANGESUM] Accumulation Loop Stats");

    #ifdef SEGTREE_VERBOSE
    value_t answer = mpc_reconstruct(tio, yield, sum);
    std::cout << "Sum = " << answer << std::endl;
    #endif
}

void SegmentTree6::Update(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, RegAS index, RegAS value) {
    auto SegTreeArray = oram.flat(tio, yield);
    auto parentArray = parent.flat(tio, yield);

    // Adjust index to be in valid range for leaf level and convert to global index
    size_t leafLevel = depth - 1;
    size_t leafStart = getLevelStart6(leafLevel);

    // Access last level Flat object and get current value
    size_t leafLength = getLevelLength6(leafLevel);
    typename Duoram<RegAS>::Flat leafLevel_flat(SegTreeArray, tio, yield, leafStart, leafLength);
    
    // --- Measure: Leaf access (read current value) ---
    STATS_PRE();
    RegAS currVal = leafLevel_flat[index];
    STATS_POST("[SEGTREE][UPDATE] Leaf Read Stats (leafLevel[adjustedIndex])");
    RegAS diff = value - currVal;

    #ifdef SEGTREE_VERBOSE
    auto recons_index = mpc_reconstruct(tio, yield, index);
    auto recons_currVal = mpc_reconstruct(tio, yield, currVal);
    auto recons_newvalue = mpc_reconstruct(tio, yield, value);
    auto recons_diff = mpc_reconstruct(tio, yield, diff);
    std::cout << "Index to be updated in the original array = " << recons_index << std::endl;
    std::cout << "Index to be updated in the segment tree array = " << leafStart + recons_index << std::endl;
    std::cout << "Current Value at index = " << recons_currVal << std::endl;
    std::cout << "New Value to be updated = " << recons_newvalue << std::endl;
    std::cout << "Diff = " << recons_diff << std::endl;
    #endif
    
    for(size_t i = 1; i <= depth; i++) {
        size_t level = depth - i;
        size_t levelStart = getLevelStart6(level);
        size_t levelLength = getLevelLength6(level);

        typename Duoram<RegAS>::Flat parentLevel(parentArray, tio, yield, levelStart, levelLength);
        typename Duoram<RegAS>::Flat segTreeLevel(SegTreeArray, tio, yield, levelStart, levelLength);
    
        STATS_PRE();
        segTreeLevel[index] += diff;
        STATS_POST("[SEGTREE][UPDATE] Level Update Stats (level=" + std::to_string(level) + ") segTreeLevel[index] += diff");

        #ifdef SEGTREE_VERBOSE
        auto recons_Index = mpc_reconstruct(tio, yield, index);
        auto recons_updated = mpc_reconstruct(tio, yield, segTreeLevel[index]);
        std::cout << "Updated Index = " << (levelStart + recons_Index) << " with value = " << recons_updated << std::endl;
        #endif

        index = parentLevel[index];
    }
}


void SegTree6(MPCIO &mpcio, const PRACOptions &opts, char **args) {
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

    // Calculate new array size: 2^d + d - 1
    address_t len = (1<<depth) + depth - 1;

    MPCTIO tio(mpcio, 0, opts.num_cpu_threads);

    run_coroutines(tio, [&tio, &mpcio, len, depth, n_updates, n_queries] (yield_t &yield) {
        
        SegmentTree6 segTree(tio.player(), len, depth);
        segTree.init(tio, yield);
        std::cout << "===== Segment Tree Init Stats =====" << std::endl;
        std::cout << "Updates: " << n_updates << ", Queries: " << n_queries << std::endl;
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        mpcio.reset_stats();
        tio.reset_lamport();

        #ifdef SEGTREE_VERBOSE
        // Print initial segment tree
        std::cout << "\n===== Initial Segment Tree =====" << std::endl;
        segTree.printSegmentTree(tio, yield);
        #endif

        // Perform updates
        for (size_t u = 0; u < n_updates; ++u) {
            std::cout << "\n===== Update " << (u + 1) << " begins =====" << std::endl;
            RegAS index;
            size_t idx_val = (u % (1 << (depth - 1))) + 1; // Cycle through leaf indices (1 to 2^(d-1))
            index.set(tio.player() == 0 ? idx_val : 0); // Cycle through leaf indices


            RegAS value;
            size_t val_to_set = (u + 1) * 50;
            value.set(tio.player() == 0 ? val_to_set : 0); // Use different values for each update

            segTree.Update(tio, mpcio, yield, index, value);
            std::cout << "Update " << (u + 1) << " ends" << std::endl;
        }

        std::cout << "===== Updates Stats =====" << std::endl;
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        mpcio.reset_stats();
        tio.reset_lamport();

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
            size_t max_leaf_idx = (1 << (depth - 1));
            size_t left_val = (q % (1 << (depth - 1))) + 1;  // Start anywhere in valid range (1 to 2^(d-1))
            size_t right_val = left_val + (q % (max_leaf_idx + 1 - left_val)); // Ensuring right > left

            left_index.set(tio.player() == 0 ? left_val : 0); // Start of range
            right_index.set(tio.player() == 0 ? right_val : 0); // End of range

            #ifdef SEGTREE_VERBOSE
            auto recons_left = mpc_reconstruct(tio, yield, left_index);
            auto recons_right = mpc_reconstruct(tio, yield, right_index);
            std::cout << "Range Sum Query [" << recons_left << ", " << recons_right << "]" << std::endl;
            #endif

            segTree.RangeSum(tio, mpcio, yield, left_index, right_index);
            std::cout << "Range Sum Query " << (q + 1) << " ends" << std::endl;
        }

        std::cout << "===== Range Sum Stats =====" << std::endl;
    });
}