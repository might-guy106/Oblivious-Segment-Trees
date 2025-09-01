#include <functional>
#include "mpcops.hpp"
#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include "rdpf.hpp"
#include "shapes.hpp"
#include "segmentTree.hpp"

#define SEGTREE_VERBOSE

/*  
    1 both parties have additive shares of l and r
    2 always include them
    3 if l is left child then add its sibling to the sum or if r is right child add its sibling to sum
    3.1 l and r can be known left child or right child using isEven array computer already (if even left)
    3.2 get the additive shares of the sibling which was also present (as we will be initialising it in init function similar to isEven) .
    3.3 based on above two conditions (and the isdone and isvaild bit) set the sibling index in bitvector.
    4 if sibling of l is r then set isdone bit to 1

*/
    
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

    auto siblingArray = sibling.flat(tio, yield);
    siblingArray.init([this] (size_t i) -> size_t {
        if (i >= 1 && i < num_items) {
            return arrayIndexToLevelPos((i % 2 == 0) ? i + 1 : i - 1); // we took 1s sibling as 0 and 0s as 1 which should not cause any issue.
        } else {
            return size_t(0);
        }
    });

    auto parentArray = parent.flat(tio, yield);
    parentArray.init([this] (size_t i) -> size_t {
        if (i >= 1 && i < num_items) {
            return arrayIndexToLevelPos(i/2);
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
    auto siblingArray = sibling.flat(tio, yield);
    auto parentArray = parent.flat(tio, yield);

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
        typename Duoram < RegAS > ::Flat siblingLevel(siblingArray, tio, yield, (1ULL << level), (1ULL << level)+1);
        typename Duoram < RegAS > ::Flat parentLevel(parentArray, tio, yield, (1ULL << level), (1ULL << level)+1);

        RegAS leftParent = parentLevel[left];
        RegAS rightParent = parentLevel[right];
        RegAS leftSibling = siblingLevel[left]; 
        RegAS rightSibling = siblingLevel[right];
        
        // if l and r are siblings it is already done before this iteration
        CDPF cdpf2 = tio.cdpf(yield);
        RegAS diff2 = leftParent - rightParent;
        auto[lt_c2, eq_c2, gt_c2] = cdpf2.compare(tio, yield, diff2, tio.aes_ops());
        mpc_or(tio, yield, isDone, eq_c2, isDone);
        
        CDPF cdpf = tio.cdpf(yield);
        RegAS diff = right - left;
        auto[lt_c, eq_c, gt_c] = cdpf.compare(tio, yield, diff, tio.aes_ops());
        RegBS valid;
        mpc_or(tio, yield, valid, eq_c, gt_c);  


        RegXS leftSiblingIncluded = bitVecLevel[leftSibling];
        RegXS rightSiblingIncluded = bitVecLevel[rightSibling];

        RegXS isEvenL = isEvenLevel[left];
        RegBS isL_leftchild = isEvenL.bitat(0);
        mpc_select(tio, yield, leftSiblingIncluded, isL_leftchild, leftSiblingIncluded, incl);
        
        RegXS isEvenR = isEvenLevel[right];
        RegBS isR_rightchild = one ^ isEvenR.bitat(0);
        mpc_select(tio, yield, rightSiblingIncluded, isR_rightchild, rightSiblingIncluded, incl);

        // checks if not done and is valid
        RegBS Check;
        mpc_and(tio, yield, Check, one ^ isDone, valid);

        mpc_select(tio, yield, leftSiblingIncluded, Check, excl, leftSiblingIncluded);
        mpc_select(tio, yield, rightSiblingIncluded, Check, excl, rightSiblingIncluded);

        bitVecLevel[leftSibling] = leftSiblingIncluded;
        bitVecLevel[rightSibling] = rightSiblingIncluded;

        if(i == 1)
        {
            bitVecLevel[left] = incl;
            bitVecLevel[right] = incl;
        }
        
        #ifdef SEGTREE_VERBOSE
        // value_t recons = mpc_reconstruct(tio, yield, valid);
        auto leftIndRecons = mpc_reconstruct(tio, yield, left);
        auto rightIndRecons = mpc_reconstruct(tio, yield, right);
        // auto leftSiblingIndRecons = mpc_reconstruct(tio, yield, leftSibling);
        // auto rightSiblingIndRecons = mpc_reconstruct(tio, yield, rightSibling);
        // auto lIsIncluded = mpc_reconstruct(tio, yield, bitVecLevel[left]);
        // auto rIsIncluded = mpc_reconstruct(tio, yield, bitVecLevel[right]);
        // auto lSiblingIncluded = mpc_reconstruct(tio, yield, bitVecLevel[leftSibling]);
        // auto rSiblingIncluded = mpc_reconstruct(tio, yield, bitVecLevel[rightSibling]);
        // auto isDoneRecon = mpc_reconstruct(tio, yield, isDone);
        // auto isLleftchild = mpc_reconstruct(tio, yield, isL_leftchild);
        // auto isRrightchild = mpc_reconstruct(tio, yield, isR_rightchild);
        // std::cout << "Level: " << level << " Left Node [" << leftIndRecons << "] (bitVec[ "<< ((1ULL << level) + leftIndRecons) << "]) isincluded: "  << lIsIncluded << " Right Node [" << rightIndRecons << "] (bitVec[" << ((1ULL << level) + rightIndRecons) << "]) isincluded: "  << rIsIncluded << " valid: " << recons << std::endl;
        std::cout << " Level: " << level << " [" << leftIndRecons << "," << rightIndRecons <<  "]" << std::endl;
        // std::cout << " isDone: " << isDoneRecon << " valid: " << recons << std::endl;
        // std::cout << " isLleftchild: " << isLleftchild << " isRrightchild: " << isRrightchild << std::endl;
        // std::cout << " bitVec["<< ((1ULL << level) + leftIndRecons) << "] isincluded: "  << lIsIncluded << " bitVec[" << ((1ULL << level) + rightIndRecons) << "] isincluded: "  << rIsIncluded << std::endl;
        // std::cout << " bitVec["<< ((1ULL << level) + leftSiblingIndRecons) << "] isincluded: "  << lSiblingIncluded << " bitVec[" << ((1ULL << level) + rightSiblingIndRecons) << "] isincluded: "  << rSiblingIncluded << std::endl;
        // std::cout << " bitVec[2]: " << mpc_reconstruct(tio, yield, bitVecArray[2]) << std::endl;
        #endif

        left = leftParent;
        right = rightParent;
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

    // SegTreeArray[index1] = value;
    for(size_t i=1; i<=depth; i++) {
        size_t level = depth - i;
        typename Duoram < RegAS > ::Flat parentLevel(parentArray, tio, yield, (1ULL << level), (1ULL << level)+1);
        typename Duoram < RegAS > ::Flat segTreeLevel(SegTreeArray, tio, yield, (1ULL << level), (1ULL << level)+1);

        segTreeLevel[index] += diff;
        
        #ifdef SEGTREE_VERBOSE
        auto recons_Index = mpc_reconstruct(tio, yield, index);
        auto recons_updated = mpc_reconstruct(tio, yield, segTreeLevel[index]);
        std::cout << "Updated Index = " << ((1ULL << (level-1)) + recons_Index) << " with value = " << recons_updated << std::endl;
        #endif

        RegAS parentIndex = parentLevel[index];
        index = parentIndex;
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
