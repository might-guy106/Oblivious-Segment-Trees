#include <functional>
#include "mpcops.hpp"
#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include "rdpf.hpp"
#include "shapes.hpp"
#include "segmentTree2.hpp"

// #define SEGTREE_VERBOSE
// To enable timing/stat instrumentation define SEGTREE_VERBOSE2
// #define SEGTREE_VERBOSE2

#ifdef SEGTREE_VERBOSE2
#define STATS_PRE() do { tio.sync_lamport(); mpcio.reset_stats(); tio.reset_lamport(); } while(0)
#define STATS_POST(MSG) do { tio.sync_lamport(); std::cout << MSG << std::endl; mpcio.dump_stats(std::cout); mpcio.reset_stats(); tio.reset_lamport(); } while(0)
#else
#define STATS_PRE() do {} while(0)
#define STATS_POST(MSG) do {} while(0)
#endif

/*  
    1 both parties have additive shares of l and r
    2 always include them
    3 if l is left child then add its sibling to the sum or if r is right child add its sibling to sum
    3.1 l and r can be known left child or right child using isEven array computer already (if even left)
    3.2 get the additive shares of the sibling which was also present (as we will be initialising it in init function similar to isEven) .
    3.3 based on above two conditions (and the isdone and isvaild bit) set the sibling index in bitvector.
    4 if sibling of l is r then set isdone bit to 1

*/

void SegmentTree2::init(MPCTIO &tio, yield_t & yield) {
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
            return (i % 2 == 0) ? (i + 1) : (i - 1); // direct global sibling index
        } else {
            return size_t(0);
        }
    });

    auto parentArray = parent.flat(tio, yield);
    parentArray.init([this] (size_t i) -> size_t {
        if (i >= 1 && i < num_items) {
            return i/2; // direct global parent index
        } else {
            return size_t(0);
        }
    });
}

void SegmentTree2::printSegmentTree(MPCTIO &tio, yield_t & yield) {
    auto SegTreeArray = oram.flat(tio, yield);
    auto SegTreeRecons = SegTreeArray.reconstruct();
    for(size_t i=1; i<num_items; i++) {
        std::cout << "SegTreeArray[" << i << "] = " << SegTreeRecons[i].share() << std::endl;
    }
}

void SegmentTree2::getBitVector(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, Duoram < RegXS > &bitVec, RegAS left, RegAS right) {

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
    RegBS one;
    one.set(tio.player()==0 ? 1 : 0);
    RegBS zero;
    zero.set(0);

    for(uint32_t i=1; i<=height; i++) {
        // level variable retained only for logging clarity
        uint32_t level = height - i;
        #ifndef SEGTREE_VERBOSE2
            (void)level; // silence unused warning when stats disabled
        #endif

        // --- ORAM base reads: parents + siblings + isEven (global indices) ---
        STATS_PRE();
        RegAS leftParent = parentArray[left];
        RegAS rightParent = parentArray[right];
        RegAS leftSibling = siblingArray[left];
        RegAS rightSibling = siblingArray[right];
        RegXS isEvenL = isEvenArray[left];
        RegXS isEvenR = isEvenArray[right];
        STATS_POST("[SEGTREE][BITVEC] ORAM Reads (parents+siblings+isEven) Stats (level=" + std::to_string(level) + ")");

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
        RegAS diff = right - left;
        auto[lt_c, eq_c, gt_c] = cdpf.compare(tio, yield, diff, tio.aes_ops());
        
        STATS_POST("[SEGTREE][BITVEC] cdpf.compare (range) Stats (level=" + std::to_string(level) + ")");

        // --- mpc_or for valid ---
        STATS_PRE();
        
        RegBS valid;
        mpc_or(tio, yield, valid, eq_c, gt_c);
        
        STATS_POST("[SEGTREE][BITVEC] mpc_or (valid) Stats (level=" + std::to_string(level) + ")");

        // --- mpc_select for leftSiblingIncluded ---
        STATS_PRE();
        
        RegXS leftSiblingIncluded;
        RegBS isL_leftchild = isEvenL.bitat(0);
        mpc_select(tio, yield, leftSiblingIncluded, isL_leftchild, leftSiblingIncluded, incl);
        
        STATS_POST("[SEGTREE][BITVEC] mpc_select (leftSiblingIncluded) Stats (level=" + std::to_string(level) + ")");

        // --- mpc_select for rightSiblingIncluded ---
        STATS_PRE();
        
        RegXS rightSiblingIncluded;
        RegBS isR_rightchild = one ^ isEvenR.bitat(0);
        mpc_select(tio, yield, rightSiblingIncluded, isR_rightchild, rightSiblingIncluded, incl);
        
        STATS_POST("[SEGTREE][BITVEC] mpc_select (rightSiblingIncluded) Stats (level=" + std::to_string(level) + ")");

        // --- mpc_and for Check ---
        STATS_PRE();
        
        RegBS Check;
        mpc_and(tio, yield, Check, one ^ isDone, valid);
        
        STATS_POST("[SEGTREE][BITVEC] mpc_and (Check) Stats (level=" + std::to_string(level) + ")");

        // --- mpc_select for leftSiblingIncluded (Check) ---
        STATS_PRE();
        
        mpc_select(tio, yield, leftSiblingIncluded, Check, excl, leftSiblingIncluded);
        
        STATS_POST("[SEGTREE][BITVEC] mpc_select (leftSiblingIncluded, Check) Stats (level=" + std::to_string(level) + ")");

        // --- mpc_select for rightSiblingIncluded (Check) ---
        STATS_PRE();
        
        mpc_select(tio, yield, rightSiblingIncluded, Check, excl, rightSiblingIncluded);
        
        STATS_POST("[SEGTREE][BITVEC] mpc_select (rightSiblingIncluded, Check) Stats (level=" + std::to_string(level) + ")");

        // --- Set bitVec entries ---
        STATS_PRE();

        bitVecArray[leftSibling] = leftSiblingIncluded;

        STATS_POST("[SEGTREE][BITVEC] Set bitVecLevel[leftSibling] entries Stats (level=" + std::to_string(level) + ")");

        // --- Set bitVec entries ---
        STATS_PRE();

        bitVecArray[rightSibling] = rightSiblingIncluded;

        STATS_POST("[SEGTREE][BITVEC] Set bitVecLevel[rightSibling] entries Stats (level=" + std::to_string(level) + ")");

        if(i == 1)
        {   
            STATS_PRE();
            
            bitVecArray[left] = incl;
            bitVecArray[right] = incl;

            STATS_POST("[SEGTREE][BITVEC] Set bitVecLevel[left/right] entries Stats (level=" + std::to_string(level) + ")");
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

void SegmentTree2::RangeSum(MPCTIO &tio,  MPCIO &mpcio, yield_t & yield, RegAS left, RegAS right) {
    
    RegAS disp;
    disp.set(tio.player() == 0 ? (1 << (depth - 1)) : 0);

    left = left + disp;
    right = right + disp;

    Duoram < RegXS > bitVec(tio.player(), num_items);
    getBitVector(tio, mpcio, yield, bitVec, left, right);

    auto bitVecArray = bitVec.flat(tio, yield);
    auto SegTreeArray = oram.flat(tio, yield);


    RegAS sum;
    sum.set(0);

    // --- Measure: RangeSum accumulation loop ---
    STATS_PRE();
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
    STATS_POST("[SEGTREE][RANGESUM] Accumulation Loop Stats");

    #ifdef SEGTREE_VERBOSE
    value_t answer = mpc_reconstruct(tio, yield, sum);
    std::cout << "Sum = " << answer << std::endl;
    #endif
}

void SegmentTree2::Update(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, RegAS index, RegAS value) {
    auto SegTreeArray = oram.flat(tio, yield);
    auto parentArray = parent.flat(tio, yield);

    RegAS disp;
    disp.set(tio.player() == 0 ? (1 << (depth - 1)) : 0);

    RegAS indexSegArr = index + disp;

    // --- Measure: Leaf access (read current value) ---
    STATS_PRE();
    RegAS currVal = SegTreeArray[indexSegArr];
    STATS_POST("[SEGTREE][UPDATE] Leaf Read Stats (SegTreeArray[indexSegArr])");
    RegAS diff = value - currVal;

    #ifdef SEGTREE_VERBOSE
    auto recons_index = mpc_reconstruct(tio, yield, index);
    auto recons_index1 = mpc_reconstruct(tio, yield, indexSegArr);
    auto recons_currVal = mpc_reconstruct(tio, yield, currVal);
    auto recons_newvalue = mpc_reconstruct(tio, yield, value);
    auto recons_diff = mpc_reconstruct(tio, yield, diff);
    std::cout << "Index to be updated in the original array = " << recons_index << std::endl;
    std::cout << "Index to be updated in the segment tree array = " << recons_index1 << std::endl;
    std::cout << "Current Value at index = " << recons_currVal << std::endl;
    std::cout << "New Value to be updated = " << recons_newvalue << std::endl;
    std::cout << "Diff = " << recons_diff << std::endl;
    #endif

    SegTreeArray[indexSegArr] += diff;
    for(size_t i=1; i<=depth-1; i++) {
        size_t level = depth - i;
        #ifndef SEGTREE_VERBOSE2
            (void)level;
        #endif

        // --- Measure: Parent access (parentLevel[index]) ---
        STATS_PRE();
        RegAS parentIndex = parentArray[indexSegArr];
        STATS_POST("[SEGTREE][UPDATE] Parent Read Stats (level=" + std::to_string(level) + ") parentLevel[index]");
        
        // --- Measure: Per-level update write (SegTreeArray[index] += diff) ---
        STATS_PRE();
        SegTreeArray[parentIndex] += diff;
        STATS_POST("[SEGTREE][UPDATE] Level Update Stats (level=" + std::to_string(level) + ") SegTreeArray[index] += diff");
        
        indexSegArr = parentIndex;

        #ifdef SEGTREE_VERBOSE
        auto recons_Index = mpc_reconstruct(tio, yield, indexSegArr);
        auto recons_updated = mpc_reconstruct(tio, yield, SegTreeArray[indexSegArr]);
        std::cout << "Updated Index = " << (recons_Index) << " with value = " << recons_updated << std::endl;
        #endif
    }
}


void SegTree2(MPCIO &mpcio, const PRACOptions &opts, char **args) {
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

    run_coroutines(tio, [&tio, &mpcio, len, depth, n_updates, n_queries] (yield_t &yield) {
        
        SegmentTree2 segTree(tio.player(), len, depth);
        segTree.init(tio, yield);
        std::cout << "===== Segment Tree 2 Init Stats =====" << std::endl;
        std::cout << "Depth: " << depth << ", Size: " << len << std::endl;
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
            size_t idx_val = u % (1 << (depth - 1));
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
            size_t max_leaf_idx = (1 << (depth - 1)) - 1;
            size_t left_val = q % (1 << (depth - 1));  // Start anywhere in valid range
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
