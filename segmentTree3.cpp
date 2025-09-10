#include <functional>
#include "mpcops.hpp"
#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include "rdpf.hpp"
#include "shapes.hpp"
#include "segmentTree3.hpp"

// uncomment to enable intermediate reconstructions and logging
#define SEGTREE_VERBOSE

/*
The segment tree data structure is stored in an array using a standard binary tree layout with starting index as 1 (and not 0).
For nodes stored at index i of the array:
- The parent is stored at i/2
- The left child is stored at 2*i  
- The right child is stored at 2*i + 1

Each internal node stores the sum of values in its subtree range.
Leaf nodes store the actual array values.

Example: For an array [100, 200, 300, 400] (4 elements), the segment tree structure is:

                                1000
                               /    \
                             300    700
                            /  \   /  \
                          100  200 300  400

The segment tree is stored in array format as:
Index:  [ 0 |  1  |  2  |  3  |  4  |  5  |  6  |  7  ]
Value:  [ 0 | 1000| 300 | 700 | 100 | 200 | 300 | 400 ]

Where:
- Index 1: Root (sum of entire array = 1000)
- Index 2: Sum of left half [100,200] = 300  
- Index 3: Sum of right half [300,400] = 700
- Indices 4-7: Leaf values [100, 200, 300, 400]

For depth d, the segment tree has:
- 2^(d-1) leaf nodes starting at index 2^(d-1)
- Total size: 2^d nodes
- Tree height: d levels (0-indexed from root)
*/

/*
RANGE SUM QUERY OPERATION

The range sum query efficiently computes the sum of elements in a given range [left, right] 

Algorithm Overview:

1. Starting from the leaf level with the range endpoints
2. Moving up the tree level by level  
3. At each level, including sibling nodes that fall within the range. if they dont fall mark the 0th node
4. Stopping when the range endpoints have same parent

Example: Range Sum Query [1, 2] on array [100, 200, 300, 400]

Initial tree:
                                1000
                               /    \
                             300    700
                            /  \   /  \
                          100  200 300  400
                          [0]  [1] [2] [3]    <- Array indices
                           4    5   6   7     <- Tree indices

Step 1: Convert array indices to tree indices
- left = 1 → tree index 5 (200)
- right = 2 → tree index 6 (300)

Step 2: Process leaf level (level 2)
- Check siblings: 
* Node 5 is right child, sibling 4 (100) is outside range
* Node 6 is left child, sibling 7 (400) is outside range
- Include both endpoints: nodes 5 and 6
- BitVector: [0,0,0,0,0,1,1,0] (include indices 5,6)

- left parent = 5/2 = 2, right parent = 6/2 = 3

Step 3: Move to parent level (level 1)  
- Since Parents of 2 and 3 is same that is 1 set isDone bit to 1
- as isDone is set we dont mark any more nodes from here
- BitVector: [0,0,0,0,0,1,1,0] (no additional siblings to include)

- left parent = 2/2 = 1, right parent = 3/2 = 1

Step 4: Move to parent level (level 0)
- Since isDone we dont mark any more nodes 
- BitVector: [0,0,0,0,0,1,1,0] (no additional siblings to include)


Final BitVector indicates which nodes to sum: nodes 5 (200) + 6 (300) = 500
*/

/*
UPDATE OPERATION

Algorithm Overview:
1. Calculate the difference between new value and current value
2. Update the target leaf node with this difference
3. Propagate the difference up to all ancestors along the path to root
4. Each ancestor's sum increases/decreases by the same difference

Example: Update array[1] from 200 to 250 (difference = +50)

Before update:
                                1000
                               /    \
                             300    700
                            /  \   /  \
                          100  200 300  400
                          [0]  [1] [2] [3]    <- Array indices
                           4    5   6   7     <- Tree indices

After update:
                                1050  ← +50
                               /    \
                             350    700  ← +50
                            /  \   /  \
                          100  250 300  400  ← +50
                          [0]  [1] [2] [3]  <- Array indices
                           4    5   6   7   <- Tree indices

The update process follows the parent chain:
- Leaf index 5 → parent 5/2 = 2 → parent 2/2 = 1 (root)
*/

void SegmentTree3::init(MPCTIO &tio, yield_t & yield) {
    
    // Create a flat reference to the main segment tree ORAM
    auto SegTreeArray = oram.flat(tio, yield);
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
        segTree[i] = segTree[2*i] + segTree[2*i + 1];
    }

    // Initialize the ORAM segment tree with computed values
    SegTreeArray.init([this, &segTree] (size_t i) -> size_t {
        if (i >= 1 && i < num_items) {
            return segTree[i];
        } else {
            return size_t(0x7fffffffffffffff);
        }
    });

    // Initialize left child sibling array
    // For each node, store its right sibling if it's a left child (even index), else 0
    auto leftChildSiblingArray = leftChildSibling.flat(tio, yield);
    leftChildSiblingArray.init([this] (size_t i) -> size_t {
        if (i >= 1 && i < num_items) {
            return (i % 2 == 0) ? (i + 1) : 0; // only even (left) keep sibling else 0
        } else { return size_t(0); }
    });

    // Initialize right child sibling array  
    // For each node, store its left sibling if it's a right child (odd index), else 0
    auto rightChildSiblingArray = rightChildSibling.flat(tio, yield);
    rightChildSiblingArray.init([this] (size_t i) -> size_t {
        if (i >= 1 && i < num_items) {
            return (i % 2 == 1) ? (i - 1) : 0; // only odd (right) keep sibling else 0
        } else { return size_t(0); }
    });

    // Initialize parent array
    // For each node, store its parent's index using standard binary tree parent formula
    auto parentArray = parent.flat(tio, yield);
    parentArray.init([this] (size_t i) -> size_t {
        if (i >= 1 && i < num_items) {
            return i/2; // direct global parent index
        } else {
            return size_t(0);
        }
    });
}

void SegmentTree3::printSegmentTree(MPCTIO &tio, yield_t & yield) {
    auto SegTreeArray = oram.flat(tio, yield);
    auto SegTreeRecons = SegTreeArray.reconstruct();
    for(size_t i=1; i<num_items; i++) {
        std::cout << "SegTreeArray[" << i << "] = " << SegTreeRecons[i].share() << std::endl;
    }
}

void SegmentTree3::getBitVector(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, Duoram < RegXS > &bitVec, RegAS left, RegAS right) {

    // Get flat references to all helper arrays
    auto bitVecArray = bitVec.flat(tio, yield);
    auto leftChildSiblingArray = leftChildSibling.flat(tio, yield);
    auto rightChildSiblingArray = rightChildSibling.flat(tio, yield);
    auto parentArray = parent.flat(tio, yield);

    // Calculate tree height for iteration
    uint32_t height = static_cast<uint32_t>(std::log2(num_items));

    // Initialize constants for bit vector operations
    RegXS incl;  // inclusion marker
    incl.set(tio.player()==0 ? 1 : 0);
    RegXS excl;  // exclusion marker
    excl.set(0);
    RegBS isDone;  // tracks if range endpoints have converged
    isDone.set(0);
    RegBS one;
    one.set(tio.player()==0 ? 1 : 0);
    RegBS zero;
    zero.set(0);
    
    // Process each level from leaves to root
    for(uint32_t i=1; i<=height; i++) {
        uint32_t level = height - i;  // current level being processed

        // Read parent indices and conditional sibling indices
        RegAS leftParent = parentArray[left];
        RegAS rightParent = parentArray[right];
        RegAS leftSibling = leftChildSiblingArray[left];   // 0 if left not a left child
        RegAS rightSibling = rightChildSiblingArray[right]; // 0 if right not a right child

        // Check if range endpoints have converged (same parent)
        CDPF cdpf2 = tio.cdpf(yield);
        RegAS diff2 = leftParent - rightParent;
        auto[lt_c2, eq_c2, gt_c2] = cdpf2.compare(tio, yield, diff2, tio.aes_ops());
        
        // Update isDone flag if parents are equal
        mpc_or(tio, yield, isDone, eq_c2, isDone);

        // Check if we have a valid range (right >= left)
        CDPF cdpf = tio.cdpf(yield);
        RegAS diff = right - left;
        auto[lt_c, eq_c, gt_c] = cdpf.compare(tio, yield, diff, tio.aes_ops());
        
        // Range is valid if right >= left
        RegBS valid;
        mpc_or(tio, yield, valid, eq_c, gt_c);

        // Initialize sibling inclusion variables
        RegXS leftSiblingIncluded;  // will become 'incl' if Check passes later
        leftSiblingIncluded.set(0);
        RegXS rightSiblingIncluded;
        rightSiblingIncluded.set(0);

        // Check if we should include siblings: range is valid and not done
        RegBS Check;
        mpc_and(tio, yield, Check, one ^ isDone, valid);
        
        // Conditionally include left sibling if check passes
        mpc_select(tio, yield, leftSiblingIncluded, Check, leftSiblingIncluded, incl);
        
        // Conditionally include right sibling if check passes
        mpc_select(tio, yield, rightSiblingIncluded, Check, rightSiblingIncluded, incl);

        // Set bit vector entries for siblings
        bitVecArray[leftSibling] = leftSiblingIncluded;
        bitVecArray[rightSibling] = rightSiblingIncluded;

        // For the first iteration (leaf level), include the range endpoints
        if(i == 1) {
            bitVecArray[left] = incl;
            bitVecArray[right] = incl;
        }
        
        #ifdef SEGTREE_VERBOSE
        auto leftIndRecons = mpc_reconstruct(tio, yield, left);
        auto rightIndRecons = mpc_reconstruct(tio, yield, right);
        std::cout << " Level: " << level << " [" << leftIndRecons << "," << rightIndRecons <<  "]" << std::endl;
        #endif

        // Move to parent level for next iteration
        left = leftParent;
        right = rightParent;
    }
}

void SegmentTree3::RangeSum(MPCTIO &tio,  MPCIO &mpcio, yield_t & yield, RegAS left, RegAS right) {
    
    // Adjust indices to segment tree array format (add offset for leaf position)
    RegAS disp;
    disp.set(tio.player() == 0 ? (1 << (depth - 1)) : 0);

    left = left + disp;
    right = right + disp;
    
    // Create bit vector to mark which nodes to include in the sum
    Duoram < RegXS > bitVec(tio.player(), num_items);
    getBitVector(tio, mpcio, yield, bitVec, left, right);

    // Get flat references to both arrays
    auto bitVecArray = bitVec.flat(tio, yield);
    auto SegTreeArray = oram.flat(tio, yield);

    // Initialize sum accumulator
    RegAS sum;
    sum.set(0);

    // Accumulate values from all nodes marked in the bit vector
    for(size_t i=1; i<num_items; i++) {
        RegXS element = bitVecArray[i];
        RegBS incl = element.bitat(0);  // check inclusion bit
        RegAS val = SegTreeArray[i];    // get segment tree value
        RegAS zero;
        zero.set(0);

        // Conditionally add value to sum if included
        RegAS sum1;
        mpc_select(tio, yield, sum1, incl, zero, val);

        sum.ashare += sum1.ashare;
    }

    #ifdef SEGTREE_VERBOSE
    value_t answer = mpc_reconstruct(tio, yield, sum);
    std::cout << "Sum = " << answer << std::endl;
    #endif
}

void SegmentTree3::Update(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, RegAS index, RegAS value) {
    
    // Get flat references to the main arrays
    auto SegTreeArray = oram.flat(tio, yield);
    auto parentArray = parent.flat(tio, yield);

    // Calculate offset to convert leaf index to segment tree array index
    RegAS disp;
    disp.set(tio.player() == 0 ? (1 << (depth - 1)) : 0);

    RegAS indexSegArr = index + disp;

    // Read current value at the leaf position
    RegAS currVal = SegTreeArray[indexSegArr];
    RegAS diff = value - currVal;  // calculate difference to propagate

    #ifdef SEGTREE_VERBOSE
    auto recons_index = mpc_reconstruct(tio, yield, index);
    auto recons_indexSegArr = mpc_reconstruct(tio, yield, indexSegArr);
    auto recons_currVal = mpc_reconstruct(tio, yield, currVal);
    auto recons_newvalue = mpc_reconstruct(tio, yield, value);
    auto recons_diff = mpc_reconstruct(tio, yield, diff);
    std::cout << "Index to be updated in the original array = " << recons_index << std::endl;
    std::cout << "Index to be updated in the segment tree array = " << recons_indexSegArr << std::endl;
    std::cout << "Current Value at index = " << recons_currVal << std::endl;
    std::cout << "New Value to be updated = " << recons_newvalue << std::endl;
    std::cout << "Diff = " << recons_diff << std::endl;
    #endif

    // Update the leaf node with the difference
    SegTreeArray[indexSegArr] += diff;
    
    // Propagate the difference up the tree to all ancestors
    for(size_t i=1; i<=depth-1; i++) {
        size_t level = depth - i;

        // Get parent index and update it
        RegAS parentIndex = parentArray[indexSegArr];
        SegTreeArray[parentIndex] += diff;
        
        // Move to the parent for next iteration
        indexSegArr = parentIndex;

        #ifdef SEGTREE_VERBOSE
        auto recons_Index = mpc_reconstruct(tio, yield, indexSegArr);
        auto recons_updated = mpc_reconstruct(tio, yield, SegTreeArray[indexSegArr]);
        std::cout << "Updated Index = " << (recons_Index) << " with value = " << recons_updated << std::endl;
        #endif
    }
}


void SegTree3(MPCIO &mpcio, const PRACOptions &opts, char **args) {
    // Parse command line arguments for configuration
    int nargs = 0;
    while (args[nargs] != nullptr) {
        ++nargs;
    }

    // Default parameters
    nbits_t depth = 5;        // tree depth 
    size_t n_updates = 1;     // number of update operations
    size_t n_queries = 1;     // number of range sum queries

    // Parse command line options
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

    // Calculate segment tree size (standard binary tree layout)
    address_t len = (1<<depth);

    // Initialize MPC communication
    MPCTIO tio(mpcio, 0, opts.num_cpu_threads);

    run_coroutines(tio, [&tio, &mpcio, len, depth, n_updates, n_queries] (yield_t &yield) {
        
        // Create and initialize the segment tree
        SegmentTree3 segTree(tio.player(), len, depth);
        segTree.init(tio, yield);
        
        // Print initialization statistics
        std::cout << "===== Segment Tree 3 Init Stats =====" << std::endl;
        std::cout << "Depth: " << depth << ", Size: " << len << std::endl;
        std::cout << "Updates: " << n_updates << ", Queries: " << n_queries << std::endl;
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        mpcio.reset_stats();
        tio.reset_lamport();

        #ifdef SEGTREE_VERBOSE
        // Print initial segment tree state for debugging
        std::cout << "\n===== Initial Segment Tree =====" << std::endl;
        segTree.printSegmentTree(tio, yield);
        #endif

        // Perform the specified number of update operations
        for (size_t u = 0; u < n_updates; ++u) {
            std::cout << "\n===== Update " << (u + 1) << " begins =====" << std::endl;
            
            // Create update parameters (cycling through leaf indices)
            RegAS index;
            size_t idx_val = u % (1 << (depth - 1));
            index.set(tio.player() == 0 ? idx_val : 0);

            RegAS value;
            size_t val_to_set = (u + 1) * 50;  // Use different values for each update
            value.set(tio.player() == 0 ? val_to_set : 0);

            // Perform the update
            segTree.Update(tio, mpcio, yield, index, value);
            std::cout << "Update " << (u + 1) << " ends" << std::endl;
        }

        // Print update operation statistics
        std::cout << "===== Updates Stats =====" << std::endl;
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        mpcio.reset_stats();
        tio.reset_lamport();

        #ifdef SEGTREE_VERBOSE
        // Print updated segment tree state for debugging
        std::cout << "\n===== Updated Segment Tree =====" << std::endl;
        segTree.printSegmentTree(tio, yield);
        #endif

        // Perform the specified number of range sum queries
        for (size_t q = 0; q < n_queries; ++q) {
            std::cout << "\n===== Range Sum Query " << (q + 1) << " begins =====" << std::endl;

            RegAS left_index, right_index;

            // Generate valid range query parameters (ensure left <= right)
            size_t max_leaf_idx = (1 << (depth - 1)) - 1;
            size_t left_val = q % (1 << (depth - 1));  // Start anywhere in valid range
            size_t right_val = left_val + (q % (max_leaf_idx + 1 - left_val)); // Ensure right >= left

            left_index.set(tio.player() == 0 ? left_val : 0);
            right_index.set(tio.player() == 0 ? right_val : 0);

            #ifdef SEGTREE_VERBOSE
            auto recons_left = mpc_reconstruct(tio, yield, left_index);
            auto recons_right = mpc_reconstruct(tio, yield, right_index);
            std::cout << "Range Sum Query [" << recons_left << ", " << recons_right << "]" << std::endl;
            #endif

            // Perform the range sum query
            segTree.RangeSum(tio, mpcio, yield, left_index, right_index);
            std::cout << "Range Sum Query " << (q + 1) << " ends" << std::endl;
        }

        // Print final statistics
        std::cout << "===== Range Sum Stats =====" << std::endl;
    });
}
