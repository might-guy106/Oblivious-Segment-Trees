#include <functional>
#include <chrono>
#include <tuple>
#include "mpcops.hpp"
#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include "rdpf.hpp"
#include "shapes.hpp"
#include "segmentTree8.hpp"

// uncomment to enable intermediate reconstructions and logging
#define SEGTREE_VERBOSE

/*
The segment tree data structure in SegmentTree8 uses a LEVEL-WISE RESTRUCTURED ARRAY layout 

ARRAY STRUCTURE:
Unlike the traditional binary tree layout, this version reorganizes the data structure where:
- Each level k has 2^k + 1 nodes (including one extra node at position 0)
- Total array size: 2^d + d - 1 nodes
- Level k starts at index: 2^k + k - 1

LEVEL-WISE LAYOUT EXAMPLE:
For depth d=3, array [100, 200, 300, 400] (4 elements):

Level 0: [0][600]           (Root level: extra_node + root_sum)
Level 1: [0][300][700]      (Level 1: extra_node + left_sum + right_sum) 
Level 2: [0][100][200][300][400]  (Leaf level: extra_node + leaf_values)

Array storage:
Index:  [0][1][2][3][4][5][6][7][8]
Value:  [0][600][0][300][700][0][100][200][300][400]
Level:   0  0   1   1    1   2   2    2    2    2

Where:
- Position 0 in each level: Extra node (unused, set to 0)
- Position 1+ in each level: Actual segment tree nodes

MATHEMATICAL FORMULAS:
- Level k start index: getLevelStart8(k) = 2^k + k - 1
- Level k length: getLevelLength8(k) = 2^k + 1  
- Total array size: 2^d + d - 1

ADVANTAGES OF LEVEL-WISE RE-STRUCTURE:
1. Efficient level-by-level processing using Flat objects
2. Reduce number of Duoram access (without this restructuring we also need to access isEven[left] and isEven[right] to check whether left is left child or right is right child)
*/
    
/*
RANGE SUM QUERY OPERATION (Level-wise Processing)

The range sum query efficiently computes the sum of elements in a given range [left, right]
using the level-wise segment tree structure. The algorithm processes entire levels at once
using level-specific Flat objects for optimal MPC performance.

Algorithm Overview:
1. Start with leaf-level indices (left, right) 
2. For each level from leaves to root:
   - Create level-specific Flat objects for bit vector, siblings, and parents
   - Read parent and sibling indices using level-relative addressing
   - Check convergence (parents equal) and validity (left <= right)
   - Mark appropriate sibling nodes in the bit vector (and also mark the endpoints themselves if it is leaf level)
   - Move to parent level using parent indices
3. Accumulate sums level-wise, skipping extra nodes (position 0 in each level)
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


// // Helper function which returns the level given the index in the global Segment Tree array  
// size_t findLevel8(size_t index, size_t depth) {
//     for (size_t level = 0; level < depth; level++) {
//         size_t levelStart = (1ULL << level) + level - 1;
//         size_t levelEnd = (1ULL << (level + 1)) + level - 1;
//         if (index >= levelStart && index <= levelEnd) {
//             return level;
//         }
//     }
//     return SIZE_MAX; // Invalid
// }

// // Helper functions which returns the start index in the global Segment Tree array for a given level
// size_t getLevelStart8(size_t level) {
//     return (1ULL << level) + level - 1;  // 2^level + level - 1
// }

// // Helper function which returns the length of a given level.
// size_t getLevelLength8(size_t level) {
//     return (1ULL << level) + 1;  // 2^level + 1
// }

// // Given an index in the global Segment Tree array, return its parent index in its level 
// size_t parentIndex8(size_t index, size_t depth) {
//     size_t currentLevel = findLevel8(index, depth);
    
//     if (currentLevel <= 0) {
//         return index;  // Root level
//     }
    
//     size_t levelStart = getLevelStart8(currentLevel);
    
//     // If it's the extra node, return itself
//     if (index == levelStart) {
//         return 0;
//     }
    
//     // For regular nodes
//     size_t positionInLevel = index - levelStart;
//     size_t parentPosition = (positionInLevel - 1) / 2 + 1;
    
//     return parentPosition;
// }

// // Given an index in the global Segment Tree array, return its left child sibling index in its level (if index is not a left child then it returns 0)
// size_t getLeftChildSiblingIndex8(size_t index, size_t depth) {
//     size_t currentLevel = findLevel8(index, depth);
    
//     if (currentLevel <= 0) {
//         return index;  // Root level, return itself
//     }
    
//     size_t levelStart = getLevelStart8(currentLevel);
    
//     // If it's the extra node, return itself
//     if (index == levelStart) {
//         return 0;
//     }
    
//     // For regular nodes
//     size_t positionInLevel = index - levelStart;
    
//     if(positionInLevel % 2 == 1) { // odd index implies left child
//         return positionInLevel + 1; // return sibling
//     } else {
//         return 0;
//     }
// }

// // Given an index in the global Segment Tree array, return its right child sibling index in its level (if index is not a right child then it returns 0)
// size_t getRightChildSiblingIndex8(size_t index, size_t depth) {
//     size_t currentLevel = findLevel8(index, depth);

//     if (currentLevel <= 0) {
//         return index;  // Root level, return itself
//     }
    
//     size_t levelStart = getLevelStart8(currentLevel);
    
//     // If it's the extra node, return itself
//     if (index == levelStart) {
//         return 0;
//     }
    
//     // For regular nodes
//     size_t positionInLevel = index - levelStart;

//     if(positionInLevel % 2 == 0) { // even index implies right child
//         return positionInLevel - 1; // return sibling
//     } else {
//         return 0;
//     }
// }

// void SegmentTree8::init(MPCTIO &tio, yield_t & yield) {
    
//     // Create and initialize  a helper array which is used to initialize the Segment Tree Duoram
//     std::vector<uint64_t> segTree(num_items, 0);
    
//     // Initialize the normal binary tree values first
//     size_t normalTree_num_leaves = 1ULL << (depth - 1);
//     size_t normalTree_size = normalTree_num_leaves * 2;
//     std::vector<uint64_t> normalTree(normalTree_size, 0);

//     // Initialize leaf nodes with values 0, 100, 200, ..., (num_leaves-1)*100
//     for (size_t i = 0; i < normalTree_num_leaves; ++i) {
//         normalTree[normalTree_num_leaves + i] = i * 100;
//     }
    
//     // Build internal nodes bottom-up by computing range sums
//     for (int i = normalTree_num_leaves - 1; i >= 1; --i) {
//         normalTree[i] = normalTree[2*i] + normalTree[2*i + 1];
//     }
    
//     // Map original tree to new level-wise structure
//     for (size_t level = 0; level < depth; level++) {
//         size_t levelStart = getLevelStart8(level);
//         size_t levelLength = getLevelLength8(level);
        
//         // Set extra node (position 0) to 0
//         segTree[levelStart] = 0;
        
//         // Map regular nodes
//         for (size_t pos = 1; pos < levelLength; pos++) {
//             size_t newIndex = levelStart + pos;
//                 size_t normalTreeIndex = (1ULL << level) + pos - 1;
//                 segTree[newIndex] = normalTree[normalTreeIndex];
//         }
//     }
    
//     // Initialize Duoram with tuples (value, parent, leftChildSibling, rightChildSibling)
//     auto SegTreeArray = oram.flat(tio, yield);
//     SegTreeArray.init([this, &segTree] (size_t i) -> std::tuple<RegAS, RegAS, RegAS, RegAS> {
//         RegAS value, parent, leftSibling, rightSibling;
        
//         if (i < num_items) {
//             value.set(segTree[i]);
//             parent.set(parentIndex8(i, depth));
//             leftSibling.set(getLeftChildSiblingIndex8(i, depth));
//             rightSibling.set(getRightChildSiblingIndex8(i, depth));
//         } else {
//             value.set(0x7fffffffffffffff);
//             parent.set(0);
//             leftSibling.set(0);
//             rightSibling.set(0);
//         }
        
//         return std::make_tuple(value, parent, leftSibling, rightSibling);
//     });
// }

// // helper function to reconstruct and print SegTreeArray
// // helper function to reconstruct and print SegTreeArray
// void SegmentTree8::printSegmentTree(MPCTIO &tio, yield_t & yield) {
//     auto SegTreeArray = oram.flat(tio, yield);
//     auto SegTreeRecons = SegTreeArray.reconstruct();
//     for(size_t i=0; i<num_items; i++) {
//         std::cout << "SegTreeArray[" << i << "] = (value: " 
//                   << std::get<0>(SegTreeRecons[i]).share() << ", parent: "
//                   << std::get<1>(SegTreeRecons[i]).share() << ", leftSib: "
//                   << std::get<2>(SegTreeRecons[i]).share() << ", rightSib: "
//                   << std::get<3>(SegTreeRecons[i]).share() << ")" << std::endl;
//     }
// }

// // Main helper function to get the bit vector for RangeSum
// void SegmentTree8::getBitVector(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, Duoram < RegXS > &bitVec, RegAS leftLevelIndex, RegAS rightLevelIndex) {

//     auto bitVecArray = bitVec.flat(tio, yield);
//     auto SegTreeArray = oram.flat(tio, yield);

//     RegBS isValid; // isValid = right >= left
//     CDPF cdpf = tio.cdpf(yield);
//     RegAS diff = rightLevelIndex - leftLevelIndex;
//     auto[lt_c1, eq_c1, gt_c1] = cdpf.compare(tio, yield, diff, tio.aes_ops());
//     mpc_or(tio, yield, isValid, gt_c1, eq_c1); // isValid = (left <= right)


//     // Step 1: Pre-compute the path indices for all levels
//     auto start_time = std::chrono::high_resolution_clock::now();
//     // Store leftLevelIndex and rightLevelIndex for each level
//     std::vector<SegTreeNode> leftPath(depth);
//     std::vector<SegTreeNode> rightPath(depth);
//     std::vector<RegAS> leftPathIndex(depth);
//     std::vector<RegAS> rightPathIndex(depth);
    
//     // Compute parent indices for all levels going up the tree
//     for(uint32_t i = 1; i <= depth; i++) {
//         size_t level = depth - i;  // Current level (leaf to root)
//         size_t levelStart = getLevelStart8(level);
//         size_t levelLength = getLevelLength8(level);
        
//         typename Duoram<SegTreeNode>::Flat segTreeLevel(SegTreeArray, tio, yield, levelStart, levelLength);
        
//         // Read parent indices (field 1) for next level
//         SegTreeNode leftTupple = segTreeLevel[leftLevelIndex];
//         SegTreeNode rightTupple = segTreeLevel[rightLevelIndex];
//         leftPath[level] = leftTupple;
//         rightPath[level] = rightTupple;
//         leftPathIndex[level] = leftLevelIndex;
//         rightPathIndex[level] = rightLevelIndex;

//         leftLevelIndex = std::get<1>(leftTupple);
//         rightLevelIndex = std::get<1>(rightTupple);
//     }

//     auto step1_end_time = std::chrono::high_resolution_clock::now();

//     // Step 2: Now parallelize the reads and computations for all levels
//     std::vector<coro_t> read_coroutines;
//     for(uint32_t i = 1; i <= depth; i++) {
//         size_t level = depth - i;

//         read_coroutines.emplace_back([&tio, level, i, this, &isValid, &leftPath, &rightPath, &leftPathIndex, &rightPathIndex, &bitVecArray, 
//                                     &SegTreeArray](yield_t &sub_yield) {
//             size_t levelStart = getLevelStart8(level);
//             size_t levelLength = getLevelLength8(level);
            
//             // Create incl fresh in each coroutine to avoid sharing issues
//             RegXS incl;
//             incl.set(tio.player()==0 ? 1 : 0);
//             RegAS one;
//             one.set(tio.player()==0 ? 1 : 0);
            
//             typename Duoram<RegXS>::Flat bitVecLevel(bitVecArray, tio, sub_yield, levelStart, levelLength);
            
//             // Read sibling values using field access
//             // Field 2 is leftChildSibling, Field 3 is rightChildSibling
//             RegAS leftIndex = leftPathIndex[level];
//             RegAS rightIndex = rightPathIndex[level];
//             RegAS leftSibling = std::get<2>(leftPath[level]);
//             RegAS rightSibling = std::get<3>(rightPath[level]);

//             // --- if left and right are adjacent or same then the marking of bitvector is already completed so we set isDone here ---
//             RegBS isNotDone;
//             CDPF cdpf1 = tio.cdpf(sub_yield);
//             RegAS diff1 = rightIndex - (leftIndex + one); // diff1 = right - (left + 1)
//             auto[lt_c1, eq_c1, gt_c1] = cdpf1.compare(tio, sub_yield, diff1, tio.aes_ops());
//             isNotDone = gt_c1; // isNotDone = (right > left + 1)

//             // --- initialize leftSiblingIncluded and rightSiblingIncluded (contains 0 initially) ---
//             RegXS leftSiblingIncluded, rightSiblingIncluded; // will become 'incl' if Check passes later

//             // --- we have to mark the leftSibling and rightSibling only if marking is not yet completed and initial query is valid ---
//             RegBS Check;
//             mpc_and(tio, sub_yield, Check, isNotDone, isValid); // Check = isNotDone AND isValid
//             mpc_select(tio, sub_yield, leftSiblingIncluded, Check, leftSiblingIncluded, incl);
//             mpc_select(tio, sub_yield, rightSiblingIncluded, Check, rightSiblingIncluded, incl);


//             // --- Setting the marks in bitVec using Levelwise array ---
//             bitVecLevel[leftSibling] = leftSiblingIncluded;
//             bitVecLevel[rightSibling] = rightSiblingIncluded;

//             // --- if its the leaf level then we also have to mark the leftLevelIndex and rightLevelIndex (this is there in the normal algorithm so it does not reveal any information) ---
//             if(i == 1)
//             {   
//                 bitVecLevel[leftIndex] = incl;
//                 bitVecLevel[rightIndex] = incl;
//             }
//         });
//     }
    
//     run_coroutines(tio, read_coroutines);

//     auto step2_end_time = std::chrono::high_resolution_clock::now();
//     auto step1_duration = std::chrono::duration_cast<std::chrono::milliseconds>(step1_end_time - start_time).count();
//     auto step2_duration = std::chrono::duration_cast<std::chrono::milliseconds>(step2_end_time - step1_end_time).count();

//     std::cout << "Step 1 (Path Computation) Time: " << step1_duration << " ms" << std::endl;
//     std::cout << "Step 2 (Level-wise Marking) Time: " << step2_duration << " ms" << std::endl;
// }

// // Main RangeSum function
// void SegmentTree8::RangeSum(MPCTIO &tio,  MPCIO &mpcio, yield_t & yield, RegAS left, RegAS right) {
//     Duoram < RegXS > bitVec(tio.player(), num_items);
//     getBitVector(tio, mpcio, yield, bitVec, left, right);

//     auto accumulation_start = std::chrono::high_resolution_clock::now();
    
//     auto bitVecArray = bitVec.flat(tio, yield);
//     auto SegTreeArray = oram.flat(tio, yield);

//     // Store partial sums for each level
//     std::vector<RegAS> levelSums(depth);
//     for(size_t level = 0; level < depth; level++) {
//         levelSums[level].set(0);
//     }

//     // Parallelize accumulation across all levels
//     std::vector<coro_t> accumulation_coroutines;
//     for(size_t level = 0; level < depth; level++) {
//         accumulation_coroutines.emplace_back([&tio, level, this, &levelSums, &bitVecArray, &SegTreeArray](yield_t &sub_yield) {
//             size_t levelStart = getLevelStart8(level);
//             size_t levelLength = getLevelLength8(level);

//             typename Duoram<RegXS>::Flat bitVecLevel(bitVecArray, tio, sub_yield, levelStart, levelLength);
//             typename Duoram<SegTreeNode>::Flat segTreeLevel(SegTreeArray, tio, sub_yield, levelStart, levelLength);
            
//             // Iterate from 1 to avoid the extra node (position 0 in each level)
//             for(size_t pos = 1; pos < levelLength; pos++) {
//                 RegXS element = bitVecLevel[pos];
//                 RegBS incl = element.bitat(0);
//                 // Access value field (field 0) of the tuple
//                 SegTreeNode curTupple = segTreeLevel[pos];
//                 RegAS val = std::get<0>(curTupple);
//                 RegAS zero;
//                 zero.set(0);

//                 RegAS sum1;
//                 mpc_select(tio, sub_yield, sum1, incl, zero, val);

//                 levelSums[level].ashare += sum1.ashare;
//             }
//         });
//     }
    
//     run_coroutines(tio, accumulation_coroutines);

//     // Combine all level sums
//     RegAS sum;
//     sum.set(0);
//     for(size_t level = 0; level < depth; level++) {
//         sum.ashare += levelSums[level].ashare;
//     }

//     auto accumulation_end = std::chrono::high_resolution_clock::now();
//     auto accumulation_duration = std::chrono::duration_cast<std::chrono::milliseconds>(accumulation_end - accumulation_start).count();
//     std::cout << "Step 3 (Sum Accumulation - Parallelized) Time: " << accumulation_duration << " ms" << std::endl;

//     #ifdef SEGTREE_VERBOSE
//     value_t answer = mpc_reconstruct(tio, yield, sum);
//     std::cout << "Sum = " << answer << std::endl;
//     #endif
// }

// // Main Update function
// void SegmentTree8::Update(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, RegAS index, RegAS value) {
//     auto update_start = std::chrono::high_resolution_clock::now();
    
//     auto SegTreeArray = oram.flat(tio, yield);

//     size_t leafLevel = depth - 1;
//     size_t leafStart = getLevelStart8(leafLevel);

//     // Access last level Flat object and get current value to calculate difference
//     size_t leafLength = getLevelLength8(leafLevel);
//     typename Duoram<SegTreeNode>::Flat leafLevel_flat(SegTreeArray, tio, yield, leafStart, leafLength);
//     // Access value field (field 0)
//     SegTreeNode leafNode = leafLevel_flat[index];
//     RegAS currVal = std::get<0>(leafNode);
//     RegAS diff = value - currVal;

//     // Debugging and intermediate reconstructions
//     #ifdef SEGTREE_VERBOSE
//     auto recons_index = mpc_reconstruct(tio, yield, index);
//     auto recons_currVal = mpc_reconstruct(tio, yield, currVal);
//     auto recons_newvalue = mpc_reconstruct(tio, yield, value);
//     auto recons_diff = mpc_reconstruct(tio, yield, diff);
//     std::cout << "Index to be updated in the original array = " << recons_index << std::endl;
//     std::cout << "Index to be updated in the segment tree array = " << leafStart + recons_index << std::endl;
//     std::cout << "Current Value at index = " << recons_currVal << std::endl;
//     std::cout << "New Value to be updated = " << recons_newvalue << std::endl;
//     std::cout << "Diff = " << recons_diff << std::endl;
//     #endif
    
//     // Phase 1: Pre-compute the path indices for all levels (similar to getBitVector)
//     std::vector<SegTreeNode> updatePath(depth);
//     std::vector<RegAS> updatePathIndex(depth);
    
//     // Compute parent indices for all levels going up the tree
//     for(size_t i = 1; i <= depth; i++) {
//         size_t level = depth - i;  // Current level (leaf to root)
//         size_t levelStart = getLevelStart8(level);
//         size_t levelLength = getLevelLength8(level);
        
//         typename Duoram<SegTreeNode>::Flat segTreeLevel(SegTreeArray, tio, yield, levelStart, levelLength);
        
//         // Read parent index (field 1) for next level
//         SegTreeNode currentNode = segTreeLevel[updatePathIndex[level]];
//         updatePath[level] = currentNode;
//         updatePathIndex[level] = index;

//         index = std::get<1>(currentNode);
//     }
    
//     auto phase1_end = std::chrono::high_resolution_clock::now();
    
//     // Phase 2: Parallelize the updates across all levels
//     std::vector<coro_t> update_coroutines;
//     for(size_t i = 1; i <= depth; i++) {
//         size_t level = depth - i;

//         update_coroutines.emplace_back([&tio, level, this, &diff, &updatePathIndex, &updatePath, &SegTreeArray](yield_t &sub_yield) {
//             size_t levelStart = getLevelStart8(level);
//             size_t levelLength = getLevelLength8(level);
            
//             typename Duoram<SegTreeNode>::Flat segTreeLevel(SegTreeArray, tio, sub_yield, levelStart, levelLength);
            
//             // Update the value field (field 0) at this level's index
//             std::get<0>(updatePath[level]) += diff;
//             segTreeLevel[updatePathIndex[level]] = updatePath[level];

//             #ifdef SEGTREE_VERBOSE
//             auto recons_Index = mpc_reconstruct(tio, sub_yield, updatePathIndex[level]);
//             SegTreeNode updated_val = segTreeLevel[updatePathIndex[level]];
//             auto recons_updated = mpc_reconstruct(tio, sub_yield, std::get<0>(updated_val));
//             std::cout << "Updated Index = " << (levelStart + recons_Index) << " with value = " << recons_updated << std::endl;
//             #endif
//         });
//     }
    
//     run_coroutines(tio, update_coroutines);
    
//     auto phase2_end = std::chrono::high_resolution_clock::now();
//     auto phase1_duration = std::chrono::duration_cast<std::chrono::milliseconds>(phase1_end - update_start).count();
//     auto phase2_duration = std::chrono::duration_cast<std::chrono::milliseconds>(phase2_end - phase1_end).count();
    
//     std::cout << "Update Phase 1 (Path Computation) Time: " << phase1_duration << " ms" << std::endl;
//     std::cout << "Update Phase 2 (Parallel Updates) Time: " << phase2_duration << " ms" << std::endl;
// }

// Main function to run Segment Tree 6 operations
void SegTree8(MPCIO &mpcio, const PRACOptions &opts, char **args) {
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

    // run_coroutines(tio, [&tio, &mpcio, len, depth, n_updates, n_queries] (yield_t &yield) {
        
    //     // Time initialization
    //     auto init_start = std::chrono::high_resolution_clock::now();
        
    //     SegmentTree8 segTree(tio.player(), len, depth);
    //     segTree.init(tio, yield);
        
    //     auto init_end = std::chrono::high_resolution_clock::now();
    //     auto init_duration = std::chrono::duration_cast<std::chrono::milliseconds>(init_end - init_start);
        
    //     std::cout << "===== Segment Tree Init Stats =====" << std::endl;
    //     std::cout << "Updates: " << n_updates << ", Queries: " << n_queries << std::endl;
    //     std::cout << "Initialization Time: " << init_duration.count() / 1000.0 << " seconds" << std::endl;
    //     tio.sync_lamport();
    //     mpcio.dump_stats(std::cout);
    //     mpcio.reset_stats();
    //     tio.reset_lamport();

    //     // Print initial segment tree
    //     #ifdef SEGTREE_VERBOSE
    //     std::cout << "\n===== Initial Segment Tree =====" << std::endl;
    //     segTree.printSegmentTree(tio, yield);
    //     #endif

    //     // Time updates
    //     auto update_start = std::chrono::high_resolution_clock::now();
        
    //     // Perform updates
    //     for (size_t u = 0; u < n_updates; ++u) {
    //         std::cout << "\n===== Update " << (u + 1) << " begins =====" << std::endl;
    //         RegAS index;
    //         size_t idx_val = (u % (1 << (depth - 1))) + 1; // Cycle through leaf indices (1 to 2^(d-1))
    //         index.set(tio.player() == 0 ? idx_val : 0); // Cycle through leaf indices


    //         RegAS value;
    //         size_t val_to_set = (u + 1) * 50;
    //         value.set(tio.player() == 0 ? val_to_set : 0); // Use different values for each update

    //         segTree.Update(tio, mpcio, yield, index, value);
    //         std::cout << "Update " << (u + 1) << " ends" << std::endl;
    //     }

    //     auto update_end = std::chrono::high_resolution_clock::now();
    //     auto update_duration = std::chrono::duration_cast<std::chrono::milliseconds>(update_end - update_start);

    //     std::cout << "===== Updates Stats =====" << std::endl;
    //     std::cout << "Total Updates Time: " << update_duration.count() / 1000.0 << " seconds" << std::endl;
    //     tio.sync_lamport();
    //     mpcio.dump_stats(std::cout);
    //     mpcio.reset_stats();
    //     tio.reset_lamport();

    //     // Print updated segment tree
    //     #ifdef SEGTREE_VERBOSE
    //     std::cout << "\n===== Updated Segment Tree =====" << std::endl;
    //     segTree.printSegmentTree(tio, yield);
    //     #endif

    //     // Time range queries
    //     auto query_start = std::chrono::high_resolution_clock::now();
        
    //     // Perform range sum queries
    //     for (size_t q = 0; q < n_queries; ++q) {    
    //         std::cout << "\n===== Range Sum Query " << (q + 1) << " begins =====" << std::endl;

    //         RegAS left_index, right_index;

    //         // Initializing valid left and right indices for range sum query
    //         size_t max_leaf_idx = (1 << (depth - 1));
    //         size_t left_val = (q % (1 << (depth - 1))) + 1;  // Start anywhere in valid range (1 to 2^(d-1))
    //         size_t right_val = left_val + (q % (max_leaf_idx + 1 - left_val)); // Ensuring right >= left (not needed as our getBitVector function handles handles invalid cases too)

    //         left_index.set(tio.player() == 0 ? left_val : 0); // Start of range
    //         right_index.set(tio.player() == 0 ? right_val : 0); // End of range

    //         #ifdef SEGTREE_VERBOSE
    //         auto recons_left = mpc_reconstruct(tio, yield, left_index);
    //         auto recons_right = mpc_reconstruct(tio, yield, right_index);
    //         std::cout << "Range Sum Query [" << recons_left << ", " << recons_right << "]" << std::endl;
    //         #endif

    //         segTree.RangeSum(tio, mpcio, yield, left_index, right_index);
    //         std::cout << "Range Sum Query " << (q + 1) << " ends" << std::endl;
    //     }

    //     auto query_end = std::chrono::high_resolution_clock::now();
    //     auto query_duration = std::chrono::duration_cast<std::chrono::milliseconds>(query_end - query_start);

    //     std::cout << "===== Range Sum Stats =====" << std::endl;
    //     std::cout << "Total Range Queries Time: " << query_duration.count() / 1000.0 << " seconds" << std::endl;
    // });
}