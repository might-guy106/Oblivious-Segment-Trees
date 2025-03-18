/**
 * The Segment tree is maintained in a way similar to the heap.
 * We wish to support two kinds of queries: Range-sum and Update single element
 * 
 * Optimized Range Sum Query:
 * We want to compute a bit vector indicating the nodes that are to be included
for covering the given range.
We start from the leaf layer and move upwards gradually. While processing
each layer, we look at the interval (l, r) yet to be covered and do the following:
• If the node corresponding to l, is a right child of its parent, include the
current node and set l = r + 1 (because that interval is now covered).
• If the node corresponding to r, is a left child of its parent, include the
current node and set r = l − 1 (because that interval is now covered).
Move upwards layer by layer as long as l ≤ r or stop if we reached the root
(include the root).

To support this algorithm, we explicitly maintain whether a node is left or right child of its parent.alignas

Optimized Update Query:
It is similar to update operation in a heap. Use a single IDPF.
 */

#include <functional>
#include <vector>
#include <cmath>
#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include "rdpf.hpp"
#include "shapes.hpp"
#include "segment_tree.hpp"
#include "mpcops.hpp"
#include "cdpf.hpp"

void SegmentTree::init(MPCTIO &tio, yield_t &yield, size_t n) {
    // auto TreeArray = oram.flat(tio, yield);

    // TreeArray.init([n](size_t i) {
    //     if (i >= 1 && i <= n) {
    //         return size_t(n/i);
    //     } else {
    //         return size_t(0x7fffffffffffffff);
    //     }
    // }); ???

    auto isLeftArray = isLeft.flat(tio, yield);
    isLeftArray.init([n](size_t i) {
            return !(i%2);
    });

    auto isRightArray = isRight.flat(tio, yield);   
    isRightArray.init([n](size_t i) {
            return (i%2);
    });

    auto nextLArray = nextL.flat(tio, yield);
    nextLArray.init([n](size_t i) {
        return size_t((i+1)/2);
    });

    auto nextRArray = nextR.flat(tio, yield);
    nextRArray.init([n](size_t i) {
        return size_t((i-1)/2);
    });

}

Duoram<RegBS> SegmentTree::computeBitVector(MPCTIO &tio, yield_t & yield, RegAS left_index, RegAS right_index) {
    auto TreeArray = oram.flat(tio, yield);
    auto isLeftArray = isLeft.flat(tio, yield);
    auto isRightArray = isRight.flat(tio, yield);
    auto nextLArray = nextL.flat(tio, yield);
    auto nextRArray = nextR.flat(tio, yield);

    size_t levels = size_t(log2(num_items)) + 1;

    Duoram<RegBS> bitVec(tio.player(), num_items + 1);
    auto bitVector = bitVec.flat(tio, yield);

    for(size_t i = levels - 1; i >= 0; i--) {
        
        CDPF cdpf = tio.cdpf(yield);
        auto[lt_c, eq_c, gt_c] = cdpf.compare(tio, yield, left_index- right_index, tio.aes_ops());

        RegBS setNode;
        mpc_or(tio, yield, &setNode, eq_c, bitVector[left_index]);
        bitVector[left_index] = setNode;

        RegBS includeLeft, includeRight;

        run_coroutines(tio, [&tio, &isRightArray, &lt_c, &left_index, &includeLeft](yield_t &yield) {
            auto iRAcoro = isRightArray.context(yield);
            mpc_and(tio, yield, includeLeft, iRAcoro[left_index], lt_c);},
            [&tio, &isLeftArray, &lt_c, &right_index, &includeRight](yield_t &yield) {
            auto iLAcoro = isLeftArray.context(yield);  
            mpc_and(tio, yield, includeRight, iLAcoro[right_index], lt_c);}); 

        run_coroutines(tio, [&tio, &includeLeft, &left_index, &bitVector, &setNode](yield_t &yield) {
            auto bVecContext = bitVector.context(yield);
            mpc_or(tio, yield, &setNode, includeLeft, bVecContext[left_index]);
            bVecContext[left_index] = setNode;},
            [&tio, &includeRight, &right_index, &bitVector, &setNode](yield_t &yield) {
            auto bVecContext = bitVector.context(yield);
            mpc_or(tio, yield, &setNode, includeRight, bVecContext[right_index]);
            bVecContext[right_index] = setNode;});
        
        run_coroutines(tio, [&tio, &nextLArray, &left_index, &left_index](yield_t &yield) {
            auto nLAcoro = nextLArray.context(yield);
            left_index = nLAcoro[left_index];},
            [&tio, &nextRArray, &right_index, &right_index](yield_t &yield) {
            auto nRAcoro = nextRArray.context(yield);
            right_index = nRAcoro[right_index];});
    }

    return bitVec;
}

void SegTree(MPCIO &mpcio, const PRACOptions &opts, char **args) {
    nbits_t depth = 5;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }
    
    address_t len = (1<<depth);

    MPCTIO tio(mpcio, 0, opts.num_cpu_threads);
    run_coroutines(tio, [&tio, len] (yield_t &yield) {
        SegmentTree segTree(tio.player(), len - 1);
        segTree.init(tio, yield, len);
        RegAS left_index, right_index;
        left_index.set(len/2);
        right_index.set(len/2 + 4);
        Duoram<RegBS> bitVec = segTree.computeBitVector(tio, yield, left_index, right_index);

        std::cout<<"Here's the required bit vector:"<< std::endl;

        auto bitArray = bitVec.flat(tio, yield);
        for(size_t i = 1; i <= len; i++) {
            std::cout<<i<<" " << mpc_reconstruct(tio, yield, bitArray[i]) << std::endl;
        }
    });
}