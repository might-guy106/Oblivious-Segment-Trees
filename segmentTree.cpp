#include <functional>
#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include "rdpf.hpp"
#include "shapes.hpp"
#include "segmentTree.hpp"

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

    uint64_t segTree[32] = {
        0,    // index 0 unused
        136,    // index 1:  sum of all leaves = 1+2+...+16
        36,  100,   // indices 2 and 3
        10,  26,  42,  58,   // indices 4,5,6,7
        3,   7,  11,  15,  19,  23,  27,  31,  // indices 8 to 15
        1,   2,   3,   4,   5,   6,   7,   8,   // indices 16 to 23 (leaves)
        9,  10,  11,  12,  13,  14,  15,  16   // indices 24 to 31 (leaves)
   };

   SegTreeArray.init((size_t i) {
    if (i >= 1 && i < num_items) {
        return segTree[i];
    } else {
        return size_t(0x7fffffffffffffff);
    }
   });

   auto isEvenArray = isEven.flat(tio, yield);
   isEvenArray.init((size_t i) {
     if (i >= 1 && i < num_items) {
          return (i % 2 == 0) ? 1 : 0;
     } else {
          return 1; // 0 is even, just for consistency
     }
    });

    auto nextLArray = nextL.flat(tio, yield);
    nextLArray.init((size_t i) {
        if (i >= 1 && i < num_items) {
            return ((arrayIndexToLevelPos(i) + 1)/2);
        } else {
            return 0;
        }
    });

    auto nextRArray = nextR.flat(tio, yield);
    nextRArray.init((size_t i) {
        if (i >= 1 && i < num_items) {
            return ((arrayIndexToLevelPos(i) - 1)/2);
        } else {
            return 0;
        }
    });
}

void SegmentTree::getBitVector(MPCTIO &tio, yield_t & yield, Duoram < RegXS > &bitVec, RegAS left, RegAS right) {

    auto bitVecArray = bitVec.flat(tio, yield);
    auto isEvenArray = isEven.flat(tio, yield);
    auto nextLArray = nextL.flat(tio, yield);
    auto nextRArray = nextR.flat(tio, yield);

    uint32_t height = static_cast<uint32_t>(std::log2(num_items));

    RegXS incl(1);
    RegXS excl(0);

    for(uint32_t i=1; i<=height; i++)
    {
        CDPF cdpf = tio.cdpf(yield);
        RegAS diff = right - left;
        auto[lt_c, eq_c, gt_c] = cdpf.compare(tio, yield, diff, tio.aes_ops());

        RegXS leftNode = bitVecArray[left];
        RegXS rightNode = bitVecArray[right];
        RegBS one(1);

        RegBS inclusionL = isEvenArray[left].bit(0);
        RegBS f = mpc_and(mpc_or(eq_c, gt_c), (inclusionL ^ one));
        mpc_select(tio, yield, &leftNode, f, excl, incl);
        bitVecArray[left] = leftNode;

        RegBS inclusionR = isEvenArray[right].bit(0);
        RegBS g = mpc_and(mpc_or(eq_c, gt_c), inclusionR);
        mpc_select(tio, yield, &rightNode, g, excl, incl);
        bitVecArray[right] = rightNode;

        left = nextLArray[left];
        right = nextRArray[right];
    }
}
    
void SegmentTree::RangeSum(MPCTIO &tio, yield_t & yield, RegAS left, RegAS right) {
    Duoram < RegXS > bitVec(tio.player(), num_items);
    getBitVector(tio, yield, bitVec, left, right);

    auto bitVecArray = bitVec.flat(tio, yield);
    auto SegTreeArray = oram.flat(tio, yield);

    RegAS sum(0);

    for(size_t i=1; i<num_items; i++) {
        RegBS incl = bitVecArray[i].bit(0);
        RegAS val = SegTreeArray[i];
        RegAS zero(0);

        RegAS sum1;
        mpc_select(tio, yield, &sum1, incl, zero, val);

        sum += sum1;
    }

    value_t answer = mpc_reconstruct(tio, yield, sum);
    std::cout << "Sum = " << answer << std::endl;
}    

void SegTree(MPCIO &mpcio, const PRACOptions &opts, char **args) {
    nbits_t depth = 5;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }
    
    address_t len = (1<<depth) + 1;

    MPCTIO tio(mpcio, 0, opts.num_cpu_threads);


    run_coroutines(tio, [&tio, len] (yield_t &yield) {
        SegmentTree segTree(tio.player(), len - 1);
        segTree.init(tio, yield, len);
        RegAS left_index, right_index;
        left_index.set(len/2);
        right_index.set(len/2 + 4);
        segTree.RangeSum(tio, yield, left_index, right_index);
    });
}