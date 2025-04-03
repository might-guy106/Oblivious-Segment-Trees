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

    for(size_t i=1; i<num_items; i++) {
        value_t recons = mpc_reconstruct(tio, yield, nextRArray[i]);
        std::cout << "nextRArray[" << i << "] = " << recons << std::endl;
    }

    auto parentArray = parent.flat(tio, yield);
    parentArray.init([this] (size_t i) -> size_t {
        if (i >= 1 && i < num_items) {
            return i / 2;
        } else {
            return size_t(0);
        }
    });
}

void SegmentTree::getBitVector(MPCTIO &tio, yield_t & yield, Duoram < RegXS > &bitVec, RegAS left, RegAS right) {

    auto bitVecArray = bitVec.flat(tio, yield);
    auto isEvenArray = isEven.flat(tio, yield);
    auto nextLArray = nextL.flat(tio, yield);
    auto nextRArray = nextR.flat(tio, yield);

    uint32_t height = static_cast<uint32_t>(std::log2(num_items));

    RegXS incl;
    incl.set(1);
    RegXS excl;
    excl.set(0);

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

        RegXS leftNode = bitVecLevel[left];
        RegXS rightNode = bitVecLevel[right];
        RegBS one;
        one.set(1);
        RegBS valid;
        mpc_or(tio, yield, valid, eq_c, gt_c);

        value_t recons = mpc_reconstruct(tio, yield, valid);
        std::cout << "Valid = " << recons << std::endl;

        RegXS isEvenL = isEvenLevel[left];
        RegBS inclusionL = isEvenL.bitat(0);
        RegBS f;
        mpc_and(tio, yield, f, valid, (inclusionL ^ one));
        mpc_select(tio, yield, leftNode, f, excl, incl);
        bitVecLevel[left] = leftNode;

        value_t recons1 = mpc_reconstruct(tio, yield, f);
        std::cout << "Left Node = " << recons1 << std::endl;

        RegXS isEvenR = isEvenLevel[right];
        RegBS inclusionR = isEvenR.bitat(0);
        RegBS g;
        mpc_and(tio, yield, g, valid, inclusionR);
        mpc_select(tio, yield, rightNode, g, excl, incl);
        bitVecLevel[right] = rightNode;

        value_t recons2 = mpc_reconstruct(tio, yield, g);
        std::cout << "Right Node = " << recons2 << std::endl;

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

        // value_t answer = mpc_reconstruct(tio, yield, sum);
        // std::cout << "Sum = " << answer << std::endl;

        // value_t recons = mpc_reconstruct(tio, yield, element);
        // std::cout << "Element = " << recons << std::endl;
    }

    value_t answer = mpc_reconstruct(tio, yield, sum);
    std::cout << "Sum = " << answer << std::endl;
}    

void SegmentTree::Update(MPCTIO &tio, yield_t & yield, RegAS index, RegAS value) {
    auto SegTreeArray = oram.flat(tio, yield);
    auto parentArray = parent.flat(tio, yield);
    
    RegAS disp;
    disp.set(16);

    RegAS index1 = index + disp;

    RegAS currVal = SegTreeArray[index1];
    RegAS diff = value - currVal;

    SegTreeArray[index1] = value;
    for(size_t i=1; i<=4; i++) {
        RegAS parentIndex = parentArray[index1];
        SegTreeArray[parentIndex] += diff;
        index1 = parentIndex;
    }

    RegAS checkParent;
    checkParent.set(12);
    RegAS parentValue = SegTreeArray[checkParent];
    value_t recons = mpc_reconstruct(tio, yield, parentValue);
    std::cout << "Parent Value = " << recons << std::endl;
}


void SegTree(MPCIO &mpcio, const PRACOptions &opts, char **args) {
    nbits_t depth = 5;

    // if (*args) {
    //     depth = atoi(*args);
    //     ++args;
    // }
    
    address_t len = (1<<depth);

    MPCTIO tio(mpcio, 0, opts.num_cpu_threads);

    run_coroutines(tio, [&tio, len] (yield_t &yield) {
        SegmentTree segTree(tio.player(), len);
        segTree.init(tio, yield);

        std::cout << "Update begins" << std::endl;
        RegAS index;
        index.set(8);
        RegAS value;
        value.set(100);
        segTree.Update(tio, yield, index, value);
        std::cout << "Update ends" << std::endl;

        std::cout << "Range Sum begins" << std::endl;
        RegAS left_index, right_index;
        left_index.set(0);
        right_index.set(9);
        segTree.RangeSum(tio, yield, left_index, right_index);
        std::cout << "Range Sum ends" << std::endl;
    });
}