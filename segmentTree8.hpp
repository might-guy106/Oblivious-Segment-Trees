#include "types.hpp"
#include "mpcio.hpp"
#include "coroutine.hpp"
#include "options.hpp"
#include "mpcops.hpp"
#include "duoram.hpp"
#include <tuple>

class SegmentTree8 {
    private:
    // Single Duoram storing tuples of (value, parent, leftChildSibling, rightChildSibling)
    using SegTreeNode = std::tuple<RegAS, RegAS, RegAS, RegAS>;
    Duoram<SegTreeNode> oram;

    void getBitVector(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, Duoram < RegXS > &bitVec, RegAS left, RegAS right);

    public:

    size_t num_items;
    size_t depth;

    SegmentTree8(int player_num, size_t size, size_t d) : oram(player_num, size) {
            num_items = size;
            depth = d;
            std:: cout << "Segment Tree of depth " << depth << " with " << num_items << " nodes created" << std::endl;
        }

        void init(MPCTIO &tio, yield_t & yield);

    void RangeSum(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, RegAS left, RegAS right);
    // Instrumented Update (global-index version): measures (1) leaf read, (2) per-level update write, (3) parent index read timings/resources
    void Update(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, RegAS index, RegAS value);
    void printSegmentTree(MPCTIO &tio, yield_t & yield);

};

void SegTree8(MPCIO &mpcio, const PRACOptions &opts, char **args);
