#include "types.hpp"
#include "mpcio.hpp"
#include "coroutine.hpp"
#include "options.hpp"
#include "mpcops.hpp"
#include "duoram.hpp"


class SegmentTree12 {
    private:
    Duoram<RegAS> TreeOram;    // Duoram to store complete tree
    Duoram<RegAS> SiblingOram; // Duoram to store sibling pointers (absolute indices)

    public:
    size_t num_items;
    size_t depth;

    SegmentTree12(int player_num, size_t size, size_t d) : TreeOram(player_num, size), SiblingOram(player_num, size) {
        num_items = size;
        depth = d;
        std::cout << "Segment Tree of depth " << depth << " with " << num_items << " nodes created" << std::endl;
    }

    void init(MPCTIO &tio, yield_t &yield);

    // Range sum over absolute indices in the tree.
    // NOTE: left and right are absolute tree indices (NOT level-relative leaf indices).
    RegAS RangeSum(MPCTIO &tio, MPCIO &mpcio, yield_t &yield, RegXS left, RegXS right);

    // Update at an absolute index in the tree.
    void Update(MPCTIO &tio, MPCIO &mpcio, yield_t &yield, RegXS index, RegAS value);

    void printSegmentTree(MPCTIO &tio, yield_t &yield);
};

void SegTree12(MPCIO &mpcio, const PRACOptions &opts, char **args);
