#include "coroutine.hpp"
#include "duoram.hpp"
#include "mpcio.hpp"
#include "mpcops.hpp"
#include "options.hpp"
#include "types.hpp"

class SegmentTreeOpt {
  private:
    Duoram<RegAS> TreeOram; // Duoram to store complete tree

  public:
    size_t num_items;
    size_t depth;

    SegmentTreeOpt(int player_num, size_t size, size_t d) : TreeOram(player_num, size) {
        num_items = size;
        depth = d;
        std::cout << "Segment Tree of depth " << depth << " with " << num_items << " nodes created" << std::endl;
    }

    void init(MPCTIO &tio, yield_t &yield);

    // Range sum over absolute indices in the tree.
    // NOTE: left and right are absolute tree indices (NOT level-relative leaf indices).
    RegAS RangeSum(MPCTIO &tio, MPCIO &mpcio, yield_t &yield, RegXS left, RegXS right);
    void Update(MPCTIO &tio, MPCIO &mpcio, yield_t &yield, RegXS index, RegAS value);
    void printSegmentTree(MPCTIO &tio, yield_t &yield);
};

void SegTreeOpt(MPCIO &mpcio, const PRACOptions &opts, char **args);
