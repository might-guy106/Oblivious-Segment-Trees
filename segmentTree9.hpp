#include "types.hpp"
#include "mpcio.hpp"
#include "coroutine.hpp"
#include "options.hpp"
#include "mpcops.hpp"
#include "duoram.hpp"

// Forward declaration
class PerformanceLogger;


class SegmentTree9 {
    private:

    Duoram<RegAS> TreeOram; // Duoram to store complete tree
    Duoram<RegAS> SiblingOram; // Duoram to store sibling indexes

    RegAS computeRangeSum(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, RegAS left, RegAS right,
                          PerformanceLogger* logger = nullptr, size_t operation_id = 0);

    public:

    size_t num_items;
    size_t depth;

    SegmentTree9(int player_num, size_t size, size_t d) : TreeOram(player_num, size) , SiblingOram(player_num, size) {
        num_items = size;
        depth = d;
        std:: cout << "Segment Tree of depth " << depth << " with " << num_items << " nodes created" << std::endl;
    }

    void init(MPCTIO &tio, yield_t & yield);

    void RangeSum(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, RegAS left, RegAS right,
                  PerformanceLogger* logger = nullptr, size_t operation_id = 0);
    void Update(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, RegAS index, RegAS value,
                PerformanceLogger* logger = nullptr, size_t operation_id = 0);
    void printSegmentTree(MPCTIO &tio, yield_t & yield);

};

void SegTree9(MPCIO &mpcio, const PRACOptions &opts, char **args);
