#include "types.hpp"
#include "mpcio.hpp"
#include "coroutine.hpp"
#include "options.hpp"
#include "mpcops.hpp"
#include "duoram.hpp"

class SegmentTree {
    private:
        Duoram < RegAS > oram;
        Duoram < RegXS > isEven;
        Duoram < RegXS > nextL;
        Duoram < RegXS > nextR;
        size_t num_items;

        void getBitVector(MPCTIO &tio, yield_t & yield, Duoram < RegXS > &bitVec, RegAS left, RegAS right);

    public:
        SegmentTree(int player_num, size_t size) : oram(player_num, size), isEven(player_num, size), nextL(player_num, size), nextR(player_num, size) {
            num_items = size;
        }

        void init(MPCTIO &tio, yield_t & yield);

        void RangeSum(MPCTIO &tio, yield_t & yield, RegAS left, RegAS right);
        // void Update(MPCTIO &tio, yield_t & yield, RegAS index, RegAS value);
    
};

void SegTree(MPCIO &mpcio, const PRACOptions &opts, char **args);