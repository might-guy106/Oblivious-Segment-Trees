#include "types.hpp"
#include "mpcio.hpp"
#include "coroutine.hpp"
#include "options.hpp"
#include "mpcops.hpp"
#include "duoram.hpp"
#include <vector>

class SegmentTree {
    private:
        Duoram < RegAS > oram;
        Duoram < RegBS > isLeft;
        Duoram < RegBS > isRight;
        Duoram < RegAS > nextL;
        Duoram < RegAS > nextR;
        size_t num_items;

    public:
    SegmentTree(int player_num, size_t size)
    : oram(player_num, size + 1),
      isLeft(player_num, size + 1),
      isRight(player_num, size + 1),
      nextL(player_num, size + 1),
      nextR(player_num, size + 1) {num_items = size; }//size is a power of 2 - 1

        void init(MPCTIO &tio, yield_t & yield, size_t n);

        Duoram<RegBS> computeBitVector(MPCTIO &tio, yield_t & yield, RegAS left_index, RegAS right_index); //left_index and right_index are in leaf layer (starting from (size+1)/2)

        RegAS rangeSum(MPCTIO &tio, yield_t & yield, RegAS left_index, RegAS right_index);

};

void SegTree(MPCIO &mpcio, const PRACOptions &opts, char **args);

