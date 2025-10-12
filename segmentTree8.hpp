#include "types.hpp"
#include "mpcio.hpp"
#include "coroutine.hpp"
#include "options.hpp"
#include "mpcops.hpp"
#include "duoram.hpp"

// Forward declaration
class PerformanceLogger;

// Segment Tree Node structure with named fields
struct SegTreeNode {
    RegAS value;
    RegAS parent;
    RegAS leftChildSibling;
    RegAS rightChildSibling;

    // Field-access macros for convenience (using unique names to avoid conflicts)
    #define SEGTREE_VALUE field(&SegTreeNode::value)
    #define SEGTREE_PARENT field(&SegTreeNode::parent)
    #define SEGTREE_LEFT_SIBLING field(&SegTreeNode::leftChildSibling)
    #define SEGTREE_RIGHT_SIBLING field(&SegTreeNode::rightChildSibling)

    // For debugging
    void dump() const {
        printf("[val=%016lx par=%016lx lsib=%016lx rsib=%016lx]", 
               value.share(), parent.share(), 
               leftChildSibling.share(), rightChildSibling.share());
    }

    // Required operations for Duoram
    inline void randomize() {
        value.randomize();
        parent.randomize();
        leftChildSibling.randomize();
        rightChildSibling.randomize();
    }

    inline SegTreeNode &operator+=(const SegTreeNode &rhs) {
        this->value += rhs.value;
        this->parent += rhs.parent;
        this->leftChildSibling += rhs.leftChildSibling;
        this->rightChildSibling += rhs.rightChildSibling;
        return *this;
    }

    inline SegTreeNode operator+(const SegTreeNode &rhs) const {
        SegTreeNode res = *this;
        res += rhs;
        return res;
    }

    inline SegTreeNode &operator-=(const SegTreeNode &rhs) {
        this->value -= rhs.value;
        this->parent -= rhs.parent;
        this->leftChildSibling -= rhs.leftChildSibling;
        this->rightChildSibling -= rhs.rightChildSibling;
        return *this;
    }

    inline SegTreeNode operator-(const SegTreeNode &rhs) const {
        SegTreeNode res = *this;
        res -= rhs;
        return res;
    }

    inline SegTreeNode operator-() const {
        SegTreeNode res;
        res.value = -this->value;
        res.parent = -this->parent;
        res.leftChildSibling = -this->leftChildSibling;
        res.rightChildSibling = -this->rightChildSibling;
        return res;
    }

    // Multiply each field by the local share of the corresponding field
    inline SegTreeNode mulshare(const SegTreeNode &rhs) const {
        SegTreeNode res = *this;
        res.value.mulshareeq(rhs.value);
        res.parent.mulshareeq(rhs.parent);
        res.leftChildSibling.mulshareeq(rhs.leftChildSibling);
        res.rightChildSibling.mulshareeq(rhs.rightChildSibling);
        return res;
    }

    // Turn DPF leaf into unit element
    template <nbits_t WIDTH>
    inline void unit(const RDPF<WIDTH> &dpf,
        typename RDPF<WIDTH>::LeafNode leaf) {
        value = dpf.unit_as(leaf);
        parent = dpf.unit_as(leaf);
        leftChildSibling = dpf.unit_as(leaf);
        rightChildSibling = dpf.unit_as(leaf);
    }

    // Perform an update on each of the fields, using field-specific
    // MemRefs constructed from the Shape shape and the index idx
    template <typename Sh, typename U>
    inline static void update(Sh &shape, yield_t &shyield, U idx,
            const SegTreeNode &M) {
        run_coroutines(shyield,
            [&shape, &idx, &M] (yield_t &yield) {
                Sh Sh_coro = shape.context(yield);
                Sh_coro[idx].SEGTREE_VALUE += M.value;
            },
            [&shape, &idx, &M] (yield_t &yield) {
                Sh Sh_coro = shape.context(yield);
                Sh_coro[idx].SEGTREE_PARENT += M.parent;
            },
            [&shape, &idx, &M] (yield_t &yield) {
                Sh Sh_coro = shape.context(yield);
                Sh_coro[idx].SEGTREE_LEFT_SIBLING += M.leftChildSibling;
            },
            [&shape, &idx, &M] (yield_t &yield) {
                Sh Sh_coro = shape.context(yield);
                Sh_coro[idx].SEGTREE_RIGHT_SIBLING += M.rightChildSibling;
            });
    }
};

// I/O operations (for sending over the network)
template <typename T>
T& operator>>(T& is, SegTreeNode &x)
{
    is >> x.value >> x.parent >> x.leftChildSibling >> x.rightChildSibling;
    return is;
}

template <typename T>
T& operator<<(T& os, const SegTreeNode &x)
{
    os << x.value << x.parent << x.leftChildSibling << x.rightChildSibling;
    return os;
}

class SegmentTree8 {
    private:
    // Single Duoram storing SegTreeNode structs
    Duoram<SegTreeNode> oram;

    RegAS computeRangeSum(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, RegAS left, RegAS right, 
                          PerformanceLogger* logger = nullptr, size_t operation_id = 0);

    public:

    size_t num_items;
    size_t depth;

    SegmentTree8(int player_num, size_t size, size_t d) : oram(player_num, size) {
            num_items = size;
            depth = d;
            std:: cout << "Segment Tree of depth " << depth << " with " << num_items << " nodes created" << std::endl;
        }

        void init(MPCTIO &tio, yield_t & yield);

    void RangeSum(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, RegAS left, RegAS right, 
                  PerformanceLogger* logger = nullptr, size_t operation_id = 0);
    // Instrumented Update (global-index version): measures (1) leaf read, (2) per-level update write, (3) parent index read timings/resources
    void Update(MPCTIO &tio, MPCIO &mpcio, yield_t & yield, RegAS index, RegAS value,
                PerformanceLogger* logger = nullptr, size_t operation_id = 0);
    void printSegmentTree(MPCTIO &tio, yield_t & yield);

};

void SegTree8(MPCIO &mpcio, const PRACOptions &opts, char **args);
