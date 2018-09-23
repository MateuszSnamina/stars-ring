#include<cassert>
#include<algorithm>

#include<stars_ring_basis/cycle_helper_functions.hpp>

namespace stars_ring_basis {

    std::tuple<RawState, unsigned> determine_representative_and_shift(RawState v) {
        assert(v.size() % 2 == 0);
        const unsigned step = 2;
        RawState representative = v;
        unsigned shift = 0;
        for (unsigned i = 1; i < v.size() / step; i++) {
            std::rotate(v.begin(), v.begin() + step, v.end());
            if (v > representative) {
                representative = v;
                shift = i;
            }
        }
        return std::make_tuple(representative, shift);
    }

    RawState determine_representative(RawState v) {
        return std::get<0>(determine_representative_and_shift(v));
    }

    unsigned determine_cycle_length(const RawState& v) {
        assert(v.size() % 2 == 0);
        const unsigned step = 2;
        RawState temp = v;
        for (unsigned i = 1;; i++) {
            std::rotate(temp.begin(), temp.begin() + step, temp.end());
            if (std::equal(temp.begin(), temp.end(), v.begin()))
                return i;
        }
    }

    std::vector<RawState> determine_cycle_states(const RawState& v) {
        assert(v.size() % 2 == 0);
        const unsigned step = 2;
        std::vector<RawState> result;
        result.push_back(v);
        while (true) {
            RawState temp = result.back();
            std::rotate(temp.begin(), temp.begin() + step, temp.end());
            if (std::equal(temp.begin(), temp.end(), result.front().begin()))
                break;
            result.push_back(temp);
        }
        return result;
    }

}