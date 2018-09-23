#include<cassert>
#include<iomanip>

#include<stars_ring_basis/raw_state.hpp>

namespace stars_ring_basis {

    std::ostream& operator<<(
            std::ostream& s,
            const RawState & el) {
        auto f = s.flags();
        s << "〖";
        for (unsigned val : el)
            std::cout << std::setw(2) << val;
        s << "〗";
        s.flags(f);
        return s;
    }

}