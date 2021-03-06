#ifndef RAW_STATE_HPP
#define RAW_STATE_HPP

#include <iostream>
#include <vector>

namespace stars_ring_basis {

typedef std::vector<unsigned> RawState;

std::ostream& operator<<(std::ostream& s, const RawState& el);

}  // namespace stars_ring_basis

#endif