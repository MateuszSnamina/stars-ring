#ifndef CYCLE_HELPER_FUNCTIONS_HPP
#define CYCLE_HELPER_FUNCTIONS_HPP

#include <tuple>
#include <vector>

#include <stars_ring_basis/raw_state.hpp>

namespace stars_ring_basis {

std::tuple<RawState, unsigned> determine_representative_and_shift(RawState v);
RawState determine_representative(RawState v);
unsigned determine_cycle_length(const RawState& v);
std::vector<RawState> determine_cycle_states(const RawState& v);

}  // namespace stars_ring_basis

#endif