#ifndef RAW_STATE_OPERATIONS_HPP
#define RAW_STATE_OPERATIONS_HPP

#include <boost/optional.hpp>
#include <vector>

#include <maths_in_physic/oscilator_realm.hpp>
#include <maths_in_physic/spin_realm.hpp>
#include <stars_ring_basis/raw_state.hpp>

namespace stars_ring_basis {
namespace operations {

inline unsigned n_stars_A(const RawState &state) {
  unsigned result = 0;
  for (unsigned i = 0; i < state.size(); i = i + 2) result += state[i];
  return result;
}

inline unsigned n_stars_B(const RawState &state) {
  unsigned result = 0;
  for (unsigned i = 1; i < state.size(); i = i + 2) result += state[i];
  return result;
}

inline unsigned n_stars(const RawState &state) {
  return std::accumulate(state.begin(), state.end(), 0.0);
}

inline RawState _crcr(unsigned i, unsigned j, RawState state) {
  state[i]++;
  state[j]++;
  return state;
}

inline RawState _anan(unsigned i, unsigned j, RawState state) {
  state[i]--;
  state[j]--;
  return state;
}

inline boost::optional<RawState> crcr(unsigned i, unsigned j,
                                      const RawState &state,
                                      unsigned max_n_stars) {
  if (state[i] + 1 < max_n_stars && state[j] + 1 < max_n_stars)
    return _crcr(i, j, state);
  return boost::optional<RawState>();
}

inline boost::optional<RawState> anan(unsigned i, unsigned j,
                                      const RawState &state) {
  if (state[i] > 0 && state[j] > 0) return _anan(i, j, state);
  return boost::optional<RawState>();
}

RawState make_ground_state(unsigned n_sites);
RawState make_one_star_state(unsigned n_sites, unsigned site);

double determine_coupling(const RawState &bra, const RawState &ket,
                          const maths_in_physic::SpinRealm &spin_realm);

double determine_coupling(
    const RawState &bra, const RawState &ket,
    const maths_in_physic::OscilatorRealm &oscilatorRealm);

}  // namespace operations
}  // namespace stars_ring_basis

#endif