#include <stars_ring_basis/raw_state_operations.hpp>

namespace stars_ring_basis {
namespace operations {

RawState make_ground_state(unsigned n_sites) {
  RawState result(n_sites);
  return result;
}

RawState make_one_star_state(unsigned n_sites, unsigned site) {
  assert(site < n_sites);
  RawState result(n_sites);
  result[site] = 1;
  return result;
}

double determine_coupling(const RawState& bra, const RawState& ket,
                          const maths_in_physic::SpinRealm& spin_realm) {
  assert(bra.size() == ket.size());
  unsigned size = bra.size();
  assert(size % 2 == 0);
  double result = 0.0;
  for (unsigned i = 0, j = 1; i < size; i++, j = (j + 1) % size) {
    {
      if (bra == ket) {
        result -=
            spin_realm.S_z(bra[i], ket[i]) * spin_realm.S_z(bra[j], ket[j]);
        continue;  // only for performance optimization.
      }
    }
    {
      const RawState ket_pp = _anan(i, j, ket);
      if (bra == ket_pp) {
        result += spin_realm.S_p(bra[i], ket[i]) *
                  spin_realm.S_p(bra[j], ket[j]) / 2.0;
        break;  // only for performance optimization.
      }
    }
    {
      const RawState ket_mm = _crcr(i, j, ket);
      if (bra == ket_mm) {
        result += spin_realm.S_m(bra[i], ket[i]) *
                  spin_realm.S_m(bra[j], ket[j]) / 2.0;
        break;  // only for performance optimization.
      }
    }
  }
  return result;
}

double determine_coupling(
    const RawState& bra, const RawState& ket,
    const maths_in_physic::OscilatorRealm& oscilatorRealm) {
  assert(false);  // not implemented
  return 0.0;
}

}  // namespace operations
}  // namespace stars_ring_basis
