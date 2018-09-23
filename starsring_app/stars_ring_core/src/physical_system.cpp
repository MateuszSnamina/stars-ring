#include <armadillo>
#include <cassert>

#include <stars_ring_core/physical_system.hpp>

namespace {
const std::complex<double> imag_unit(0, 1);
}

namespace stars_ring_core {

PhysicalSystem::PhysicalSystem(unsigned n_cells, unsigned max_n_stars)
    : _n_sites(n_cells * 2),
      _n_cells(n_cells),
      _n_sites_in_cell(2),
      _max_n_stars(max_n_stars),
      _k_unit1(2 * arma::datum::pi / _n_sites),
      _k_unit2(2 * arma::datum::pi / _n_cells),
      _chi(std::exp(imag_unit * _k_unit2)) {
  assert(n_cells > 0);
}

PhysicalSystem::Builder::Builder() : _n_cells(4), _max_n_stars(2) {}

PhysicalSystem::Builder& PhysicalSystem::Builder::n_cells(unsigned n_cells) {
  assert(n_cells > 0);
  _n_cells = n_cells;
  return *this;
}

PhysicalSystem::Builder& PhysicalSystem::Builder::max_n_stars(
    unsigned max_n_stars) {
  _max_n_stars = max_n_stars;
  return *this;
}

std::unique_ptr<PhysicalSystem> PhysicalSystem::Builder::build() {
  return std::unique_ptr<PhysicalSystem>(
      new PhysicalSystem(_n_cells, _max_n_stars));
}

SettledInPhysicalSystem::SettledInPhysicalSystem(
    std::shared_ptr<PhysicalSystem> physical_system)
    : _physical_system(physical_system) {
  assert(physical_system);
}

std::ostream& operator<<(std::ostream& s,
                         const PhysicalSystem& physical_system) {
  s << "PhysicalSystem("
    << "n_cells = " << physical_system.n_cells() << ", "
    << "max_n_stars = " << physical_system.max_n_stars() << ")";
  return s;
}

}  // namespace stars_ring_core