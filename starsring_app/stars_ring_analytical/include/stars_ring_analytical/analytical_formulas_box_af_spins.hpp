#ifndef ANALYTICAL_FORMULAS_BOX_AF_SPINS_HPP
#define ANALYTICAL_FORMULAS_BOX_AF_SPINS_HPP

#include <stars_ring_analytical/analytical_formulas_box_af_oscilators.hpp>
#include <stars_ring_analytical/standard_calculator.hpp>
#include <stars_ring_core/physical_system.hpp>

namespace stars_ring_analytical {

class AnalyticalFormulasBoxAfSpins : public AnalyticalFormulasBox {
 public:
  AnalyticalFormulasBoxAfSpins(
      std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
      unsigned multiplicity);
  double ground_state_classical_energy() const override;
  double ground_state_correlation_energy() const override;
  double exc_state_relative_energy(unsigned nk) const override;

 private:
  const unsigned _multiplicity;
  const AnalyticalFormulasBoxAfOscylators
      _analytical_formulas_box_af_oscylators;
};

}  // namespace stars_ring_analytical

#endif
