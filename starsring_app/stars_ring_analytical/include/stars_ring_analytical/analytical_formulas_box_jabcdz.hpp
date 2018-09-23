#ifndef ANALYTICAL_FORMULAS_BOX_JABCDZ_HPP
#define ANALYTICAL_FORMULAS_BOX_JABCDZ_HPP

#include <memory>

#include <stars_ring_analytical/analytical_formulas_box_af_spins.hpp>
#include <stars_ring_analytical/standard_calculator.hpp>
#include <stars_ring_core/physical_system.hpp>

namespace stars_ring_analytical {

class AnalyticalFormulasBoxJABCDZ : public AnalyticalFormulasBox {
 public:
  AnalyticalFormulasBoxJABCDZ(
      std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
      unsigned multiplicity, double A, double B, double C, double D, double J,
      double Ez, double theta);
  AnalyticalFormulasBoxJABCDZ(
      std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
      unsigned multiplicity, double A, double B, double C, double D, double J,
      double Ez);

  double ground_state_classical_energy() const override;
  double ground_state_correlation_energy() const override;
  double exc_state_relative_energy(unsigned nk) const override;
  double theta() const;
  double mean_orbital_operator() const;
  double J_spin() const;
  double theta_opt_nell() const;
  double theta_opt_corrected_nell() const;

 private:
  AnalyticalFormulasBoxAfSpins _analytical_formulas_box_af_spins;
  const unsigned _multiplicity;
  const double _A, _B, _C, _D;
  const double _J, _Ez;
  const double _theta;
};

}  // namespace stars_ring_analytical
#endif
