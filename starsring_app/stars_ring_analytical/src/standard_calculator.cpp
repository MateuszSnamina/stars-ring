#include <stars_ring_analytical/standard_calculator.hpp>

#include <cassert>

namespace stars_ring_analytical {

AnalyticalFormulasBox::AnalyticalFormulasBox(
    std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system)
    : stars_ring_core::SettledInPhysicalSystem(physical_system) {}

double AnalyticalFormulasBox::ground_state_energy() const {
  return ground_state_classical_energy() + ground_state_correlation_energy();
}

StandardCalculator::StandardCalculator(
    std::shared_ptr<AnalyticalFormulasBox> formulas_box)
    : _formulas_box(formulas_box) {}

void StandardCalculator::calculate() {
  _ground_state_classical_energy =
      _formulas_box->ground_state_classical_energy();
  _ground_state_correlation_energy =
      _formulas_box->ground_state_correlation_energy();
  _ground_state_energy = _formulas_box->ground_state_energy();
  const unsigned max_tantamount_kn =
      _formulas_box->physical_system()->n_cells() / 2 + 1;
  _exc_state_relative_energy.resize(max_tantamount_kn);
  for (unsigned nk = 0; nk < max_tantamount_kn; ++nk)
    _exc_state_relative_energy[nk] =
        _formulas_box->exc_state_relative_energy(nk);
}

void StandardCalculator::formulas_box(
    std::shared_ptr<AnalyticalFormulasBox> formulas_box) {
  _formulas_box = formulas_box;
}

std::shared_ptr<AnalyticalFormulasBox> StandardCalculator::formulas_box()
    const {
  return _formulas_box;
}

double StandardCalculator::ground_state_classical_energy() const {
  return _ground_state_classical_energy;
}

double StandardCalculator::ground_state_correlation_energy() const {
  return _ground_state_correlation_energy;
}

double StandardCalculator::ground_state_energy() const {
  return _ground_state_energy;
}

unsigned StandardCalculator::n_exc_states() const {
  return _exc_state_relative_energy.size();
}

double StandardCalculator::exc_state_relative_energy(unsigned idx) const {
  assert(idx < _exc_state_relative_energy.size());
  return _exc_state_relative_energy[idx];
}

double StandardCalculator::exc_state_absolute_energy(unsigned idx) const {
  assert(idx < _exc_state_relative_energy.size());
  return _ground_state_energy + _exc_state_relative_energy[idx];
}

}  // namespace stars_ring_analytical