#include <cassert>
#include <cmath>
#include <iomanip>

#include <stars_ring_basis/cycle_helper_functions.hpp>
#include <stars_ring_basis/raw_state_operations.hpp>
#include <stars_ring_basis/translation_unique_element.hpp>

namespace stars_ring_basis {

TranslationUniqueElement::TranslationUniqueElement(
    std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
    const RawState& state, unsigned perturbation_order)
    : SettledInPhysicalSystem(physical_system),
      _raw_representative_state(determine_representative(state)),
      //_raw_equivalent_states(determine_cycle_states(raw_representative_state())),
      //_cycle_length(raw_equivalent_states().size()),
      _cycle_length(determine_cycle_length(state)),
      _perturbation_order(perturbation_order),
      _n_stars_A(operations::n_stars_A(_raw_representative_state)),
      _n_stars_B(operations::n_stars_B(_raw_representative_state)),
      _n_stars(operations::n_stars(_raw_representative_state)) {
  assert(physical_system);
  assert(state.size() == physical_system->n_sites());
}

const std::vector<RawState> TranslationUniqueElement::raw_equivalent_states()
    const {
  return determine_cycle_states(raw_representative_state());
  // return _raw_equivalent_states;
}

bool TranslationUniqueElement::support_nk(unsigned nk) const {
  assert(raw_representative_state().size() % 2 == 0);
  assert(raw_representative_state().size() == physical_system()->n_sites());
  assert(nk < physical_system()->n_cells());
  assert(physical_system()->n_cells() % cycle_length() == 0);
  return nk % (physical_system()->n_cells() / cycle_length()) == 0;
}

double TranslationUniqueElement::norm(unsigned nk) const {
  if (!support_nk(nk)) return 0;
  return physical_system()->n_cells() / std::sqrt(cycle_length());
}

std::ostream& operator<<(std::ostream& s, const TranslationUniqueElement& el) {
  auto f = s.flags();
  s << el.raw_representative_state();
  s << " ";
  s << "〘" << std::setw(2) << el.perturbation_order() << "⚡〙";
  s << "〘" << std::setw(2) << el.cycle_length() << "♲〙";
  s.flags(f);
  return s;
}

std::unique_ptr<TranslationUniqueElement> make_translation_unique_element(
    std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
    const RawState& state, unsigned perturbation_order) {
  assert(physical_system);
  assert(state.size() == physical_system->n_sites());
  TranslationUniqueElement* raw_ptr =
      new TranslationUniqueElement(physical_system, state, perturbation_order);
  return std::unique_ptr<TranslationUniqueElement>(raw_ptr);
}

}  // namespace stars_ring_basis