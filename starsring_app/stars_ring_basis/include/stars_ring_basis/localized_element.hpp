#ifndef LOCALIZED_ELEMENT_HPP
#define LOCALIZED_ELEMENT_HPP

#include <iostream>
#include <memory>

#include <stars_ring_basis/raw_state.hpp>
#include <stars_ring_core/physical_system.hpp>

namespace stars_ring_basis {

class LocalizedElement : public stars_ring_core::SettledInPhysicalSystem {
 public:
  // Physical properties getters:
  const RawState& raw_state() const;
  unsigned perturbation_order() const;
  unsigned n_stars_A() const;
  unsigned n_stars_B() const;
  unsigned n_stars() const;
  // Factory:
  friend std::unique_ptr<LocalizedElement> make_localized_element(
      std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
      const RawState& state, unsigned perturbation_order);
  // Member function that enables
  // a class instance to be stored in VecMap
  typedef const RawState& KeyT;
  const KeyT key() const;

 private:
  LocalizedElement(
      std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
      const RawState& state, unsigned perturbation_order);
  const RawState _raw_state;
  const unsigned _perturbation_order;
  const unsigned _n_stars_A;
  const unsigned _n_stars_B;
  const unsigned _n_stars;
};

std::unique_ptr<LocalizedElement> make_localized_element(
    std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
    const RawState& state, unsigned perturbation_order);

std::ostream& operator<<(std::ostream& s, const LocalizedElement& el);

// Inline functions implementations:

inline const RawState& LocalizedElement::raw_state() const {
  return _raw_state;
}

inline unsigned LocalizedElement::perturbation_order() const {
  return _perturbation_order;
}

inline typename LocalizedElement::KeyT LocalizedElement::key() const {
  return _raw_state;
}

inline unsigned LocalizedElement::n_stars_A() const { return _n_stars_A; }

inline unsigned LocalizedElement::n_stars_B() const { return _n_stars_B; }

inline unsigned LocalizedElement::n_stars() const { return _n_stars; }

}  // namespace stars_ring_basis

#endif
