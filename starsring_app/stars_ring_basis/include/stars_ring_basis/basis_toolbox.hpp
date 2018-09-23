#ifndef BASIS_TOOLBOX_HPP
#define BASIS_TOOLBOX_HPP

#include <stars_ring_basis/basis_typedefs.hpp>
#include <stars_ring_basis/raw_state.hpp>
#include <stars_ring_basis/raw_state_coupled_elements_generator.hpp>
#include <stars_ring_core/physical_system.hpp>

namespace stars_ring_basis {

class BasisBox : public stars_ring_core::SettledInPhysicalSystem {
 public:
  const LocalizedBasis& localized_basis() const;
  const TranslationUniqueBasis& translation_unique_basis() const;
  const KsubspaceBasis& ksubspace_basis(unsigned nk) const;
  const std::vector<KsubspaceBasis>& ksubspace_basis() const;
  class Builder;

 private:
  BasisBox(std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
           LocalizedBasis localized_basis,
           TranslationUniqueBasis translation_uniqueBasis,
           std::vector<KsubspaceBasis> ksubspace_basis);
  const LocalizedBasis _localized_basis;
  const TranslationUniqueBasis _translation_unique_basis;
  const std::vector<KsubspaceBasis> _ksubspace_basis;
};

class BasisBox::Builder {
 public:
  BasisBox build(
      std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
      std::vector<RawState> init_states,
      std::shared_ptr<CoupledRawStatesGenerator> coupled_raw_states_generator,
      unsigned perturbation_level);
};

}  // namespace stars_ring_basis

#endif