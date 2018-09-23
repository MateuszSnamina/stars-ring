#ifndef BUILD_PERTURBATION_TRANSLATION_UNIQUE_BASIS_HPP
#define BUILD_PERTURBATION_TRANSLATION_UNIQUE_BASIS_HPP

#include <stars_ring_basis/basis_typedefs.hpp>
#include <stars_ring_basis/perturbation_basis_maker.hpp>
#include <stars_ring_basis/raw_state.hpp>
#include <stars_ring_core/physical_system.hpp>

namespace stars_ring_basis {

TranslationUniqueBasis build_translation_unique_basis(
    std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
    std::vector<RawState> init_states,
    std::shared_ptr<CoupledElementsGenerator<TranslationUniqueElement>>
        coupled_translation_unique_elements_generator,
    unsigned perturbation_level);
}

#endif