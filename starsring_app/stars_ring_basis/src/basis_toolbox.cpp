#include<cassert>
#include<armadillo>

#include<stars_ring_basis/transform_translation_unique_basis_to_ksubspace_basis.hpp>
#include<stars_ring_basis/perturbation_basis_maker.hpp>
#include<stars_ring_basis/localized_coupled_elements_generator.hpp>
#include<stars_ring_basis/translation_unique_coupled_elements_generator.hpp>
#include<stars_ring_basis/basis_toolbox.hpp>
#include<stars_ring_basis/perturbation_basis_maker.hpp>
#include<stars_ring_basis/build_perturbation_localized_basis.hpp>
#include<stars_ring_basis/build_perturbation_translation_unique_basis.hpp>

namespace stars_ring_basis {

    BasisBox::BasisBox(
            std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
            LocalizedBasis localized_basis,
            TranslationUniqueBasis translation_unique_basis,
            std::vector<KsubspaceBasis> ksubspace_basis) :
    SettledInPhysicalSystem(physical_system),
    _localized_basis(localized_basis),
    _translation_unique_basis(translation_unique_basis),
    _ksubspace_basis(ksubspace_basis) {
        assert(physical_system);
        assert(localized_basis.physical_system());
        assert(translation_unique_basis.physical_system());
        for (unsigned i = 0; i < ksubspace_basis.size(); i++)
            assert(ksubspace_basis[i].physical_system());
        assert(localized_basis.physical_system() == physical_system);
        assert(translation_unique_basis.physical_system() == physical_system);
        for (unsigned i = 0; i < ksubspace_basis.size(); i++)
            assert(ksubspace_basis[i].physical_system() == physical_system);
        assert(ksubspace_basis.size() == physical_system->n_cells());
    }

    BasisBox BasisBox::Builder::build(
            std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
            std::vector<RawState> init_states,
            std::shared_ptr<CoupledRawStatesGenerator> coupled_raw_states_generator,
            unsigned perturbation_level) {
        assert(physical_system);
        assert(coupled_raw_states_generator);        
        // Build localized_basis.
        std::shared_ptr<CoupledElementsGenerator < LocalizedElement>> coupled_localized_elements_generator
                = std::make_shared<CoupledLocalizedElementsGenerator>(coupled_raw_states_generator);
        const LocalizedBasis localized_basis =
                build_localized_basis(physical_system, init_states, coupled_localized_elements_generator, perturbation_level);
        // Build translation_unique_basis.
        std::shared_ptr<CoupledElementsGenerator < TranslationUniqueElement>> coupled_translation_unique_elements_generator
                = std::make_shared<CoupledTranslationUniqueElementsGenerator>(coupled_raw_states_generator);
        const TranslationUniqueBasis translation_unique_basis =
                build_translation_unique_basis(physical_system, init_states, coupled_translation_unique_elements_generator, perturbation_level);
        // Build ksubspace_basis.        
        const std::vector<KsubspaceBasis> ksubspace_basis =
                transform_translation_unique_basis_to_ksubspace_basis(physical_system, translation_unique_basis);
        // Return all the basis:
        return BasisBox(physical_system, localized_basis, translation_unique_basis, ksubspace_basis);
    }

    const LocalizedBasis& BasisBox::localized_basis() const {
        return _localized_basis;
    }

    const TranslationUniqueBasis& BasisBox::translation_unique_basis() const {
        return _translation_unique_basis;
    }

    const KsubspaceBasis& BasisBox::ksubspace_basis(unsigned nk) const {
        assert(nk < _ksubspace_basis.size());
        return _ksubspace_basis[nk];
    }

    const std::vector<KsubspaceBasis>& BasisBox::ksubspace_basis() const {
        return _ksubspace_basis;
    }

}