#ifndef TRANSLATION_UNIQUE_BASIS_TO_KSUBSPACE_BASIS_HPP
#define TRANSLATION_UNIQUE_BASIS_TO_KSUBSPACE_BASIS_HPP

#include<stars_ring_basis/basis_typedefs.hpp>

namespace stars_ring_basis {

    KsubspaceBasis transform_translation_unique_basis_to_ksubspace_basis(
            const TranslationUniqueBasis& translation_unique_basis,
            unsigned nk);

    std::vector<KsubspaceBasis> transform_translation_unique_basis_to_ksubspace_basis(
            std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
            const TranslationUniqueBasis& translation_unique_basis);

}

#endif