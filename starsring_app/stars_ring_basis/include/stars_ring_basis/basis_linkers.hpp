#ifndef BASIS_LINKERS_HPP
#define BASIS_LINKERS_HPP

#include<complex>
#include<vector>

#include<stars_ring_basis/linear_combination.hpp>
#include<stars_ring_basis/basis_typedefs.hpp>
#include<stars_ring_basis/basis_toolbox.hpp>

namespace stars_ring_basis {

    std::vector<unsigned> // the first idx is a KsubspaceElement index.
    construct_single_ksubspace_basis_to_translation_unique_basis_linker(
            const KsubspaceBasis & ksubsspace_basis,
            const TranslationUniqueBasis & translation_unique_basis);

    std::vector<std::vector<unsigned>> // the first idx is a TranslationUniqueElement index.
    construct_translation_unique_basis_to_localized_basis_linker(
            const TranslationUniqueBasis & translation_unique_basis,
            const LocalizedBasis & localized_basis);

    std::vector<std::vector<LinearCombination<std::complex<double>, unsigned>>> // the first idx corresponds to nk, the second is an KsubspaceElement index.
    construct_ksubspace_basis_to_localized_basis_linker(
            const BasisBox & basis_box);

// NOT Used any more:    
//    std::vector<std::vector<LinearCombination<std::complex<double>, unsigned>>>
//    construct_ksubspace_basis_to_localized_basis_linker(
//            const std::vector<KsubspaceBasis> & ksubspace_basis,
//            const LocalizedBasis & localized_basis);
    
}

#endif