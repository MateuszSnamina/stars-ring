//#include<armadillo>
#include<vector>

#include<stars_ring_basis/transform_translation_unique_basis_to_ksubspace_basis.hpp>

namespace stars_ring_basis {

    KsubspaceBasis transform_translation_unique_basis_to_ksubspace_basis(
            const TranslationUniqueBasis& translation_unique_basis,
            unsigned nk) {
        KsubspaceBasis result(translation_unique_basis.physical_system());
        for (std::shared_ptr<TranslationUniqueElement> translation_unique_element : translation_unique_basis.vec_index()) {
            assert(translation_unique_element);
            if (translation_unique_element->support_nk(nk)) {
                result.add_element(make_ksubspace_element(translation_unique_element, nk));
            }
        }
        return result;
    }

    std::vector<KsubspaceBasis> transform_translation_unique_basis_to_ksubspace_basis(
            std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
            const TranslationUniqueBasis& translation_unique_basis) {
        //std::cout << "[INFO   ] [PROGRESS] Program is about to construct k-subspace-basis" << std::endl;
        //arma::wall_clock timer;
        //timer.tic();
        // -------------------------------
        assert(physical_system);
        // -------------------------------
        std::vector<KsubspaceBasis> ksubspace_basis;
        ksubspace_basis.reserve(physical_system->n_cells());
        for (unsigned nk = 0; nk < physical_system->n_cells(); nk++)
            ksubspace_basis.push_back(transform_translation_unique_basis_to_ksubspace_basis(translation_unique_basis, nk));
        // -------------------------------        
        //const double calculations_time = timer.toc();
        //std::cout << "[INFO   ] [PROGRESS] Program has constructed k-subspace-basis"
        //        << " (in " << calculations_time << "s)." << std::endl;
        //---------------------        
        return ksubspace_basis;
    }

}