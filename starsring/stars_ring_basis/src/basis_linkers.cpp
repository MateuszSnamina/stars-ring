
#include<armadillo>
#include<stars_ring_basis/basis_linkers.hpp>

namespace stars_ring_basis {

    std::vector<unsigned>
    construct_single_ksubspace_basis_to_translation_unique_basis_linker(
            const KsubspaceBasis & ksubsspace_basis,
            const TranslationUniqueBasis & translation_unique_basis) {
        //---------------------
        assert(translation_unique_basis.physical_system());
        assert(ksubsspace_basis.physical_system());
        assert(ksubsspace_basis.physical_system() == translation_unique_basis.physical_system());
        //---------------------
        std::cout << "[INFO   ] [PROGRESS] Program is about to construct k-subspace-basis-to-translation-unique-basis-linker." << std::endl;
        arma::wall_clock timer;
        timer.tic();
        //---------------------
        std::vector<unsigned> single_ksubspace_basis_to_translation_unique_basis_linker;
        single_ksubspace_basis_to_translation_unique_basis_linker.reserve(ksubsspace_basis.size());
        for (std::shared_ptr<KsubspaceElement> ksubsspace_element : ksubsspace_basis.vec_index()) {
            assert(ksubsspace_element);
            const RawState raw_ksubsspace_element = ksubsspace_element->translation_unique_element()->raw_representative_state();
            boost::optional<unsigned> raw_ksubsspace_element_idx = translation_unique_basis.find_element_and_get_its_ra_index(raw_ksubsspace_element);
            assert(raw_ksubsspace_element_idx);
            assert(*raw_ksubsspace_element_idx < translation_unique_basis.size());
            single_ksubspace_basis_to_translation_unique_basis_linker.push_back(*raw_ksubsspace_element_idx);
        }
        //---------------------        
        const double calculations_time = timer.toc();
        std::cout << "[INFO   ] [PROGRESS] Program has constructed k-subspace-basis-to-translation-unique-basis-linker"
                << " (in " << calculations_time << "s)." << std::endl;
        //---------------------        
        return single_ksubspace_basis_to_translation_unique_basis_linker;
    }

    std::vector<std::vector<unsigned>>
    construct_translation_unique_basis_to_localized_basis_linker(
            const TranslationUniqueBasis & translation_unique_basis,
            const LocalizedBasis & localized_basis) {
        //---------------------
        assert(localized_basis.physical_system());
        assert(translation_unique_basis.physical_system());
        assert(translation_unique_basis.physical_system() == localized_basis.physical_system());
        //---------------------
        std::cout << "[INFO   ] [PROGRESS] Program is about to construct translation-unique-basis-to-localized-basis-linker." << std::endl;
        arma::wall_clock timer;
        timer.tic();
        //---------------------
        std::vector<std::vector<unsigned>> translation_unique_basis_to_localized_basis_linker;
        translation_unique_basis_to_localized_basis_linker.reserve(translation_unique_basis.size());
        for (std::shared_ptr<TranslationUniqueElement> tranlation_unique_element : translation_unique_basis.vec_index()) {
            assert(tranlation_unique_element);
            std::vector<unsigned> equivalent_states_indices;
            equivalent_states_indices.reserve(tranlation_unique_element->cycle_length());
            for (const RawState & raw_localized_element : tranlation_unique_element->raw_equivalent_states()) {
                boost::optional<unsigned> raw_localized_element_idx = localized_basis.find_element_and_get_its_ra_index(raw_localized_element);
                assert(raw_localized_element_idx);
                assert(*raw_localized_element_idx < localized_basis.size());
                equivalent_states_indices.push_back(*raw_localized_element_idx);
            }
            translation_unique_basis_to_localized_basis_linker.push_back(equivalent_states_indices);
        }
        //---------------------        
        const double calculations_time = timer.toc();
        std::cout << "[INFO   ] [PROGRESS] Program has constructed translation-unique-basis-to-localized-basis-linker"
                << " (in " << calculations_time << "s)." << std::endl;
        //---------------------        
        return translation_unique_basis_to_localized_basis_linker;
    }

    std::vector<std::vector<LinearCombination<std::complex<double>, unsigned>>>
    construct_ksubspace_basis_to_localized_basis_linker(
            const BasisBox & basis_box) {
        const std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system = basis_box.physical_system();
        //---------------------
        std::cout << "[INFO   ] [PROGRESS] Program is about to construct k-subspace-basis-to-localized-basis-linker." << std::endl;
        arma::wall_clock timer;
        timer.tic();
        //---------------------
        std::vector<std::vector<unsigned>> translation_unique_basis_to_localized_basis_linker =
                construct_translation_unique_basis_to_localized_basis_linker(
                basis_box.translation_unique_basis(),
                basis_box.localized_basis());

        std::vector<std::vector<LinearCombination < std::complex<double>, unsigned>>> ksubspace_basis_to_localized_basis_linker(physical_system->n_cells()); // the first index is for kn, the second index tantamount to ra index KsubspaceBasis.
        for (unsigned nk = 0; nk < physical_system->n_cells(); nk++) {
            std::vector<unsigned> single_ksubspace_basis_to_translation_unique_basis_linker =
                    construct_single_ksubspace_basis_to_translation_unique_basis_linker(
                    basis_box.ksubspace_basis(nk),
                    basis_box.translation_unique_basis());
            std::complex<double> amplitude_roration = std::exp((nk * physical_system->k_unit2()) * std::complex<double>(0, 1));
            ksubspace_basis_to_localized_basis_linker[nk].reserve(basis_box.ksubspace_basis(nk).size());
            for (unsigned k_subspace_element_index = 0; k_subspace_element_index < basis_box.ksubspace_basis(nk).size(); k_subspace_element_index++) {
                const unsigned translation_unique_element_index = single_ksubspace_basis_to_translation_unique_basis_linker[k_subspace_element_index];
                const std::vector<unsigned> localized_elements_indices = translation_unique_basis_to_localized_basis_linker[translation_unique_element_index];
                const unsigned number_of_localized_elements = localized_elements_indices.size();
                std::complex<double> amplitude = 1.0 / std::sqrt(number_of_localized_elements);
                LinearCombination < std::complex<double>, unsigned> linear_combination;
                for (unsigned localized_element_index : localized_elements_indices) {
                    assert(localized_element_index < basis_box.localized_basis().size());
                    LinearCombinationIngredient< std::complex<double>, unsigned> linear_combination_ingredient
                            = {amplitude, localized_element_index};
                    linear_combination.push_back(linear_combination_ingredient);
                    amplitude *= amplitude_roration;
                }
                assert(std::arg(amplitude) < 1e-5);
                ksubspace_basis_to_localized_basis_linker[nk].push_back(linear_combination);
            }
        }
        //---------------------        
        const double calculations_time = timer.toc();
        std::cout << "[INFO   ] [PROGRESS] Program has constructed k-subspace-basis-to-localized-basis-linker"
                << " (in " << calculations_time << "s)." << std::endl;
        //---------------------        
        return ksubspace_basis_to_localized_basis_linker;
    }

    /*
    std::vector<std::vector<LinearCombination<std::complex<double>, unsigned>>>
    construct_ksubspace_basis_to_localized_basis_linker(
            const std::vector<KsubspaceBasis> & ksubspace_basis,
            const LocalizedBasis & localized_basis) {
        //---------------------
        assert(localized_basis.physical_system());
        for (unsigned nk = 0; nk < ksubspace_basis.size(); nk++) {
            assert(ksubspace_basis[nk].physical_system());
            assert(ksubspace_basis[nk].physical_system() == localized_basis.physical_system());
        }
        const std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system = localized_basis.physical_system();
        assert(ksubspace_basis.size() == physical_system->n_cells());
        //---------------------
        std::cout << "[INFO   ] [PROGRESS] Program is about to construct k-subspace-basis-to-localized-basis-linker [OLD]. " << std::endl;
        arma::wall_clock timer;
        timer.tic();
        std::vector<std::vector<LinearCombination < std::complex<double>, unsigned>>> ksubspace_basis_to_localized_basis_linker(physical_system->n_cells()); // the first index is for kn, the second index tantamount to ra index KsubspaceBasis.
        for (unsigned nk = 0; nk < physical_system->n_cells(); nk++) {
            std::complex<double> amplitude_roration = std::exp((nk * physical_system->k_unit2()) * std::complex<double>(0, 1));
            ksubspace_basis_to_localized_basis_linker[nk].reserve(ksubspace_basis[nk].size());
            for (std::shared_ptr<KsubspaceElement> k_subspace_element : ksubspace_basis[nk].vec_index()) {
                std::complex<double> amplitude = 1.0 / std::sqrt(k_subspace_element->translation_unique_element()->cycle_length());
                assert(k_subspace_element);
                LinearCombination < std::complex<double>, unsigned> linear_combination;
                for (const RawState & raw_localized_element : k_subspace_element->translation_unique_element()->raw_equivalent_states()) {
                    boost::optional<unsigned> raw_localized_element_idx = localized_basis.find_element_and_get_its_ra_index(raw_localized_element);
                    assert(raw_localized_element_idx);
                    assert(*raw_localized_element_idx < localized_basis.size());
                    LinearCombinationIngredient< std::complex<double>, unsigned> linear_combination_ingredient
                            = {amplitude, *raw_localized_element_idx};
                    linear_combination.push_back(linear_combination_ingredient);
                    amplitude *= amplitude_roration;
                }
                assert(std::arg(amplitude) < 1e-5);
                ksubspace_basis_to_localized_basis_linker[nk].push_back(linear_combination);
            }
        }
        //---------------------        
        const double calculations_time = timer.toc();
        std::cout << "[INFO   ] [PROGRESS] Program has constructed k-subspace-basis-to-localized-basis-linker [OLD]"
                << " (in " << calculations_time << "s)." << std::endl;
        //---------------------        
        return ksubspace_basis_to_localized_basis_linker;
    }
     */


}