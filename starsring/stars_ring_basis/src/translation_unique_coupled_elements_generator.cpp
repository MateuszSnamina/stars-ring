#include<cassert>

#include<stars_ring_basis/translation_unique_coupled_elements_generator.hpp>
#include<stars_ring_basis/raw_state_operations.hpp>

namespace stars_ring_basis {

    CoupledTranslationUniqueElementsGenerator::CoupledTranslationUniqueElementsGenerator(
            std::shared_ptr<CoupledRawStatesGenerator> coupled_raw_states_generator) :
    _coupled_raw_states_generator(coupled_raw_states_generator) {
        assert(coupled_raw_states_generator);
    }

    std::vector<std::shared_ptr<TranslationUniqueElement>>
    CoupledTranslationUniqueElementsGenerator::generate(std::shared_ptr<TranslationUniqueElement> element) const {
        assert(element);
        assert(element->physical_system());
        assert(_coupled_raw_states_generator);
        const std::vector<RawState> raw_coupled_states =
                _coupled_raw_states_generator->generate(element->raw_representative_state());
        std::vector<std::shared_ptr < TranslationUniqueElement>> result;
        result.reserve(raw_coupled_states.size());
        for (const RawState & raw_coupled_state : raw_coupled_states) {
            result.push_back(stars_ring_basis::make_translation_unique_element(
                    element->physical_system(),
                    raw_coupled_state,
                    element->perturbation_order() + 1));
        }
        return result;
    }

}