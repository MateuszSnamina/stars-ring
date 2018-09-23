#include<cassert>

#include<stars_ring_basis/localized_coupled_elements_generator.hpp>
#include<stars_ring_basis/raw_state_operations.hpp>

namespace stars_ring_basis {

    CoupledLocalizedElementsGenerator::CoupledLocalizedElementsGenerator(
            std::shared_ptr<CoupledRawStatesGenerator> coupled_raw_states_generator) :
    _coupled_raw_states_generator(coupled_raw_states_generator) {
        assert(coupled_raw_states_generator);
    }

    std::vector<std::shared_ptr<LocalizedElement>>
    CoupledLocalizedElementsGenerator::generate(std::shared_ptr<LocalizedElement> element) const {
        assert(element);
        assert(element->physical_system());
        assert(_coupled_raw_states_generator);
        const std::vector<RawState> raw_coupled_states =
                _coupled_raw_states_generator->generate(element->raw_state());
        std::vector<std::shared_ptr < LocalizedElement>> result;
        result.reserve(raw_coupled_states.size());
        for (const RawState & raw_coupled_state : raw_coupled_states) {
            result.push_back(stars_ring_basis::make_localized_element(
                    element->physical_system(),
                    raw_coupled_state,
                    element->perturbation_order() + 1));
        }
        return result;
    }

}