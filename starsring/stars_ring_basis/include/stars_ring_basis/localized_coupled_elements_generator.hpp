#ifndef LOCALIZED_COUPLED_STATES_GENERATOR_HPP
#define LOCALIZED_COUPLED_STATES_GENERATOR_HPP

#include<stars_ring_basis/perturbation_basis_maker.hpp>
#include<stars_ring_basis/raw_state_coupled_elements_generator.hpp>

namespace stars_ring_basis {

    class CoupledLocalizedElementsGenerator : public CoupledElementsGenerator<LocalizedElement> {
    public:
        explicit CoupledLocalizedElementsGenerator(
                std::shared_ptr<CoupledRawStatesGenerator> coupled_raw_states_generator);
        std::vector<std::shared_ptr<LocalizedElement>> generate(
                std::shared_ptr<LocalizedElement>) const;
    private:
        std::shared_ptr<CoupledRawStatesGenerator> _coupled_raw_states_generator;
    };

}

#endif