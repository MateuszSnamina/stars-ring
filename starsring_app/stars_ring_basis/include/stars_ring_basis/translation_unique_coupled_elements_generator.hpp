#ifndef TRANSLATION_UNIQUE_COUPLED_STATES_GENERATOR_HPP
#define TRANSLATION_UNIQUE_COUPLED_STATES_GENERATOR_HPP

#include<stars_ring_basis/perturbation_basis_maker.hpp>
#include<stars_ring_basis/raw_state_coupled_elements_generator.hpp>

namespace stars_ring_basis {

    class CoupledTranslationUniqueElementsGenerator : public CoupledElementsGenerator<TranslationUniqueElement> {
    public:
        explicit CoupledTranslationUniqueElementsGenerator(
                std::shared_ptr<CoupledRawStatesGenerator> coupled_raw_states_generator);
        std::vector<std::shared_ptr<TranslationUniqueElement>> generate(
                std::shared_ptr<TranslationUniqueElement>) const;
    private:
        std::shared_ptr<CoupledRawStatesGenerator> _coupled_raw_states_generator;
    };

}
#endif